#!/usr/bin/env python3
"""
V13 GT 专用预测器 — 基于 V12 消融结果的针对性优化
========================================================

V12 → V13 核心变更 (仅 GT):

┌─────────────────── V12 GT 的问题 ─────────────────────────┐
│  1. 9/10 辅助通道有害 (ΔER 均为负):                       │
│     最有害: v_r_mean (-0.205), speed_mean (-0.150)        │
│     唯一有益: omega_mean (+0.064)                          │
│  2. 根本原因: 统计阶次不匹配                              │
│     GT = ⟨(v-⟨v⟩)²⟩ 是二阶中心矩                        │
│     辅助通道 (speed_mean, v_r_mean...) 是一阶量            │
│     一阶量→二阶量的跨尺度映射在 F4→F2 和 F2→F1 间不一致  │
│  3. 深层 FiLM 调制反而有害 (ΔER_wo_film_concat=-0.147)   │
│  Baseline ER = 0.615, only_torque ER = 0.479              │
│  wo_film_input_concat ER = 0.468                          │
└──────────────────────────────────────────────────────────┘

┌─────────────────── V13 的解决 ─────────────────────────────┐
│  1. 辅助通道精选 (物理动机):                               │
│     - 移除全部一阶有害通道                                 │
│     - 保留 torque + omega_mean                             │
│     - 新增二阶特征 (v_r_std, v_theta_std 等) — 可选       │
│  2. 架构微调 (消融验证):                                   │
│     - 用 input_concat 替代深层 FiLM (ΔER=-0.147)          │
│     - 移除二阶差分 (对 GT 有害, ΔER=-0.053)               │
│     - 保留频率分解和多尺度局部特征 (有益)                  │
│  3. 架构本体不变: 仍使用 PredictorV12 类                   │
│  4. 支持两种模式:                                          │
│     - minimal: torque + omega_mean (3 ch, 无需新数据)      │
│     - enhanced: +4个二阶特征 (7 ch, 需要新comparison文件)  │
└──────────────────────────────────────────────────────────┘

论文叙述:
  "Based on ablation analysis revealing statistical order mismatch as the
   root cause of auxiliary channel harm, we design a GT-specific input
   configuration that replaces first-order kinematic auxiliaries with
   second-order fluctuation features matching GT's variance nature.
   Architecture settings are optimized per ablation: shallow condition
   encoding (input concatenation) replaces deep FiLM modulation, while
   the model class remains identical to the unified framework."
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional


# =============================================================================
# GT V13 辅助通道定义
# =============================================================================
"""
V13 GT 辅助通道设计原则:
  - 移除所有一阶运动学量 (speed_mean, v_r_mean, v_theta_mean)
  - 移除所有积分量 (KE_total, KE_rot_sum)
  - 移除冗余自身量 (T_total 作为辅助)
  - 保留 torque (驱动条件) + omega_mean (旋转校正)
  - 新增二阶涨落特征 (与 GT 统计阶次匹配)

Minimal 模式 (2 aux channels):
  [0] torque
  [1] omega_mean

Enhanced 模式 (6 aux channels):
  [0] torque
  [1] omega_mean
  [2] v_r_std        — 径向速度涨落 (二阶)
  [3] v_theta_std    — 切向速度涨落 (二阶)
  [4] v_z_cyl_std    — 轴向速度涨落 (二阶)
  [5] speed_std      — 速率涨落 (二阶)
"""

# =============================================================================
# GT V13 辅助通道定义 [V3.2 MODIFIED]
# =============================================================================
# Minimal 模式辅助通道
AUX_CHANNELS_MINIMAL = [
    'torque',       # [0] 驱动条件基准
    'omega_mean',   # [1] 整体旋转角速度 (GT_theta 校正所需)
]

# Enhanced 模式辅助通道 (移除有缩放风险的动能、引入各向异性)
AUX_CHANNELS_ENHANCED = [
    'torque',           # [0] 驱动条件基准
    'omega_mean',       # [1] 整体旋转角速度
    'v_r_std',          # [2] 径向速度涨落强度 (质量加权二阶)
    'v_theta_std',      # [3] 切向速度涨落强度 (质量加权二阶)
    'v_z_cyl_std',      # [4] 轴向速度涨落强度 (质量加权二阶)
    'GT_aniso',         # [5] 颗粒温度各向异性比 (跨尺度衰减比例的重要指导)
]

# 辅助通道数据源映射 (文件路径)
AUX_SOURCES_GT = {
    'torque':           {'type': 'sta', 'file': 'mixer_torque_comparison.csv'},
    'omega_mean':       {'type': 'ext', 'file': 'comparison_omega_mean.csv'},
    'v_r_std':          {'type': 'ext', 'file': 'comparison_v_r_std.csv'},
    'v_theta_std':      {'type': 'ext', 'file': 'comparison_v_theta_std.csv'},
    'v_z_cyl_std':      {'type': 'ext', 'file': 'comparison_v_z_cyl_std.csv'},
    'GT_aniso':         {'type': 'ext', 'file': 'comparison_GT_aniso.csv'},
}


def get_aux_channels(mode: str = 'auto') -> List[str]:
    """
    获取 GT V13 辅助通道列表.
    
    Args:
        mode: 'minimal', 'enhanced', 或 'auto' (自动检测文件)
    """
    if mode == 'minimal':
        return list(AUX_CHANNELS_MINIMAL)
    elif mode == 'enhanced':
        return list(AUX_CHANNELS_ENHANCED)
    else:
        return list(AUX_CHANNELS_MINIMAL)  # auto 默认 minimal, 运行时检测


# =============================================================================
# GT V13 配置
# =============================================================================
@dataclass
class GTConfigV13:
    """
    GT V13 专用配置.
    
    架构与 V12 UnifiedConfig 保持一致, 但针对 GT 消融结果优化:
      - 辅助通道精选 (minimal/enhanced)
      - use_film=False + input_concat (浅层条件编码)
      - use_2nd_order_diff=False (GT 有害)
      - weight_transient=0 (GT 瞬态段无意义)
    """
    
    # --- GT 专有: 辅助通道模式 ---
    aux_mode: str = 'auto'  # 'minimal', 'enhanced', 'auto'
    
    # --- 物理必要差异 ---
    n_decoder_heads: int = 1
    head_names: Tuple[str, ...] = ('gt',)
    correction_mode: str = 'multiplicative'
    max_ratios: Tuple[float, ...] = (0.8,)  # GT 需要更大的修正范围
    head_weights: Tuple[float, ...] = (1.0,)
    
    # --- 输入通道 (由 aux_mode 运行时确定) ---
    n_aux_channels: int = 2    # 运行时更新
    n_target_channels: int = 1  # GT 始终 1 个目标头
    
    @property
    def input_channels(self) -> int:
        return self.n_aux_channels + self.n_target_channels
    
    @property
    def target_ch_indices(self) -> Tuple[int, ...]:
        return (self.n_aux_channels,)
    
    # --- 条件编码 (与 V12 一致) ---
    condition_base_dim: int = 14
    n_rel_amp_inputs: int = 1
    condition_dim: int = 64
    feature_dim: int = 64
    
    # --- 架构: V12 消融优化设置 ---
    n_film_layers: int = 4
    
    # GT 消融结论: FiLM 有害, input_concat 有益
    use_film: bool = False
    film_ablation_mode: str = 'input_concat'
    
    # GT 消融结论: 频率分解有益 (+0.085), 局部特征有益 (+0.082)
    use_freq_decomp: bool = True
    use_local_features: bool = True
    
    # GT 消融结论: 二阶差分有害 (-0.053)
    use_2nd_order_diff: bool = False
    
    # 保持常开的组件 (与 V12 一致)
    use_direct_decoder: bool = True
    direct_res_scale: float = 1.0
    use_channel_attention: bool = True
    use_high_coh_branch: bool = True
    
    # --- 训练 (与 V12 一致) ---
    batch_size: int = 8
    learning_rate: float = 1e-3
    n_epochs: int = 600
    weight_decay: float = 1e-4
    
    # --- 损失权重 ---
    weight_cycle: float = 0.3
    weight_shape: float = 0.2
    weight_steady: float = 1.0
    weight_transient: float = 0.0  # GT: 瞬态段无意义
    weight_ratio_reg: float = 0.05
    weight_direct: float = 0.3
    weight_cohesion_adapt: float = 0.2
    
    # --- 频率分解 ---
    rel_amp_threshold: float = 0.025
    osc_decay_factor: float = 15.0
    osc_decay_min: float = 0.3
    cutoff_freq: float = 3.0
    filter_order: int = 4
    
    # --- 时间 ---
    n_time_points: int = 200
    time_range: Tuple[float, float] = (0.0, 3.0)
    steady_start: float = 0.6
    eval_start_time: float = 0.6
    
    # --- 输入消融 (V13 中一般不使用, 但保持兼容) ---
    channel_mask: Optional[Tuple[bool, ...]] = None
    
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    
    def set_aux_channels(self, channels: List[str]):
        """设置辅助通道并更新相关参数"""
        self.n_aux_channels = len(channels)
        self._aux_channel_names = channels
    
    def get_aux_channel_names(self) -> List[str]:
        if hasattr(self, '_aux_channel_names'):
            return self._aux_channel_names
        return get_aux_channels(self.aux_mode)
    
    def validate(self):
        assert len(self.max_ratios) == self.n_decoder_heads
        assert len(self.head_names) == self.n_decoder_heads
        assert len(self.head_weights) == self.n_decoder_heads
        assert self.correction_mode in ('multiplicative', 'additive')
        assert self.aux_mode in ('minimal', 'enhanced', 'auto')
        if self.channel_mask is not None:
            assert len(self.channel_mask) == self.n_aux_channels
    
    def get_active_aux_channels(self) -> List[str]:
        names = self.get_aux_channel_names()
        if self.channel_mask is None:
            return names
        return [n for n, m in zip(names, self.channel_mask) if m]
    
    def get_masked_aux_channels(self) -> List[str]:
        names = self.get_aux_channel_names()
        if self.channel_mask is None:
            return []
        return [n for n, m in zip(names, self.channel_mask) if not m]
    
    def get_arch_ablation_desc(self) -> List[str]:
        ablations = []
        if not self.use_film:
            if self.film_ablation_mode == 'input_concat':
                ablations.append('input_concat (GT-optimized)')
            else:
                ablations.append('w/o FiLM')
        if not self.use_freq_decomp:
            ablations.append('w/o freq decomp')
        if not self.use_local_features:
            ablations.append('w/o local features')
        elif not self.use_2nd_order_diff:
            ablations.append('w/o 2nd-order diff (GT-optimized)')
        return ablations


# =============================================================================
# 网络组件 (与 V12 完全一致, 复制以保持独立性)
# =============================================================================
def apply_directional_decay(ratio: torch.Tensor, decay_scale: torch.Tensor) -> torch.Tensor:
    """方向性衰减: 正ratio衰减, 负ratio保持"""
    positive_mask = (ratio > 0).float()
    if decay_scale.dim() == 0:
        decayed = ratio * decay_scale
    elif decay_scale.dim() == 1 and ratio.dim() == 2:
        decayed = ratio * decay_scale.unsqueeze(1)
    else:
        decayed = ratio * decay_scale
    return positive_mask * decayed + (1 - positive_mask) * ratio


class ConditionEncoder(nn.Module):
    """条件编码器 (含 High-Coh Branch, 与 V12 一致)"""
    def __init__(self, base_dim: int, n_rel_amp: int, hid: int, out: int):
        super().__init__()
        total_in = base_dim + n_rel_amp + 1
        self.net = nn.Sequential(
            nn.Linear(total_in, hid), nn.LayerNorm(hid), nn.GELU(),
            nn.Linear(hid, hid), nn.LayerNorm(hid), nn.GELU(),
            nn.Linear(hid, out))
        self.high_coh_branch = nn.Sequential(
            nn.Linear(total_in, hid // 2), nn.GELU(),
            nn.Linear(hid // 2, out))
    
    def forward(self, base, *extras):
        x = torch.cat([base] + list(extras), dim=1)
        out = self.net(x)
        gate = base[:, -1:].clamp(0, 1)
        out = out + gate * self.high_coh_branch(x)
        return out


class LocalFeatureExtractor(nn.Module):
    """多尺度局部特征 (与 V12 一致)"""
    def __init__(self, in_ch=1, out_ch=16, use_2nd_order_diff=True):
        super().__init__()
        self.use_2nd_order_diff = use_2nd_order_diff
        if use_2nd_order_diff:
            q = out_ch // 4
            self.conv_s = nn.Conv1d(in_ch, q, 3, padding=1)
            self.conv_m = nn.Conv1d(in_ch, q, 7, padding=3)
            self.conv_l = nn.Conv1d(in_ch, q, 15, padding=7)
            self.conv_d = nn.Conv1d(in_ch, q, 3, padding=1)
            self.register_buffer('diff_kernel', torch.tensor([[[1.0, -2.0, 1.0]]]))
        else:
            q1 = out_ch // 3
            q2 = out_ch // 3
            q3 = out_ch - q1 - q2
            self.conv_s = nn.Conv1d(in_ch, q1, 3, padding=1)
            self.conv_m = nn.Conv1d(in_ch, q2, 7, padding=3)
            self.conv_l = nn.Conv1d(in_ch, q3, 15, padding=7)
    
    def forward(self, x):
        if self.use_2nd_order_diff:
            diff = F.conv1d(x, self.diff_kernel, padding=1)
            return torch.cat([F.gelu(self.conv_s(x)), F.gelu(self.conv_m(x)),
                              F.gelu(self.conv_l(x)), F.gelu(self.conv_d(diff))], dim=1)
        else:
            return torch.cat([F.gelu(self.conv_s(x)), F.gelu(self.conv_m(x)),
                              F.gelu(self.conv_l(x))], dim=1)


class CurveEncoder(nn.Module):
    """曲线编码器 (与 V12 一致)"""
    def __init__(self, in_ch: int, out_ch: int = 64,
                 target_ch_indices: Tuple[int, ...] = (2,),
                 n_local_per_target: int = 16,
                 use_local_features: bool = True,
                 use_2nd_order_diff: bool = True):
        super().__init__()
        self.use_local_features = use_local_features
        self.conv = nn.Sequential(
            nn.Conv1d(in_ch, 32, 7, padding=3), nn.GroupNorm(8, 32), nn.GELU(),
            nn.Conv1d(32, 48, 5, padding=2), nn.GroupNorm(8, 48), nn.GELU(),
            nn.Conv1d(48, out_ch, 3, padding=1))
        
        if use_local_features:
            self.target_ch_indices = target_ch_indices
            self.locals = nn.ModuleList([
                LocalFeatureExtractor(1, n_local_per_target,
                                       use_2nd_order_diff=use_2nd_order_diff)
                for _ in target_ch_indices
            ])
            total_local = n_local_per_target * len(target_ch_indices)
            self.fusion = nn.Conv1d(out_ch + total_local, out_ch, 1)
    
    def forward(self, x):
        main = self.conv(x)
        if self.use_local_features:
            local_feats = [self.locals[i](x[:, idx:idx+1, :])
                           for i, idx in enumerate(self.target_ch_indices)]
            return self.fusion(torch.cat([main] + local_feats, dim=1))
        else:
            return main


class FiLMLayer(nn.Module):
    def __init__(self, feat: int, cond: int):
        super().__init__()
        self.gamma = nn.Linear(cond, feat)
        self.beta = nn.Linear(cond, feat)
        nn.init.ones_(self.gamma.bias); nn.init.zeros_(self.gamma.weight)
        nn.init.zeros_(self.beta.weight); nn.init.zeros_(self.beta.bias)
    
    def forward(self, f, c):
        return self.gamma(c).unsqueeze(-1) * f + self.beta(c).unsqueeze(-1)


class FiLMResBlock(nn.Module):
    """FiLM 残差块 (含 Channel Attention, 与 V12 一致)"""
    def __init__(self, ch: int, cond: int, use_film: bool = True):
        super().__init__()
        self.use_film = use_film
        self.conv1 = nn.Conv1d(ch, ch, 3, padding=1)
        self.conv2 = nn.Conv1d(ch, ch, 3, padding=1)
        if use_film:
            self.film = FiLMLayer(ch, cond)
        self.norm1 = nn.GroupNorm(8, ch)
        self.norm2 = nn.GroupNorm(8, ch)
        self.ca = nn.Sequential(
            nn.AdaptiveAvgPool1d(1), nn.Flatten(),
            nn.Linear(ch, ch // 4), nn.GELU(),
            nn.Linear(ch // 4, ch), nn.Sigmoid())
    
    def forward(self, x, c):
        r = x
        x = F.gelu(self.norm1(x)); x = self.conv1(x)
        x = F.gelu(self.norm2(x)); x = self.conv2(x)
        if self.use_film:
            x = self.film(x, c)
        x = x * self.ca(x).unsqueeze(-1)
        return x + r


class RatioDecoder(nn.Module):
    def __init__(self, in_ch: int, max_ratio: float):
        super().__init__()
        self.max_ratio = max_ratio
        self.conv = nn.Sequential(
            nn.Conv1d(in_ch, 32, 5, padding=2), nn.GroupNorm(8, 32), nn.GELU(),
            nn.Conv1d(32, 16, 3, padding=1), nn.GELU(),
            nn.Conv1d(16, 1, 3, padding=1))
    
    def forward(self, x):
        return self.max_ratio * torch.tanh(self.conv(x).squeeze(1))


class DirectDecoder(nn.Module):
    def __init__(self, in_ch: int):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(in_ch, 32, 5, padding=2), nn.GroupNorm(8, 32), nn.GELU(),
            nn.Conv1d(32, 16, 3, padding=1), nn.GELU(),
            nn.Conv1d(16, 1, 3, padding=1))
    
    def forward(self, x):
        return self.conv(x).squeeze(1)


# =============================================================================
# GT V13 预测器 (架构与 PredictorV12 完全一致)
# =============================================================================
class GTPredictorV13(nn.Module):
    """
    GT V13 专用预测器.
    
    架构与 PredictorV12 完全一致:
      [ChannelMask] → [CurveEncoder] → [FiLM/NonFiLM+CA Blocks] → features
          ├→ [RatioDecoder]     + [RatioDecoder_low]*
          └→ [DirectDecoder]
    
    区别仅在配置层面:
      - 辅助通道精选 (2 或 6 个, 而非 V12 的 8 个)
      - use_film=False, input_concat 模式
      - use_2nd_order_diff=False
    """
    
    def __init__(self, config: GTConfigV13):
        super().__init__()
        config.validate()
        self.config = config
        
        # --- Channel mask ---
        if config.channel_mask is not None:
            full_mask = list(config.channel_mask) + [True] * config.n_target_channels
            self.register_buffer('channel_mask',
                                 torch.tensor(full_mask, dtype=torch.bool))
        else:
            self.register_buffer('channel_mask', None)
        
        # --- 条件编码器 ---
        self.cond_enc = ConditionEncoder(
            base_dim=config.condition_base_dim,
            n_rel_amp=config.n_rel_amp_inputs,
            hid=config.condition_dim,
            out=config.condition_dim,
        )
        
        # --- Input concat 模式 ---
        self._use_input_concat = (not config.use_film and
                                  config.film_ablation_mode == 'input_concat')
        if self._use_input_concat:
            curve_in_ch = config.input_channels + config.condition_dim
        else:
            curve_in_ch = config.input_channels
        
        # --- 曲线编码器 ---
        self.curve_enc = CurveEncoder(
            in_ch=curve_in_ch,
            out_ch=config.feature_dim,
            target_ch_indices=config.target_ch_indices,
            use_local_features=config.use_local_features,
            use_2nd_order_diff=config.use_2nd_order_diff,
        )
        
        # --- FiLM/非FiLM 残差块 ---
        self.blocks = nn.ModuleList([
            FiLMResBlock(config.feature_dim, config.condition_dim,
                         use_film=config.use_film)
            for _ in range(config.n_film_layers)
        ])
        
        # --- Ratio 解码器 ---
        self.ratio_decoders = nn.ModuleList([
            RatioDecoder(config.feature_dim, config.max_ratios[0])
        ])
        
        # --- 低频 Ratio 解码器 ---
        if config.use_freq_decomp:
            self.ratio_decoders_low = nn.ModuleList([
                RatioDecoder(config.feature_dim, config.max_ratios[0])
            ])
        else:
            self.ratio_decoders_low = None
        
        # --- 直接残差解码器 ---
        self.direct_decoder = DirectDecoder(config.feature_dim)
    
    def apply_channel_mask(self, x: torch.Tensor) -> torch.Tensor:
        if self.channel_mask is not None:
            mask = self.channel_mask.view(1, -1, 1)
            return x * mask.float()
        return x
    
    def encode(self, curves, cond_base, *rel_amps_and_scale):
        curves = self.apply_channel_mask(curves)
        cond = self.cond_enc(cond_base, *rel_amps_and_scale)
        
        if self._use_input_concat:
            B, C, T = curves.shape
            cond_broadcast = cond.unsqueeze(-1).expand(B, -1, T)
            curves = torch.cat([curves, cond_broadcast], dim=1)
        
        f = self.curve_enc(curves)
        for blk in self.blocks:
            f = blk(f, cond)
        return f, cond
    
    def decode_ratios(self, features):
        ratio = self.ratio_decoders[0](features)
        ratio_low = (self.ratio_decoders_low[0](features)
                     if self.ratio_decoders_low is not None else None)
        return ratio, ratio_low
    
    def decode_direct(self, features) -> torch.Tensor:
        return self.direct_decoder(features)
    
    def apply_correction(self, ratio, ratio_low,
                         use_freq, osc_decay,
                         raw, low, high,
                         direct_res=None):
        B = raw.shape[0]
        pred = torch.zeros_like(raw)
        eff_ratio_low = torch.zeros_like(ratio if ratio_low is None else ratio_low)
        
        freq_enabled = self.config.use_freq_decomp
        
        if freq_enabled:
            uf = use_freq.squeeze(-1) if use_freq.dim() > 1 else use_freq
            ds = osc_decay.squeeze(-1) if osc_decay.dim() > 1 else osc_decay
            
            for i in range(B):
                if uf[i]:
                    eff = apply_directional_decay(ratio_low[i:i+1], ds[i:i+1]).squeeze(0)
                    eff_ratio_low[i] = eff
                    pred[i] = low[i] * (1 + eff) + high[i]
                else:
                    pred[i] = raw[i] * (1 + ratio[i])
                    eff_ratio_low[i] = ratio_low[i] if ratio_low is not None else ratio[i]
        else:
            pred = raw * (1 + ratio)
            eff_ratio_low = ratio
        
        if direct_res is not None:
            pred = pred + direct_res * self.config.direct_res_scale
        
        return pred, eff_ratio_low
    
    def print_architecture(self):
        c = self.config
        n_params = sum(p.numel() for p in self.parameters())
        aux_names = c.get_aux_channel_names()
        masked = c.get_masked_aux_channels()
        mask_str = f"MASKED: {masked}" if masked else "ALL ACTIVE"
        arch_abl = c.get_arch_ablation_desc()
        arch_str = ', '.join(arch_abl) if arch_abl else "NONE"
        
        print(f"\n  ┌─── GTPredictorV13 Architecture ──────────────────┐")
        print(f"  │  Version:      V13 (GT-optimized)")
        print(f"  │  Aux mode:     {c.aux_mode} ({len(aux_names)} channels)")
        print(f"  │  Aux channels: {aux_names}")
        print(f"  │  Correction:   {c.correction_mode}")
        print(f"  │  max_ratio:    {c.max_ratios}")
        print(f"  │  Input ch:     {c.input_channels} (aux={c.n_aux_channels} + target={c.n_target_channels})")
        if self._use_input_concat:
            print(f"  │  Input concat: +{c.condition_dim} cond → {c.input_channels + c.condition_dim} total")
        print(f"  │  FiLM:         {'ON' if c.use_film else 'OFF → input_concat'}")
        print(f"  │  Freq decomp:  {'ON' if c.use_freq_decomp else 'OFF'}")
        print(f"  │  Local feats:  {'ON' if c.use_local_features else 'OFF'}"
              f"{'  (no 2nd-order diff)' if c.use_local_features and not c.use_2nd_order_diff else ''}")
        print(f"  │  GT-specific:  wo_2nd_order_diff, input_concat, curated_aux")
        print(f"  │  ─────────────────────────────────────────────── │")
        print(f"  │  Arch opts:    {arch_str}")
        print(f"  │  Channel mask: {mask_str}")
        print(f"  │  Total params: {n_params:,}")
        print(f"  └─────────────────────────────────────────────────┘")


# =============================================================================
# V13 消融实验工厂 (GT 专用)
# =============================================================================
def list_gt_v13_ablation_experiments(aux_channels: List[str]) -> Dict[str, Dict]:
    """
    V13 GT 专用消融实验.
    
    与 V12 消融不同, V13 消融用于验证:
    1. minimal vs enhanced 模式的差异
    2. 各二阶特征的个体贡献
    3. omega_mean 在 V13 中是否仍然有益
    """
    experiments = {}
    
    # Baseline: 全部保留
    experiments['baseline'] = {
        'remove': [],
        'mask': tuple([True] * len(aux_channels)),
        'hypothesis': 'V13 full model baseline',
    }
    
    # 逐一移除每个辅助通道
    for i, ch in enumerate(aux_channels):
        mask = [True] * len(aux_channels)
        mask[i] = False
        experiments[f'remove_{ch}'] = {
            'remove': [ch],
            'mask': tuple(mask),
            'hypothesis': f'Test contribution of {ch}',
        }
    
    # 仅保留 torque (最小输入)
    if len(aux_channels) > 1:
        mask = [False] * len(aux_channels)
        mask[0] = True  # torque
        experiments['only_torque'] = {
            'remove': aux_channels[1:],
            'mask': tuple(mask),
            'hypothesis': 'Torque-only baseline (minimal information)',
        }
    
    # Enhanced 模式: 移除所有二阶特征 [V3.2 MODIFIED]
    second_order = ['v_r_std', 'v_theta_std', 'v_z_cyl_std', 'GT_aniso']
    so_present = [ch for ch in second_order if ch in aux_channels]
    if so_present:
        mask = [True] * len(aux_channels)
        for ch in so_present:
            mask[aux_channels.index(ch)] = False
        experiments['remove_all_2nd_order'] = {
            'remove': so_present,
            'mask': tuple(mask),
            'hypothesis': 'Remove all 2nd-order features, test if they help GT',
        }
    return experiments
