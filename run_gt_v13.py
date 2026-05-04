#!/usr/bin/env python3
"""
V13 GT 预测运行脚本 — 基于 V12 消融结果的 GT 专用优化
========================================================

用法:
  # 自动模式 (检测二阶特征文件, 自动选择 minimal/enhanced)
  python run_gt_v13.py
  
  # 强制 minimal 模式 (仅 torque + omega_mean)
  python run_gt_v13.py --mode minimal
  
  # 强制 enhanced 模式 (含二阶涨落特征)
  python run_gt_v13.py --mode enhanced
  
  # 与 V12 baseline 对比运行
  python run_gt_v13.py --compare_v12
  
  # V13 内部消融实验
  python run_gt_v13.py --ablation_sweep

V12 → V13 数据管线变更:
  - 辅助通道从 8 个精简为 2-6 个 (基于消融结果)
  - 新增二阶涨落特征 (v_r_std, v_theta_std 等) 的加载支持
  - 条件编码不变 (14-dim)
  - 其余数据管线与 V12 一致
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional, List, Any
import torch
import torch.nn as nn
import torch.nn.functional as F_torch
from torch.utils.data import Dataset, DataLoader
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt
import pickle
import argparse
import json
import random
import warnings
warnings.filterwarnings('ignore')


def set_seed(seed: int = 42):
    """固定全局随机种子"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


from gt_predictor_v13 import (
    GTPredictorV13, GTConfigV13,
    apply_directional_decay,
    AUX_CHANNELS_MINIMAL, AUX_CHANNELS_ENHANCED, AUX_SOURCES_GT,
    get_aux_channels, list_gt_v13_ablation_experiments,
)


# =============================================================================
# [1] 配置与常量
# =============================================================================
MATERIALS = ('standard', 'DualDensity', 'SizeDiff')
COHESIONS = ('000', '005', '010', '020')
COHESION_VALUES = {'000': 0.00, '005': 0.005, '010': 0.01, '020': 0.02}
COHESION_WEIGHTS = {'000': 1.0, '005': 1.0, '010': 1.2, '020': 1.5}


# =============================================================================
# [2] 通用工具
# =============================================================================
def safe_normalize(arr):
    return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)


def compute_osc_decay_scale(rel_amp, decay_factor=15.0, min_scale=0.3):
    return max(min_scale, 1.0 - rel_amp * decay_factor)


def compute_relative_amplitude(signal, steady_mask=None):
    sig = signal[steady_mask] if steady_mask is not None else signal
    sig = sig[~np.isnan(sig)] if hasattr(sig, '__len__') else sig
    if len(sig) == 0:
        return 0.0
    return float(np.std(sig) / (np.mean(np.abs(sig)) + 1e-10))


class FrequencyDecomposer:
    def __init__(self, cutoff_freq, sample_rate, order=4):
        nyquist = sample_rate / 2
        if cutoff_freq >= nyquist:
            cutoff_freq = nyquist * 0.9
        self.b_low, self.a_low = butter(order, cutoff_freq / nyquist, btype='low')
        self.order = order

    def decompose(self, signal):
        signal = safe_normalize(signal)
        try:
            padlen = min(len(signal) - 1, 3 * self.order)
            low = filtfilt(self.b_low, self.a_low, signal, padlen=padlen)
        except:
            low = np.convolve(signal, np.ones(15)/15, mode='same')
        return low.astype(np.float32), (signal - low).astype(np.float32)


class AdaptiveCohesionEncoder:
    LOG_SCALE = 100

    def __init__(self, cohesion_values):
        sorted_vals = sorted(set(cohesion_values.values()))
        nonzero = [v for v in sorted_vals if v > 0]
        self.max_coh = max(sorted_vals) if max(sorted_vals) > 0 else 1.0
        self.log_norm = np.log1p(self.max_coh * self.LOG_SCALE)
        self.low_threshold, self.low_scale = self._sp(nonzero, 'low')
        self.high_threshold, self.high_scale = self._sp(nonzero, 'high')

    @staticmethod
    def _sp(nz, mode):
        if len(nz) < 2:
            val = nz[0] if nz else 1.0
            return val / 2, 10.0 / val
        a, b = (nz[0], nz[1]) if mode == 'low' else (nz[-2], nz[-1])
        return np.sqrt(a * b), 2 * np.log(7) / (b - a)

    def encode(self, cv):
        return np.array([
            cv / self.max_coh,
            np.log1p(cv * self.LOG_SCALE) / self.log_norm,
            np.sqrt(cv / self.max_coh),
            1 / (1 + np.exp(-self.low_scale * (cv - self.low_threshold))),
            1 / (1 + np.exp(-self.high_scale * (cv - self.high_threshold))),
        ], dtype=np.float32)

    @property
    def dim(self):
        return 5


class MetricsCalculator:
    def __init__(self, eval_mask=None):
        self.eval_mask = eval_mask

    def compute(self, pred, target, baseline):
        pred, target, baseline = safe_normalize(pred), safe_normalize(target), safe_normalize(baseline)
        if self.eval_mask is not None:
            pred, target, baseline = pred[self.eval_mask], target[self.eval_mask], baseline[self.eval_mask]
        pr = np.sqrt(np.mean((pred - target)**2))
        br = np.sqrt(np.mean((baseline - target)**2))
        tm = np.mean(np.abs(target)) + 1e-10
        er = pr / (br + 1e-10)
        corr = np.corrcoef(pred, target)[0, 1] if len(pred) > 1 else 0
        corr = 0.0 if np.isnan(corr) else corr
        pred_rel = pr / tm
        base_rel = br / tm
        small = base_rel < 0.03
        success = er < 1.0 or pred_rel < 0.05 or (small and pred_rel < base_rel + 0.01)
        return {
            'rmse': pr, 'baseline_rmse': br, 'error_ratio': er,
            'pred_rel': pred_rel * 100, 'baseline_rel': base_rel * 100,
            'correlation': corr, 'success': success,
        }


# =============================================================================
# [3] 二阶特征文件检测 [V3.2 MODIFIED]
# =============================================================================
def detect_second_order_features(ext_dir: str) -> List[str]:
    """
    检测二阶涨落特征 comparison 文件是否存在.
    """
    ext_path = Path(ext_dir)
    # 仅检测真正量纲安全且有物理意义的特征
    second_order_features = [
        'v_r_std', 'v_theta_std', 'v_z_cyl_std', 'GT_aniso',
    ]
    available = []
    
    for feat in second_order_features:
        file_name = f"comparison_{feat}.csv"
        all_exist = True
        some_exist = False
        
        for mat in MATERIALS:
            for coh in COHESIONS:
                config_dir = ext_path / f"{mat}_se{coh}"
                if config_dir.exists():
                    if (config_dir / file_name).exists():
                        some_exist = True
                    else:
                        all_exist = False
        
        if all_exist and some_exist:
            available.append(feat)
    
    return available


def auto_detect_mode(ext_dir: str) -> Tuple[str, List[str]]:
    """
    自动检测辅助通道模式.
    
    Returns:
        (mode_name, channel_list)
    """
    available_2nd = detect_second_order_features(ext_dir)
    
    if len(available_2nd) >= 3:  # 至少 3 个二阶特征可用
        channels = list(AUX_CHANNELS_MINIMAL) + available_2nd
        mode = 'enhanced'
        print(f"\n[Auto-detect] Found {len(available_2nd)} second-order features: {available_2nd}")
        print(f"[Auto-detect] Mode: ENHANCED ({len(channels)} aux channels)")
    else:
        channels = list(AUX_CHANNELS_MINIMAL)
        mode = 'minimal'
        if available_2nd:
            print(f"\n[Auto-detect] Only {len(available_2nd)} second-order features found: {available_2nd}")
            print(f"[Auto-detect] Need ≥3 for enhanced mode. Using MINIMAL ({len(channels)} aux channels)")
        else:
            print(f"\n[Auto-detect] No second-order features found")
            print(f"[Auto-detect] Mode: MINIMAL ({len(channels)} aux channels)")
            print(f"[Auto-detect] To enable enhanced mode, run sta_other_info_v2.py with")
            print(f"              expanded key_features including v_r_std, v_theta_std, etc.")
    
    return mode, channels


# =============================================================================
# [4] GT V13 数据加载
# =============================================================================
def load_gt_data(aux_channels: List[str],
                 sta_dir: str = 'sta_results',
                 ext_dir: str = 'extended_features_results') -> Dict:
    """
    加载 GT V13 数据.
    
    与 V12 的区别: 仅加载 GT 需要的辅助通道 (2-6 个而非 8 个)
    """
    sta_path, ext_path = Path(sta_dir), Path(ext_dir)
    all_data = {}
    
    print(f"\n[Loading GT V13 data]")
    print(f"  Aux channels: {aux_channels}")
    
    for mat in MATERIALS:
        for coh in COHESIONS:
            key = f"{mat}_se{coh}"
            config_sta = sta_path / key
            config_ext = ext_path / key
            
            # --- 必须有 torque ---
            torque_file = config_sta / "mixer_torque_comparison.csv"
            if not torque_file.exists():
                continue
            try:
                torque_df = pd.read_csv(torque_file)
            except:
                continue
            
            entry = {
                'material': mat, 'cohesion': coh,
                'time_torque': torque_df['Time(s)'].values,
                'torque': {s: safe_normalize(torque_df[f'f{s}'].values) for s in [1, 2, 4]},
            }
            
            # --- GT 目标量 ---
            gt_file = config_ext / "comparison_T_total.csv"
            if not gt_file.exists():
                continue
            try:
                gt_df = pd.read_csv(gt_file)
                tc = 'Time(s)' if 'Time(s)' in gt_df.columns else 'time'
                entry['gt'] = {
                    s: {'values': gt_df[f'f{s}'].values, 'time': gt_df[tc].values}
                    for s in [1, 2, 4]
                }
            except:
                continue
            
            # --- 辅助通道 ---
            entry['aux_channels'] = {}
            for aux_name in aux_channels:
                src = AUX_SOURCES_GT.get(aux_name)
                if src is None:
                    # 未知辅助通道, 尝试从 ext 目录加载
                    comp_file = config_ext / f"comparison_{aux_name}.csv"
                    if comp_file.exists():
                        try:
                            df = pd.read_csv(comp_file)
                            tc_a = 'Time(s)' if 'Time(s)' in df.columns else 'time'
                            if all(f"f{s}" in df.columns for s in [1, 2, 4]):
                                entry['aux_channels'][aux_name] = {
                                    'data': {s: df[f"f{s}"].values for s in [1, 2, 4]},
                                    'time': {s: df[tc_a].values for s in [1, 2, 4]},
                                }
                        except:
                            pass
                elif src['type'] == 'sta':
                    if aux_name == 'torque':
                        entry['aux_channels'][aux_name] = {
                            'data': entry['torque'],
                            'time': {s: entry['time_torque'] for s in [1, 2, 4]},
                        }
                elif src['type'] == 'ext':
                    comp_file = config_ext / src['file']
                    if comp_file.exists():
                        try:
                            df = pd.read_csv(comp_file)
                            tc_a = 'Time(s)' if 'Time(s)' in df.columns else 'time'
                            if all(f"f{s}" in df.columns for s in [1, 2, 4]):
                                entry['aux_channels'][aux_name] = {
                                    'data': {s: df[f"f{s}"].values for s in [1, 2, 4]},
                                    'time': {s: df[tc_a].values for s in [1, 2, 4]},
                                }
                        except:
                            pass
                
                if aux_name not in entry['aux_channels']:
                    entry['aux_channels'][aux_name] = None
            
            all_data[key] = entry
            n_loaded = sum(1 for v in entry['aux_channels'].values() if v is not None)
            print(f"  ✓ {key} (aux: {n_loaded}/{len(aux_channels)})")
    
    print(f"\nLoaded {len(all_data)} configurations")
    return all_data


# =============================================================================
# [5] 预处理器
# =============================================================================
class GTPreprocessor:
    def __init__(self, config: GTConfigV13, aux_channels: List[str]):
        self.config = config
        self.aux_channels = aux_channels
        self.target_time = np.linspace(config.time_range[0], config.time_range[1],
                                       config.n_time_points)
        self.steady_mask = self.target_time >= config.steady_start
        dt = (config.time_range[1] - config.time_range[0]) / max(config.n_time_points - 1, 1)
        self.decomposer = FrequencyDecomposer(config.cutoff_freq, 1.0 / dt, config.filter_order)
        self.stats = {}

    def fit(self, all_data):
        collectors: Dict[str, list] = {}
        for data in all_data.values():
            for scale in [2, 4]:
                for aux_name in self.aux_channels:
                    aux = data['aux_channels'].get(aux_name)
                    if aux is not None and scale in aux['data']:
                        collectors.setdefault(f'aux_{aux_name}', []).append(aux['data'][scale])
                hdata = data['gt'][scale]
                collectors.setdefault('head_gt', []).append(hdata['values'])
        
        for key, arrs in collectors.items():
            v = np.concatenate(arrs)
            v = v[np.isfinite(v)]
            if len(v) > 0:
                self.stats[key] = (float(np.mean(v)), float(max(np.std(v), 1e-8)))
        
        print(f"\n[Normalization] (F4+F2 only, {self.config.n_time_points} pts)")
        if 'head_gt' in self.stats:
            print(f"  GT: mean={self.stats['head_gt'][0]:.6e}, std={self.stats['head_gt'][1]:.6e}")
        for aux_name in self.aux_channels:
            k = f'aux_{aux_name}'
            if k in self.stats:
                print(f"  {aux_name}: mean={self.stats[k][0]:.6e}, std={self.stats[k][1]:.6e}")

    def align_and_normalize(self, time, values, feat_type):
        values = safe_normalize(values)
        mask = np.isfinite(values)
        if not np.all(mask):
            time, values = time[mask], values[mask]
        if len(time) < 2:
            return np.zeros(len(self.target_time), dtype=np.float32)
        aligned = safe_normalize(
            interp1d(time, values, kind='linear', bounds_error=False,
                     fill_value='extrapolate')(self.target_time))
        if feat_type in self.stats:
            m, s = self.stats[feat_type]
            return ((aligned - m) / s).astype(np.float32)
        return aligned.astype(np.float32)

    def align_raw(self, time, values):
        values = safe_normalize(values)
        mask = np.isfinite(values)
        if not np.all(mask):
            time, values = time[mask], values[mask]
        if len(time) < 2:
            return np.zeros(len(self.target_time), dtype=np.float32)
        return safe_normalize(
            interp1d(time, values, kind='linear', bounds_error=False,
                     fill_value='extrapolate')(self.target_time))

    def decompose(self, signal):
        return self.decomposer.decompose(signal)

    def compute_rel_amp(self, signal):
        mask = self.steady_mask if self.config.time_range[0] < self.config.steady_start else None
        return compute_relative_amplitude(signal, mask)

    def get_gt_std(self):
        return self.stats.get('head_gt', (0, 1.0))[1]


# =============================================================================
# [6] 样本构建器
# =============================================================================
MATERIAL_TO_IDX = {'standard': 0, 'DualDensity': 1, 'SizeDiff': 2}


class GTPairBuilder:
    def __init__(self, config: GTConfigV13, aux_channels: List[str],
                 prep: GTPreprocessor):
        self.config = config
        self.aux_channels = aux_channels
        self.prep = prep
        self.coh_encoder = AdaptiveCohesionEncoder(COHESION_VALUES)
    
    def _build_condition(self, mat, coh):
        mat_oh = np.zeros(3, dtype=np.float32)
        mat_oh[MATERIAL_TO_IDX[mat]] = 1.0
        coh_enc = self.coh_encoder.encode(COHESION_VALUES[coh])
        inter_low = mat_oh * coh_enc[3]
        inter_high = mat_oh * coh_enc[4]
        return np.concatenate([mat_oh, coh_enc, inter_low, inter_high])
    
    def _process_gt(self, data, coarse_s, fine_s):
        t_c = data['gt'][coarse_s]['time']
        t_f = data['gt'][fine_s]['time']
        v_c = data['gt'][coarse_s]['values']
        v_f = data['gt'][fine_s]['values']
        
        raw_c = self.prep.align_raw(t_c, v_c)
        raw_f = self.prep.align_raw(t_f, v_f)
        low_c, high_c = self.prep.decompose(raw_c)
        low_f, high_f = self.prep.decompose(raw_f)
        rel_amp = self.prep.compute_rel_amp(raw_f)
        use_freq = rel_amp >= self.config.rel_amp_threshold
        osc_decay = compute_osc_decay_scale(rel_amp, self.config.osc_decay_factor,
                                             self.config.osc_decay_min)
        eps = 1e-8
        ratio = ((raw_f - raw_c) / (np.abs(raw_c) + eps)).astype(np.float32)
        ratio_low = ((low_f - low_c) / (np.abs(low_c) + eps)).astype(np.float32)
        
        return {
            'raw_c': raw_c, 'low_c': low_c, 'high_c': high_c,
            'raw_f': raw_f, 'low_f': low_f, 'high_f': high_f,
            'ratio': ratio, 'ratio_low': ratio_low,
            'rel_amp': rel_amp, 'use_freq': use_freq, 'osc_decay': osc_decay,
        }
    
    def _build_input_channels(self, data, scale):
        """
        构建 V13 GT 输入通道.
        
        前 N 通道: 辅助通道 (按 aux_channels 顺序)
        最后 1 通道: GT 目标通道
        """
        channels = []
        
        # 辅助通道
        for aux_name in self.aux_channels:
            aux = data['aux_channels'].get(aux_name)
            if aux is not None and scale in aux['data']:
                channels.append(self.prep.align_and_normalize(
                    aux['time'][scale], aux['data'][scale], f'aux_{aux_name}'))
            else:
                channels.append(np.zeros(self.config.n_time_points, dtype=np.float32))
        
        # GT 目标通道
        gt = data['gt'][scale]
        channels.append(self.prep.align_and_normalize(
            gt['time'], gt['values'], 'head_gt'))
        
        return np.stack(channels).astype(np.float32)
    
    def build_training_samples(self, all_data):
        samples = []
        for key, data in all_data.items():
            mat, coh = data['material'], data['cohesion']
            gt_proc = self._process_gt(data, coarse_s=4, fine_s=2)
            samples.append({
                'material': mat, 'cohesion': coh,
                'input_f4': self._build_input_channels(data, 4),
                'input_f2': self._build_input_channels(data, 2),
                'condition_base': self._build_condition(mat, coh),
                'weight': COHESION_WEIGHTS.get(coh, 1.0),
                'gt': gt_proc,
            })
        print(f"\n[Training samples] {len(samples)} pairs (F4→F2)")
        return samples
    
    def build_inference_samples(self, all_data):
        samples = []
        for key, data in all_data.items():
            mat, coh = data['material'], data['cohesion']
            gt_proc = self._process_gt(data, coarse_s=2, fine_s=1)
            # F4 raw for Richardson
            t4, v4 = data['gt'][4]['time'], data['gt'][4]['values']
            gt_proc['f4_raw'] = self.prep.align_raw(t4, v4)
            samples.append({
                'material': mat, 'cohesion': coh,
                'input_f2': self._build_input_channels(data, 2),
                'condition_base': self._build_condition(mat, coh),
                'gt': gt_proc,
            })
        print(f"\n[Inference samples] {len(samples)} pairs (F2→F1)")
        return samples


# =============================================================================
# [7] Dataset
# =============================================================================
class GTTrainDataset(Dataset):
    def __init__(self, samples):
        self.samples = samples
    
    def __len__(self):
        return len(self.samples)
    
    def __getitem__(self, idx):
        s = self.samples[idx]
        g = s['gt']
        return {
            'input_f4': torch.FloatTensor(s['input_f4']),
            'input_f2': torch.FloatTensor(s['input_f2']),
            'condition_base': torch.FloatTensor(s['condition_base']),
            'weight': torch.FloatTensor([s['weight']]),
            'material': s['material'], 'cohesion': s['cohesion'],
            'raw_c': torch.FloatTensor(g['raw_c']),
            'low_c': torch.FloatTensor(g['low_c']),
            'high_c': torch.FloatTensor(g['high_c']),
            'raw_f': torch.FloatTensor(g['raw_f']),
            'low_f': torch.FloatTensor(g['low_f']),
            'high_f': torch.FloatTensor(g['high_f']),
            'ratio': torch.FloatTensor(g['ratio']),
            'ratio_low': torch.FloatTensor(g['ratio_low']),
            'rel_amp': torch.FloatTensor([g['rel_amp']]),
            'use_freq': torch.BoolTensor([g['use_freq']]),
            'osc_decay': torch.FloatTensor([g['osc_decay']]),
        }


class GTInferDataset(Dataset):
    def __init__(self, samples):
        self.samples = samples
    
    def __len__(self):
        return len(self.samples)
    
    def __getitem__(self, idx):
        s = self.samples[idx]
        g = s['gt']
        return {
            'input_f2': torch.FloatTensor(s['input_f2']),
            'condition_base': torch.FloatTensor(s['condition_base']),
            'material': s['material'], 'cohesion': s['cohesion'],
            'raw_c': torch.FloatTensor(g['raw_c']),
            'low_c': torch.FloatTensor(g['low_c']),
            'high_c': torch.FloatTensor(g['high_c']),
            'raw_f': torch.FloatTensor(g['raw_f']),
            'low_f': torch.FloatTensor(g['low_f']),
            'high_f': torch.FloatTensor(g['high_f']),
            'ratio': torch.FloatTensor(g['ratio']),
            'ratio_low': torch.FloatTensor(g['ratio_low']),
            'f4_raw': torch.FloatTensor(g['f4_raw']),
            'rel_amp': torch.FloatTensor([g['rel_amp']]),
            'use_freq': torch.BoolTensor([g['use_freq']]),
            'osc_decay': torch.FloatTensor([g['osc_decay']]),
        }


# =============================================================================
# [8] 训练器
# =============================================================================
class GTTrainer:
    def __init__(self, config: GTConfigV13, model: GTPredictorV13):
        self.config = config
        self.model = model.to(config.device)
        self.opt = torch.optim.AdamW(model.parameters(), lr=config.learning_rate,
                                      weight_decay=config.weight_decay)
        self.sched = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
            self.opt, T_0=100, T_mult=2)
        
        time = np.linspace(config.time_range[0], config.time_range[1], config.n_time_points)
        self.steady_mask = torch.FloatTensor(
            (time >= config.steady_start).astype(float)).to(config.device)
        self.transient_mask = torch.FloatTensor(
            (time < config.steady_start).astype(float)).to(config.device)
        self.history = {'total': [], 'gt': [], 'cycle': []}
    
    def train_step(self, batch):
        self.model.train()
        dev = self.config.device
        cfg = self.config
        
        f4 = batch['input_f4'].to(dev)
        f2 = batch['input_f2'].to(dev)
        cb = batch['condition_base'].to(dev)
        weights = batch['weight'].to(dev)
        B = f4.shape[0]
        
        rel_amp = batch['rel_amp'].to(dev)
        use_freq = batch['use_freq'].to(dev)
        osc_decay = batch['osc_decay'].to(dev)
        raw_c = batch['raw_c'].to(dev)
        low_c = batch['low_c'].to(dev)
        high_c = batch['high_c'].to(dev)
        raw_f = batch['raw_f'].to(dev)
        gt_ratio = batch['ratio'].to(dev)
        gt_ratio_low = batch['ratio_low'].to(dev)
        
        sr_05 = torch.full((B, 1), 0.5, device=dev)
        sr_20 = torch.full((B, 1), 2.0, device=dev)
        
        # --- Forward (F4 input) ---
        feats_fwd, _ = self.model.encode(f4, cb, rel_amp, sr_05)
        ratio_fwd, ratio_low_fwd = self.model.decode_ratios(feats_fwd)
        direct_fwd = self.model.decode_direct(feats_fwd)
        dr = direct_fwd * cfg.direct_res_scale
        
        pred_f, _ = self.model.apply_correction(
            ratio_fwd, ratio_low_fwd, use_freq, osc_decay,
            raw_c, low_c, high_c, dr)
        
        # --- Head loss ---
        w = weights.squeeze(-1)
        loss_r = torch.tensor(0.0, device=dev)
        freq_on = cfg.use_freq_decomp
        
        for i in range(B):
            if freq_on and use_freq[i]:
                loss_r += w[i] * F_torch.mse_loss(ratio_low_fwd[i], gt_ratio_low[i])
            else:
                loss_r += w[i] * F_torch.mse_loss(ratio_fwd[i], gt_ratio[i])
        loss_r = loss_r / B
        
        loss_recon = (weights * ((pred_f - raw_f)**2).mean(1, keepdim=True)).mean()
        loss_steady = (self.steady_mask * (pred_f - raw_f)**2).mean()
        loss_trans = (self.transient_mask * (pred_f - raw_f)**2).mean()
        loss_shape = F_torch.mse_loss(ratio_fwd[:, 1:] - ratio_fwd[:, :-1],
                                       gt_ratio[:, 1:] - gt_ratio[:, :-1])
        loss_reg = ratio_fwd.abs().mean()
        
        h_loss = (loss_r + loss_recon
                  + cfg.weight_steady * loss_steady
                  + cfg.weight_transient * loss_trans
                  + cfg.weight_shape * loss_shape
                  + cfg.weight_ratio_reg * loss_reg)
        
        total_loss = h_loss
        
        # --- Direct decoder regularization ---
        if cfg.weight_direct > 0:
            total_loss += cfg.weight_direct * direct_fwd.abs().mean()
        
        # --- Cycle consistency ---
        loss_cycle = torch.tensor(0.0, device=dev)
        if cfg.weight_cycle > 0:
            feats_cyc, _ = self.model.encode(f2, cb, rel_amp, sr_20)
            ratio_cyc, ratio_low_cyc = self.model.decode_ratios(feats_cyc)
            dr_cyc = self.model.decode_direct(feats_cyc) * cfg.direct_res_scale
            
            pred_back, _ = self.model.apply_correction(
                ratio_cyc, ratio_low_cyc, use_freq, osc_decay,
                pred_f.detach(), low_c, high_c, dr_cyc)
            loss_cycle = F_torch.mse_loss(pred_back, raw_c)
            total_loss += cfg.weight_cycle * loss_cycle
        
        # --- Cohesion-adaptive loss ---
        if cfg.weight_cohesion_adapt > 0:
            loss_coh = torch.tensor(0.0, device=dev)
            n_high = 0
            for i in range(B):
                if weights[i] > 2.0:
                    raw_fn = (raw_f[i] - raw_f[i].mean()) / (raw_f[i].std() + 1e-8)
                    pred_fn = (pred_f[i] - pred_f[i].mean()) / (pred_f[i].std() + 1e-8)
                    loss_coh += F_torch.mse_loss(pred_fn, raw_fn)
                    n_high += 1
            if n_high > 0:
                loss_coh = loss_coh / n_high
            total_loss += cfg.weight_cohesion_adapt * loss_coh
        
        if torch.isnan(total_loss):
            return {'total': 0.0, 'gt': 0.0, 'cycle': 0.0}
        
        self.opt.zero_grad()
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
        self.opt.step()
        
        return {'total': total_loss.item(), 'gt': h_loss.item(),
                'cycle': loss_cycle.item()}
    
    def train(self, train_ds):
        cfg = self.config
        print(f"\n{'='*60}")
        print(f"Training V13 GT Predictor")
        self.model.print_architecture()
        print(f"{'='*60}")
        
        loader = DataLoader(train_ds, batch_size=cfg.batch_size, shuffle=True)
        
        for ep in range(cfg.n_epochs):
            losses = {'total': 0.0, 'gt': 0.0, 'cycle': 0.0}
            for batch in loader:
                l = self.train_step(batch)
                for k, v in l.items():
                    losses[k] += v
            n_batches = len(loader)
            for k in losses:
                losses[k] /= max(n_batches, 1)
                self.history[k].append(losses[k])
            self.sched.step()
            
            if (ep + 1) % 50 == 0:
                print(f"  Epoch {ep+1:3d}: loss={losses['total']:.6f}, "
                      f"gt={losses['gt']:.6f}")
        
        print(f"  Final: {self.history['total'][-1]:.6f}")
    
    def save(self, path):
        torch.save({'model': self.model.state_dict(), 'history': self.history,
                     'config': self.config}, path)
        print(f"  Saved model to {path}")


# =============================================================================
# [9] 评估器
# =============================================================================
class GTEvaluator:
    def __init__(self, config: GTConfigV13, model: GTPredictorV13,
                 prep: GTPreprocessor):
        self.config = config
        self.model = model.to(config.device).eval()
        self.prep = prep
        eval_mask = None
        if config.time_range[0] < config.steady_start:
            time = np.linspace(config.time_range[0], config.time_range[1],
                               config.n_time_points)
            eval_mask = time >= config.eval_start_time
        self.metrics = MetricsCalculator(eval_mask)
    
    def evaluate(self, loader):
        results = []
        dev = self.config.device
        
        with torch.no_grad():
            for batch in loader:
                f2 = batch['input_f2'].to(dev)
                cb = batch['condition_base'].to(dev)
                B = f2.shape[0]
                sr = torch.full((B, 1), 0.5, device=dev)
                rel_amp = batch['rel_amp'].to(dev)
                
                feats, _ = self.model.encode(f2, cb, rel_amp, sr)
                ratio, ratio_low = self.model.decode_ratios(feats)
                direct = self.model.decode_direct(feats)
                
                for i in range(B):
                    uf = batch['use_freq'].to(dev)
                    od = batch['osc_decay'].to(dev)
                    raw_c = batch['raw_c'][i:i+1].to(dev)
                    low_c = batch['low_c'][i:i+1].to(dev)
                    high_c = batch['high_c'][i:i+1].to(dev)
                    
                    dr_i = direct[i:i+1] * self.config.direct_res_scale
                    pred_i, _ = self.model.apply_correction(
                        ratio[i:i+1], ratio_low[i:i+1] if ratio_low is not None else None,
                        uf[i:i+1], od[i:i+1],
                        raw_c, low_c, high_c, dr_i)
                    
                    pred_np = pred_i.squeeze(0).cpu().numpy()
                    f1_np = batch['raw_f'][i].numpy()
                    f2_np = batch['raw_c'][i].numpy()
                    f4_np = batch['f4_raw'][i].numpy()
                    
                    # Richardson extrapolation baseline
                    richardson = 2.0 * f2_np - f4_np
                    
                    m = self.metrics.compute(pred_np, f1_np, richardson)
                    m.update({
                        'material': batch['material'][i],
                        'cohesion': batch['cohesion'][i],
                        'config': f"{batch['material'][i]}_se{batch['cohesion'][i]}",
                        'pred': pred_np,
                        'f1': f1_np,
                        'f2': f2_np,
                        'f4': f4_np,
                        'richardson': richardson,
                    })
                    results.append(m)
        
        return results


# =============================================================================
# [10] LOOCV 运行函数
# =============================================================================
def run_single(config: GTConfigV13, aux_channels: List[str],
               sta_dir: str, ext_dir: str, output_dir: str):
    """执行 12-fold LOOCV 并返回结果"""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    # 加载数据
    all_data = load_gt_data(aux_channels, sta_dir, ext_dir)
    if not all_data:
        print("[ERROR] No data loaded!")
        return None
    
    configs = list(all_data.keys())
    all_results = []
    
    print(f"\n{'='*60}")
    print(f"  12-fold LOOCV — GT V13")
    print(f"  Aux mode: {config.aux_mode} ({len(aux_channels)} channels)")
    print(f"  Configs: {len(configs)}")
    print(f"{'='*60}")
    
    for fold_idx, test_key in enumerate(configs):
        print(f"\n--- Fold {fold_idx+1}/{len(configs)}: test={test_key} ---")
        
        set_seed()
        
        # Split
        train_data = {k: v for k, v in all_data.items() if k != test_key}
        test_data = {test_key: all_data[test_key]}
        
        # Preprocessor
        prep = GTPreprocessor(config, aux_channels)
        prep.fit(train_data)
        config.direct_res_scale = prep.get_gt_std()
        
        # Build samples
        builder = GTPairBuilder(config, aux_channels, prep)
        train_samples = builder.build_training_samples(train_data)
        infer_samples = builder.build_inference_samples(test_data)
        
        train_ds = GTTrainDataset(train_samples)
        infer_ds = GTInferDataset(infer_samples)
        infer_loader = DataLoader(infer_ds, batch_size=1, shuffle=False)
        
        # Model
        model = GTPredictorV13(config)
        if fold_idx == 0:
            model.print_architecture()
        
        # Train
        trainer = GTTrainer(config, model)
        trainer.train(train_ds)
        
        # Evaluate
        evaluator = GTEvaluator(config, model, prep)
        fold_results = evaluator.evaluate(infer_loader)
        
        for r in fold_results:
            er = r['error_ratio']
            corr = r['correlation']
            status = "✓" if r['success'] else "✗"
            print(f"  {status} {r['config']}: ER={er:.3f}, corr={corr:.3f}")
        
        all_results.extend(fold_results)
        
        # Save fold model
        fold_dir = out_path / f"fold_{fold_idx:02d}_{test_key}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        trainer.save(fold_dir / "model.pt")
    
    # Summary
    ers = [r['error_ratio'] for r in all_results]
    corrs = [r['correlation'] for r in all_results]
    n_success = sum(1 for r in all_results if r['success'])
    
    print(f"\n{'='*60}")
    print(f"  GT V13 LOOCV Summary ({config.aux_mode} mode)")
    print(f"  Avg ER:    {np.mean(ers):.4f} ± {np.std(ers):.4f}")
    print(f"  Avg Corr:  {np.mean(corrs):.4f}")
    print(f"  Success:   {n_success}/{len(all_results)}")
    print(f"{'='*60}")
    
    # Save summary
    summary = {
        'mode': config.aux_mode,
        'aux_channels': aux_channels,
        'avg_er': float(np.mean(ers)),
        'std_er': float(np.std(ers)),
        'avg_corr': float(np.mean(corrs)),
        'n_success': n_success,
        'n_total': len(all_results),
        'per_config': {
            r['config']: {
                'error_ratio': float(r['error_ratio']),
                'correlation': float(r['correlation']),
                'success': bool(r['success']), # 保险起见，把 success 也转为 Python 原生 bool
                'rmse': float(r['rmse']),
            }
            for r in all_results
        }
    }
    with open(out_path / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    return all_results


# =============================================================================
# [11] V12 Baseline 对比运行
# =============================================================================
def run_v12_gt_baseline(sta_dir, ext_dir, output_dir):
    """
    运行 V12 GT baseline (全 8 通道, 完整架构) 作为对比.
    通过导入 V12 模块实现.
    """
    try:
        from unified_predictor_v12 import PredictorV12, make_gt_config
        from run_predictor_v12 import (
            TARGET_REGISTRY, load_all_data, UnifiedPreprocessor,
            UnifiedPairBuilder, UnifiedTrainDataset, UnifiedInferDataset,
            UnifiedTrainer, UnifiedEvaluator
        )
        
        spec = TARGET_REGISTRY['gt']
        config = make_gt_config()
        
        print(f"\n{'#'*60}")
        print(f"# V12 GT Baseline (full 8 aux channels)")
        print(f"{'#'*60}")
        
        all_data = load_all_data(spec, sta_dir, ext_dir)
        if not all_data:
            return None
        
        configs = list(all_data.keys())
        all_results = []
        
        for fold_idx, test_key in enumerate(configs):
            set_seed()
            train_data = {k: v for k, v in all_data.items() if k != test_key}
            test_data = {test_key: all_data[test_key]}
            
            prep = UnifiedPreprocessor(config, spec)
            prep.fit(train_data)
            config.direct_res_scale = prep.get_target_std('gt')
            
            builder = UnifiedPairBuilder(config, spec, prep)
            train_samples = builder.build_training_samples(train_data)
            infer_samples = builder.build_inference_samples(test_data)
            
            train_ds = UnifiedTrainDataset(train_samples, spec)
            infer_ds = UnifiedInferDataset(infer_samples, spec)
            infer_loader = DataLoader(infer_ds, batch_size=1, shuffle=False)
            
            model = PredictorV12(config)
            trainer = UnifiedTrainer(config, spec, model)
            trainer.train(train_ds)
            
            evaluator = UnifiedEvaluator(config, spec, model, prep)
            fold_results = evaluator.evaluate(infer_loader)
            all_results.extend(fold_results)
        
        ers = [r['heads']['gt']['error_ratio'] for r in all_results]
        print(f"\n  V12 Baseline: Avg ER = {np.mean(ers):.4f}")
        return all_results
    
    except ImportError:
        print("[WARN] Cannot import V12 modules for baseline comparison")
        return None


# =============================================================================
# [12] 消融汇总
# =============================================================================
def summarize_ablation(results_dict: Dict[str, list], aux_channels: List[str]):
    """打印消融实验结果汇总"""
    print(f"\n{'='*70}")
    print(f"  GT V13 Ablation Summary")
    print(f"{'='*70}")
    print(f"{'Experiment':<30} {'Avg ER':>8} {'Std':>8} {'ΔER':>8} {'Succ':>6}")
    print(f"{'-'*70}")
    
    baseline_er = None
    for exp_name, results in sorted(results_dict.items()):
        ers = [r['error_ratio'] for r in results]
        avg_er = np.mean(ers)
        std_er = np.std(ers)
        n_succ = sum(1 for r in results if r['success'])
        
        if exp_name == 'baseline':
            baseline_er = avg_er
            delta = "  ---"
        else:
            delta = f"{avg_er - baseline_er:+.4f}" if baseline_er is not None else "  N/A"
        
        print(f"  {exp_name:<28} {avg_er:8.4f} {std_er:8.4f} {delta:>8} {n_succ:3d}/{len(results)}")
    
    print(f"{'='*70}")


# =============================================================================
# [13] Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description='V13 GT-Optimized Predictor')
    parser.add_argument('--mode', type=str, default='auto',
                        choices=['minimal', 'enhanced', 'auto'],
                        help='Auxiliary channel mode')
    parser.add_argument('--sta_dir', type=str, default='sta_results')
    parser.add_argument('--ext_dir', type=str, default='extended_features_results')
    parser.add_argument('--output_dir', type=str, default='gt_output_v13')
    parser.add_argument('--compare_v12', action='store_true',
                        help='Also run V12 baseline for comparison')
    parser.add_argument('--ablation_sweep', action='store_true',
                        help='Run V13 internal ablation experiments')
    parser.add_argument('--n_epochs', type=int, default=None,
                        help='Override number of training epochs')
    
    args = parser.parse_args()
    
    # === 确定辅助通道模式 ===
    if args.mode == 'auto':
        mode, aux_channels = auto_detect_mode(args.ext_dir)
    elif args.mode == 'enhanced':
        # 检查文件是否存在
        available_2nd = detect_second_order_features(args.ext_dir)
        if len(available_2nd) < len(AUX_CHANNELS_ENHANCED) - 2:
            print(f"\n[WARNING] Enhanced mode requested but only {len(available_2nd)} "
                  f"second-order features available")
            print(f"  Available: {available_2nd}")
            print(f"  Missing features will be zero-padded")
        aux_channels = list(AUX_CHANNELS_ENHANCED)
        mode = 'enhanced'
    else:
        aux_channels = list(AUX_CHANNELS_MINIMAL)
        mode = 'minimal'
    
    # === 创建配置 ===
    config = GTConfigV13(aux_mode=mode)
    config.set_aux_channels(aux_channels)
    if args.n_epochs is not None:
        config.n_epochs = args.n_epochs
    config.validate()
    
    print(f"\n{'='*70}")
    print(f"  GT Predictor V13 — {mode.upper()} mode")
    print(f"  Aux channels ({len(aux_channels)}): {aux_channels}")
    print(f"  Architecture: input_concat, wo_2nd_order_diff, freq_decomp ON")
    print(f"  Output: {args.output_dir}")
    print(f"  Device: {config.device}")
    print(f"{'='*70}")
    
    # === V12 对比 ===
    if args.compare_v12:
        v12_output = f"{args.output_dir}/v12_baseline"
        run_v12_gt_baseline(args.sta_dir, args.ext_dir, v12_output)
    
    # === 消融模式 ===
    if args.ablation_sweep:
        experiments = list_gt_v13_ablation_experiments(aux_channels)
        all_exp_results = {}
        
        for exp_name, exp_info in experiments.items():
            print(f"\n{'#'*60}")
            print(f"# Ablation: {exp_name}")
            print(f"# Hypothesis: {exp_info['hypothesis']}")
            print(f"{'#'*60}")
            
            exp_config = GTConfigV13(aux_mode=mode)
            exp_config.set_aux_channels(aux_channels)
            if exp_info['mask'] != tuple([True] * len(aux_channels)):
                exp_config.channel_mask = exp_info['mask']
            if args.n_epochs is not None:
                exp_config.n_epochs = args.n_epochs
            exp_config.validate()
            
            exp_dir = f"{args.output_dir}/ablation_{exp_name}"
            set_seed()
            results = run_single(exp_config, aux_channels,
                                 args.sta_dir, args.ext_dir, exp_dir)
            if results is not None:
                all_exp_results[exp_name] = results
        
        summarize_ablation(all_exp_results, aux_channels)
        
        # Save summary
        summary = {}
        for exp_name, results in all_exp_results.items():
            ers = [r['error_ratio'] for r in results]
            summary[exp_name] = {
                'avg_er': float(np.mean(ers)),
                'std_er': float(np.std(ers)),
                'n_success': sum(1 for r in results if r['success']),
            }
        with open(f"{args.output_dir}/ablation_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        return
    
    # === 标准运行 ===
    set_seed()
    results = run_single(config, aux_channels,
                         args.sta_dir, args.ext_dir, args.output_dir)
    
    if results is not None:
        print(f"\n{'='*70}")
        print(f"  GT V13 Complete — {mode} mode")
        print(f"  Average ER: {np.mean([r['error_ratio'] for r in results]):.4f}")
        print(f"{'='*70}")


if __name__ == "__main__":
    main()
