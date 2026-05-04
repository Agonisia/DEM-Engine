#!/usr/bin/env python3
"""
V14 统一可视化模块 (学术报告风格)
==========================================

V13 → V14 变更:
  1. 采用学术报告风格: Times New Roman, serif, markers with markevery, hatched bars
  2. 全英文标签 (无 CJK 依赖)
  3. 仅输出 PNG (600 dpi)
  4. 对全部材料 (standard, SizeDiff, DualDensity) 出图
  5. 移除 freq_decomp 和 input_concat 相关消融变体
  6. 使用科学计数法 y 轴

图表清单:
  1. predictions_{hname}/       - 各配置预测曲线 (标题显示 RMSE)
  2. by_material_{hname}/       - Surface Energy vs RMSE
  3. hierarchy_{hname}.png      - 多尺度层级总览
  4. summary_{hname}.png        - 汇总面板 (RMSE/MAE/SMAPE)
  5. csv_{hname}/               - 全量 CSV 导出
  6. ablation_bar_{hname}.png   - 输入消融 ΔRMSE 条状图
  7. ablation_groups_{hname}.png- 分组消融对比
  8. ablation_heatmap.png       - 跨目标热力图 (ΔRMSE)
  9. arch_ablation_*            - 架构消融 (ΔRMSE)
"""

import numpy as np
import pandas as pd
import pickle
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from scipy.ndimage import uniform_filter1d
import argparse
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

# =============================================================================
# 全局学术风格
# =============================================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'serif'],
    'mathtext.fontset': 'stix',
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 12,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'lines.linewidth': 2.0,
    'axes.unicode_minus': False,
})

# 经典学术配色 + marker
COLOR_F1 = 'k'        # 黑: target (F1)
COLOR_F2 = 'b'        # 蓝: input (F2)
COLOR_F4 = 'g'        # 绿: coarse (F4)
COLOR_ML = 'r'        # 红: ML prediction
COLOR_RE = '#9467bd'   # 紫: Richardson

MARK_EVERY = 10
MARKER_SIZE = 6
DPI = 600

# =============================================================================
# 目标类型元信息
# =============================================================================
TARGET_META = {
    'ke':      {'head_names': ('ke_total',), 'n_pts': 200, 'trange': (0.0, 3.0),
                'unit': 'J', 'ylabel': 'Kinetic Energy [J]', 'label': 'KE Total',
                'head_labels': {'ke_total': 'Kinetic Energy (Total)'}},
    'ke_comp': {'head_names': ('trans', 'rot'), 'n_pts': 200, 'trange': (0.0, 3.0),
                'unit': 'J', 'ylabel': 'Kinetic Energy [J]', 'label': 'KE Components',
                'head_labels': {'trans': 'Translational Kinetic Energy',
                                'rot': 'Rotational Kinetic Energy'}},
    'gt':      {'head_names': ('gt',), 'n_pts': 200, 'trange': (0.0, 3.0),
                'unit': r'm$^2$/s$^2$', 'ylabel': r'Granular Temperature [m$^2$/s$^2$]',
                'label': 'Granular Temperature',
                'head_labels': {'gt': 'Granular Temperature'}},
}

# 架构消融变体 (V14: 移除了 wo_freq_decomp 和 wo_film_input_concat)
ARCH_VARIANT_COLORS = {
    'full_model':        '#2c3e50',
    'wo_film':           '#e74c3c',
    'wo_local_features': '#9b59b6',
    'wo_2nd_order_diff': '#1abc9c',
}
ARCH_VARIANT_LABELS = {
    'full_model':        'Full model',
    'wo_film':           'w/o FiLM',
    'wo_local_features': 'w/o local features',
    'wo_2nd_order_diff': 'w/o 2nd-order diff',
}

# 输入通道英文标签
CHANNEL_LABELS = {
    'torque': 'Torque', 'KE_total': 'Total KE',
    'v_theta_mean': 'Mean tangential vel.', 'speed_mean': 'Mean speed',
    'omega_mean': 'Mean angular vel.', 'T_total': 'Granular temp.',
    'KE_rot_sum': 'Rotational KE',
    'all_kinematic': 'All kinematic', 'all_energy': 'All energy',
}
GROUP_LABELS = {
    'baseline': 'Full model\n(baseline)',
    'remove_all_kinematic': 'Remove all\nkinematic',
    'remove_all_energy': 'Remove all\nenergy',
    'only_torque': 'Torque only',
}
HEAD_LABELS = {
    'ke_total': 'Total KE', 'trans': 'Translational KE',
    'rot': 'Rotational KE', 'gt': 'Granular Temperature',
}


# =============================================================================
# 配置
# =============================================================================
@dataclass
class VizConfig:
    steady_start: float = 0.6
    steady_end: float = 3.0
    cohesion_values: Dict[str, float] = field(default_factory=lambda: {
        '000': 0.0, '005': 0.05, '010': 0.1, '015': 0.15, '020': 0.2
    })
    smooth_window: int = 5
    csv_nth: int = 15


# =============================================================================
# 辅助函数 (学术风格)
# =============================================================================
def _smooth(arr, w=5):
    if arr is None or w <= 1:
        return arr
    return uniform_filter1d(arr, size=w, mode='nearest')


def _apply_sci_yaxis(ax):
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_fontsize(13)


def _snap_ylim_to_top_tick(ax):
    fig = ax.get_figure()
    fig.canvas.draw()
    ticks = ax.get_yticks()
    ylo, yhi = ax.get_ylim()
    above = [t for t in ticks if t > yhi - 1e-15]
    if above:
        ax.set_ylim(ylo, above[0])
    else:
        visible = [t for t in ticks if t >= ylo]
        if visible:
            ax.set_ylim(ylo, visible[-1])


def _academic_style(ax):
    ax.grid(True, linestyle=':', alpha=0.5, linewidth=0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(1.2)


def _decimate(arr, nth):
    if arr is None or nth <= 1:
        return arr
    out = np.full_like(arr, np.nan, dtype=float)
    out[::nth] = arr[::nth]
    return out


def _save(fig, path):
    fig.savefig(str(path), bbox_inches='tight', dpi=DPI)
    plt.close(fig)


# =============================================================================
# 核心可视化器
# =============================================================================
class VisualizerV14:
    """V14 可视化器: 学术报告风格, 全英文, PNG only."""

    def __init__(self, target: str, output_dir: str, cfg: Optional[VizConfig] = None):
        self.target = target
        self.meta = TARGET_META[target]
        self.out = Path(output_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.cfg = cfg or VizConfig()

        t0, t1 = self.meta['trange']
        n = self.meta['n_pts']
        self.full_time = np.linspace(t0, t1, n)
        self.steady_mask = (self.full_time >= self.cfg.steady_start) & \
                           (self.full_time <= self.cfg.steady_end)
        steady_t = self.full_time[self.steady_mask]
        self.display_time = steady_t - steady_t[0]
        self.display_xmax = self.display_time[-1] if len(self.display_time) > 0 else 2.4
        self.x_ticks = np.arange(0.0, self.display_xmax + 0.01, 0.4)

    def _crop(self, arr):
        if arr is None:
            return None
        n_steady = int(self.steady_mask.sum())
        if len(arr) == n_steady:
            return arr
        if len(arr) == len(self.full_time):
            return arr[self.steady_mask]
        return None

    def _sc(self, arr):
        return _smooth(self._crop(arr), self.cfg.smooth_window)

    def load(self, input_dir: str) -> List[Dict]:
        for name in ['infer_results.pkl', 'infer_results_stable.pkl',
                      'infer_results_steady.pkl', 'infer_results_full.pkl']:
            p = Path(input_dir) / name
            if p.exists():
                with open(p, 'rb') as f:
                    results = pickle.load(f)
                print(f"  loaded {p} ({len(results)} results)")
                return results
        print(f"  No pkl found in {input_dir}")
        return []

    def _compute_ylim(self, results, hname):
        ymin, ymax = float('inf'), float('-inf')
        for r in results:
            h = r['heads'].get(hname)
            if h is None:
                continue
            for key in ['target', 'input', 'pred', 'f4']:
                arr = self._crop(h.get(key))
                if arr is not None:
                    c = arr[np.isfinite(arr)]
                    if len(c) > 0:
                        ymin = min(ymin, np.min(c))
                        ymax = max(ymax, np.max(c))
        if ymin == float('inf'):
            return (0, 1)
        margin = (ymax - ymin) * 0.05
        bottom = 0.0 if ymin >= 0 else ymin - margin
        return (bottom, ymax + margin)

    # =========================================================================
    # 1. 各配置预测曲线 — 学术报告风格
    # =========================================================================
    def _plot_single(self, h, ax, ylim=None):
        t = self.display_time
        # F4 (coarse)
        arr = self._sc(h.get('f4'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_F4, linestyle='-.', marker='s',
                    markevery=MARK_EVERY, markersize=MARKER_SIZE)
        # F2 (input)
        arr = self._sc(h.get('input'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_F2, linestyle='--', marker='o',
                    markevery=MARK_EVERY, markersize=MARKER_SIZE)
        # F1 (target)
        arr = self._sc(h.get('target'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_F1, linestyle='-')
        # Richardson
        arr = self._sc(h.get('richardson_pred'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_RE, linestyle=':', marker='d',
                    markevery=MARK_EVERY, markersize=MARKER_SIZE - 1)
        # ML prediction
        arr = self._sc(h.get('pred'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_ML, linestyle='-', marker='v',
                    markevery=MARK_EVERY, markersize=MARKER_SIZE)

        ax.set_xlim(0, self.display_xmax)
        ax.set_xticks(self.x_ticks)
        if ylim:
            ax.set_ylim(ylim)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=6, min_n_ticks=5))
        _apply_sci_yaxis(ax)
        _academic_style(ax)

    def _save_curve_legend(self, path):
        import matplotlib.lines as mlines
        fig, ax = plt.subplots(figsize=(10, 1.2))
        ax.set_axis_off()
        handles = [
            mlines.Line2D([], [], color=COLOR_F1, linestyle='-',
                          label='Scale factor = 1 (Target)'),
            mlines.Line2D([], [], color=COLOR_ML, linestyle='-', marker='v',
                          markersize=MARKER_SIZE, label='ML Prediction'),
            mlines.Line2D([], [], color=COLOR_F2, linestyle='--', marker='o',
                          markersize=MARKER_SIZE, label='Scale factor = 2 (Input)'),
            mlines.Line2D([], [], color=COLOR_F4, linestyle='-.', marker='s',
                          markersize=MARKER_SIZE, label='Scale factor = 4 (Coarse)'),
            mlines.Line2D([], [], color=COLOR_RE, linestyle=':', marker='d',
                          markersize=MARKER_SIZE - 1, label='Richardson Extrap.'),
        ]
        ax.legend(handles=handles, loc='center', ncol=5,
                  frameon=True, edgecolor='black')
        _save(fig, path)

    def plot_predictions(self, results):
        for hname in self.meta['head_names']:
            tag = f"predictions_{hname}"
            print(f"\n[{tag}]")
            subdir = self.out / tag
            subdir.mkdir(parents=True, exist_ok=True)

            ylim = self._compute_ylim(results, hname)

            for r in results:
                h = r['heads'].get(hname)
                if h is None:
                    continue
                mat, coh = r['material'], r['cohesion']
                se = self.cfg.cohesion_values.get(coh, 0.0)
                rmse = h['rmse']
                hlabel = self.meta.get('head_labels', {}).get(hname, hname)

                fig, ax = plt.subplots(figsize=(8, 5))
                self._plot_single(h, ax, ylim=ylim)
                ax.set_title(
                    f'{hlabel} | Surface Energy = {se:.3f} J/m'
                    r'$^2$'
                    f' | RMSE = {rmse:.6f}',
                    pad=20)
                ax.set_xlabel('Time [s]')
                ax.set_ylabel(self.meta['ylabel'])
                _snap_ylim_to_top_tick(ax)
                plt.tight_layout()
                fname = f"pred_{mat}_{coh}.png"
                _save(fig, subdir / fname)
                print(f"  {fname}")

            self._save_curve_legend(subdir / "legend.png")
            print(f"  legend.png  ({len(results)} plots)")

    # =========================================================================
    # 2. Surface Energy 对比图 — RMSE
    # =========================================================================
    def plot_by_material(self, results):
        for hname in self.meta['head_names']:
            tag = f"by_material_{hname}"
            print(f"\n[{tag}]")
            subdir = self.out / tag
            subdir.mkdir(parents=True, exist_ok=True)

            materials = sorted(set(r['material'] for r in results))

            all_rmse = []
            for r in results:
                h = r['heads'].get(hname)
                if h is None:
                    continue
                for k in ['baseline_rmse', 'rmse', 'richardson_rmse']:
                    v = h.get(k)
                    if v is not None and np.isfinite(v):
                        all_rmse.append(v)
            y_max = max(all_rmse) * 1.15 if all_rmse else 1.0

            all_se = [self.cfg.cohesion_values.get(r['cohesion'], 0) for r in results]
            se_max = max(all_se) if all_se else 1.0
            x_margin = se_max * 0.05

            for mat in materials:
                mat_r = sorted(
                    [r for r in results if r['material'] == mat],
                    key=lambda r: self.cfg.cohesion_values.get(r['cohesion'], 0))

                se_vals = [self.cfg.cohesion_values.get(r['cohesion'], 0) for r in mat_r]
                base_rmse = [r['heads'][hname]['baseline_rmse'] for r in mat_r]
                pred_rmse = [r['heads'][hname]['rmse'] for r in mat_r]
                rich_rmse = [r['heads'][hname].get('richardson_rmse', np.nan) for r in mat_r]

                fig, ax = plt.subplots(figsize=(8, 5.5))
                ax.plot(se_vals, base_rmse, 'o--', color='#7f7f7f', ms=8, lw=2, label='F2 Baseline')
                ax.plot(se_vals, pred_rmse, 's-', color=COLOR_ML, ms=8, lw=2, label='ML Prediction')
                ax.plot(se_vals, rich_rmse, '^-.', color=COLOR_RE, ms=8, lw=2, label='Richardson')

                hlabel = self.meta.get('head_labels', {}).get(hname, hname)
                ax.set_title(f'{hlabel} | {mat}: Surface Energy vs RMSE', pad=20)
                ax.set_xlabel(r'Surface Energy [J/m$^2$]')
                ax.set_ylabel(f'RMSE [{self.meta["unit"]}]')
                ax.set_xlim(-x_margin, se_max + x_margin)
                ax.set_ylim(bottom=0, top=y_max)
                ax.legend(loc='best')
                _apply_sci_yaxis(ax)
                _academic_style(ax)
                plt.tight_layout()
                _save(fig, subdir / f"se_{mat}.png")
                print(f"  se_{mat}.png")

                # CSV export
                df = pd.DataFrame({
                    'cohesion': [r['cohesion'] for r in mat_r],
                    'surface_energy': se_vals,
                    'baseline_rmse': base_rmse,
                    'ml_rmse': pred_rmse,
                    'richardson_rmse': rich_rmse,
                })
                df.to_csv(subdir / f"se_{mat}.csv", index=False, float_format='%.6g')

    # =========================================================================
    # 3. 多尺度层级总览
    # =========================================================================
    def plot_hierarchy(self, results):
        for hname in self.meta['head_names']:
            tag = f"hierarchy_{hname}"
            print(f"\n[{tag}]")

            n = len(results)
            if n == 0:
                continue
            nc = min(n, 4)
            nr = (n + nc - 1) // nc

            fig, axes = plt.subplots(nr, nc, figsize=(5*nc, 3.8*nr))
            if n == 1:
                axes = np.array([axes])
            axes = axes.flatten()

            ylim = self._compute_ylim(results, hname)
            t = self.display_time

            for idx, r in enumerate(results):
                ax = axes[idx]
                h = r['heads'].get(hname)
                if h is None:
                    continue
                for key, clr, ls, mk in [
                    ('f4', COLOR_F4, '-.', 's'),
                    ('input', COLOR_F2, '--', 'o'),
                    ('target', COLOR_F1, '-', None),
                    ('richardson_pred', COLOR_RE, ':', 'd'),
                    ('pred', COLOR_ML, '-', 'v'),
                ]:
                    arr = self._sc(h.get(key))
                    if arr is not None and len(arr) == len(t):
                        kw = dict(color=clr, linestyle=ls, linewidth=1.4)
                        if mk:
                            kw['marker'] = mk
                            kw['markevery'] = MARK_EVERY * 2
                            kw['markersize'] = 4
                        ax.plot(t, arr, **kw)
                ax.set_ylim(ylim)
                ax.set_xlim(0, self.display_xmax)
                ax.set_xticks(self.x_ticks)
                mat = r['material'][:3]
                coh = r['cohesion']
                rmse = h['rmse']
                mark = "\u2713" if h['success'] else "\u2717"
                ax.set_title(f"{mark} {mat}_{coh} RMSE={rmse:.5f}", fontsize=10, pad=8)
                ax.grid(True, linestyle=':', alpha=0.4)
                ax.set_xlabel('Time [s]', fontsize=9)

            for idx in range(n, len(axes)):
                axes[idx].set_visible(False)

            hlabel = self.meta.get('head_labels', {}).get(hname, hname)
            plt.suptitle(f'Scale Hierarchy: F4 \u2192 F2 \u2192 F1 | {hlabel}',
                         fontsize=14, fontweight='bold', y=1.02)
            plt.tight_layout()
            _save(fig, self.out / f"hierarchy_{hname}.png")
            print(f"  hierarchy_{hname}.png")

        self._save_curve_legend(self.out / "hierarchy_legend.png")

    # =========================================================================
    # 4. 汇总面板 — RMSE / MAE / SMAPE
    # =========================================================================
    def plot_summary(self, results):
        for hname in self.meta['head_names']:
            tag = f"summary_{hname}"
            print(f"\n[{tag}]")

            n = len(results)
            if n == 0:
                continue

            ml_rmses = [r['heads'][hname]['rmse'] for r in results]
            re_rmses = [r['heads'][hname]['richardson_rmse'] for r in results]
            base_rmses = [r['heads'][hname]['baseline_rmse'] for r in results]

            fig = plt.figure(figsize=(16, 9))
            gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

            mats = [r['material'][:3] for r in results]
            cohs = [r['cohesion'] for r in results]
            labels = [f"{m}_{c}" for m, c in zip(mats, cohs)]

            # Panel 1: RMSE bar (hatched)
            ax1 = fig.add_subplot(gs[0, 0])
            y_pos = np.arange(n)
            ax1.barh(y_pos - 0.2, ml_rmses, 0.2, color=COLOR_ML, alpha=0.8,
                     edgecolor='k', lw=0.5, label='ML')
            ax1.barh(y_pos, re_rmses, 0.2, color=COLOR_RE, alpha=0.8,
                     edgecolor='k', lw=0.5, hatch='//', label='Richardson')
            ax1.barh(y_pos + 0.2, base_rmses, 0.2, color='#7f7f7f', alpha=0.6,
                     edgecolor='k', lw=0.5, hatch='\\\\', label='F2 Baseline')
            ax1.set_yticks(y_pos)
            ax1.set_yticklabels(labels, fontsize=8)
            ax1.set_xlabel(f'RMSE [{self.meta["unit"]}]')
            ax1.set_title('RMSE Comparison')
            ax1.legend(fontsize=7, loc='lower right')
            ax1.grid(axis='x', alpha=0.3, linestyle=':')

            # Panel 2: Boxplot RMSE by cohesion
            ax2 = fig.add_subplot(gs[0, 1])
            cohesions = sorted(set(cohs))
            data_by_coh = [[r['heads'][hname]['rmse']
                            for r in results if r['cohesion'] == cc] for cc in cohesions]
            if any(len(d) > 0 for d in data_by_coh):
                bp = ax2.boxplot(data_by_coh,
                                 labels=[f'SE={c}' for c in cohesions], patch_artist=True)
                for patch in bp['boxes']:
                    patch.set_facecolor('lightblue')
                    patch.set_edgecolor('black')
            ax2.set_ylabel(f'RMSE [{self.meta["unit"]}]')
            ax2.set_title('RMSE by Cohesion')
            ax2.grid(axis='y', alpha=0.3, linestyle=':')

            # Panel 3: Success counts
            ax3 = fig.add_subplot(gs[0, 2])
            n_ml = sum(1 for r in results if r['heads'][hname]['success'])
            n_re = sum(1 for r in results
                       if r['heads'][hname]['richardson_rmse'] < r['heads'][hname]['baseline_rmse'])
            n_wins = sum(1 for r in results
                          if r['heads'][hname]['rmse'] < r['heads'][hname]['richardson_rmse'])
            cats = ['ML<F2', 'RE<F2', 'ML<RE', 'Total']
            vals = [n_ml, n_re, n_wins, n]
            bc = [COLOR_ML, COLOR_RE, '#27ae60', '#95a5a6']
            ax3.bar(cats, vals, color=bc, edgecolor='black', lw=0.5)
            for i, v in enumerate(vals):
                ax3.text(i, v + 0.2, str(v), ha='center', fontsize=10, fontweight='bold')
            ax3.set_ylabel('Count')
            ax3.set_title(f'ML beats F2: {n_ml}/{n}')
            ax3.grid(axis='y', alpha=0.3, linestyle=':')

            # Panel 4: Best
            ax4 = fig.add_subplot(gs[1, 0])
            best = min(results, key=lambda r: r['heads'][hname]['rmse'])
            self._mini_curves(ax4, best, hname, 'Best')

            # Panel 5: Worst
            ax5 = fig.add_subplot(gs[1, 1])
            worst = max(results, key=lambda r: r['heads'][hname]['rmse'])
            self._mini_curves(ax5, worst, hname, 'Worst')

            # Panel 6: Metrics summary table
            ax6 = fig.add_subplot(gs[1, 2])
            ax6.axis('off')
            hlabel = self.meta.get('head_labels', {}).get(hname, hname)
            be = best['heads'][hname]['rmse']
            we = worst['heads'][hname]['rmse']
            txt = (
                f"{self.meta['label']} -- {hlabel}\n"
                f"{'='*30}\n"
                f"Configs:       {n}\n"
                f"ML < F2:       {n_ml}/{n} ({100*n_ml/n:.0f}%)\n"
                f"RE < F2:       {n_re}/{n} ({100*n_re/n:.0f}%)\n"
                f"ML < RE:       {n_wins}/{n}\n"
                f"{'-'*30}\n"
                f"Avg ML  RMSE:  {np.mean(ml_rmses):.6f}\n"
                f"Avg RE  RMSE:  {np.mean(re_rmses):.6f}\n"
                f"Avg F2  RMSE:  {np.mean(base_rmses):.6f}\n"
                f"{'-'*30}\n"
                f"Avg ML  MAE:   {np.mean([r['heads'][hname]['mae'] for r in results]):.6f}\n"
                f"Avg ML  SMAPE: {np.mean([r['heads'][hname]['smape'] for r in results]):.3f}%\n"
                f"{'-'*30}\n"
                f"Best:  {best['material'][:3]}_{best['cohesion']} ({be:.6f})\n"
                f"Worst: {worst['material'][:3]}_{worst['cohesion']} ({we:.6f})"
            )
            ax6.text(0.5, 0.5, txt, transform=ax6.transAxes, fontsize=10,
                     family='monospace', va='center', ha='center',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            plt.suptitle(f'{self.meta["label"]} -- {hlabel}',
                         fontsize=14, fontweight='bold', y=1.02)
            plt.tight_layout()
            _save(fig, self.out / f"summary_{hname}.png")
            print(f"  summary_{hname}.png")

    def _mini_curves(self, ax, r, hname, label):
        h = r['heads'][hname]
        t = self.display_time
        arr = self._sc(h.get('target'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_F1, lw=1.8)
        arr = self._sc(h.get('input'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_F2, lw=1.3, ls='--')
        arr = self._sc(h.get('pred'))
        if arr is not None and len(arr) == len(t):
            ax.plot(t, arr, color=COLOR_ML, lw=1.3)
        rmse = h['rmse']
        ax.set_title(f"{label}: {r['material'][:3]}_{r['cohesion']} (RMSE={rmse:.5f})",
                     fontsize=10, pad=8)
        ax.set_xlim(0, self.display_xmax)
        ax.set_xticks(self.x_ticks)
        ax.grid(True, linestyle=':', alpha=0.4)
        ax.set_xlabel('Time [s]', fontsize=9)

    # =========================================================================
    # 5. CSV 导出
    # =========================================================================
    def export_csv(self, results):
        for hname in self.meta['head_names']:
            tag = f"csv_{hname}"
            print(f"\n[{tag}]  (nth={self.cfg.csv_nth})")
            subdir = self.out / tag
            subdir.mkdir(parents=True, exist_ok=True)

            nth = self.cfg.csv_nth
            t = self.display_time

            for r in results:
                h = r['heads'].get(hname)
                if h is None:
                    continue
                mat, coh = r['material'], r['cohesion']
                df = pd.DataFrame({'Time[s]': t})

                target_c = self._crop(h.get('target'))
                input_c = self._crop(h.get('input'))
                f4_c = self._crop(h.get('f4'))
                pred_c = self._crop(h.get('pred'))
                rich_c = self._crop(h.get('richardson_pred'))

                if target_c is not None and len(target_c) == len(t):
                    df['f1'] = target_c
                if input_c is not None and len(input_c) == len(t):
                    df['f2'] = _decimate(input_c, nth)
                if f4_c is not None and len(f4_c) == len(t):
                    df['f4'] = _decimate(f4_c, nth)
                if pred_c is not None and len(pred_c) == len(t):
                    df['predicted'] = pred_c
                if rich_c is not None and len(rich_c) == len(t):
                    df['richardson'] = _decimate(rich_c, nth)

                fname = f"{self.meta['label'].lower()}_{mat}_{coh}_{hname}.csv"
                df.to_csv(subdir / fname, index=False, float_format='%.8g')
                print(f"  {fname}")

    # =========================================================================
    # 6. 输入消融 ΔRMSE 条状图
    # =========================================================================
    def plot_ablation_bar(self, all_ablation_results: Dict[str, List[Dict]]):
        if 'baseline_full' not in all_ablation_results:
            print("  [WARN] No baseline_full, skipping ablation bar")
            return

        for hname in self.meta['head_names']:
            tag = f"ablation_bar_{hname}"
            print(f"\n[{tag}]")

            base_rmses = [r['heads'][hname]['rmse']
                          for r in all_ablation_results['baseline_full']]
            base_avg = np.mean(base_rmses)

            channels, delta_rmses, delta_stds = [], [], []
            for exp_name, exp_results in all_ablation_results.items():
                if not exp_name.startswith('remove_'):
                    continue
                ch_name = exp_name.replace('remove_', '')
                rmses = [r['heads'][hname]['rmse'] for r in exp_results]
                channels.append(ch_name)
                delta_rmses.append(np.mean(rmses) - base_avg)
                delta_stds.append(np.std(rmses))

            if not channels:
                continue

            order = np.argsort(delta_rmses)[::-1]
            channels = [channels[i] for i in order]
            delta_rmses = [delta_rmses[i] for i in order]
            delta_stds = [delta_stds[i] for i in order]

            en_labels = [CHANNEL_LABELS.get(ch, ch) for ch in channels]

            fig, ax = plt.subplots(figsize=(10, max(5.5, len(channels) * 0.6)))
            colors = ['#e74c3c' if d > 0 else '#27ae60' if d < 0 else '#95a5a6'
                       for d in delta_rmses]
            hatches = ['//' if d > 0 else '\\\\' if d < 0 else '' for d in delta_rmses]

            bars = ax.barh(range(len(en_labels)), delta_rmses,
                           color=colors, edgecolor='k', alpha=0.85,
                           height=0.65, zorder=3)
            for bar, h in zip(bars, hatches):
                bar.set_hatch(h)

            for i, (d, s) in enumerate(zip(delta_rmses, delta_stds)):
                ax.plot([d - s, d + s], [i, i], color='#555555',
                        linewidth=0.8, alpha=0.4, zorder=2)

            ax.set_yticks(range(len(en_labels)))
            ax.set_yticklabels(en_labels)
            ax.axvline(0, color='#333333', lw=0.6, zorder=1)

            hlabel = HEAD_LABELS.get(hname, hname)
            ax.set_title(
                f'Input Channel Ablation -- {hlabel}\n'
                f'(Baseline RMSE = {base_avg:.6f})', pad=15)
            ax.set_xlabel(f'$\\Delta$RMSE [{self.meta["unit"]}] (positive = channel is useful)')
            ax.invert_yaxis()
            _academic_style(ax)

            for i, d in enumerate(delta_rmses):
                pct = d / base_avg * 100 if base_avg != 0 else 0
                offset = max(abs(d) * 0.08, 2e-6)
                x = d + offset if d >= 0 else d - offset
                ha = 'left' if d >= 0 else 'right'
                ax.text(x, i, f'{pct:+.1f}%', va='center', ha=ha,
                        fontsize=10, fontweight='bold',
                        color='#c0392b' if d > 0 else '#1e8449')

            plt.tight_layout()
            _save(fig, self.out / f"ablation_bar_{hname}.png")
            print(f"  ablation_bar_{hname}.png")

            pd.DataFrame({
                'channel': channels,
                'delta_rmse': [f'{d:.8f}' for d in delta_rmses],
                'delta_rmse_pct': [f'{d/base_avg*100:+.2f}%' if base_avg != 0 else '0'
                                   for d in delta_rmses],
            }).to_csv(self.out / f"ablation_bar_{hname}.csv", index=False)

    # =========================================================================
    # 7. 分组消融对比
    # =========================================================================
    def plot_ablation_groups(self, all_ablation_results: Dict[str, List[Dict]]):
        group_keys = [k for k in all_ablation_results
                      if k.startswith('remove_all') or k.startswith('only_')]
        if not group_keys or 'baseline_full' not in all_ablation_results:
            return

        for hname in self.meta['head_names']:
            tag = f"ablation_groups_{hname}"
            print(f"\n[{tag}]")

            base_avg = np.mean([r['heads'][hname]['rmse']
                                for r in all_ablation_results['baseline_full']])

            labels = ['baseline'] + group_keys
            avg_rmses = [base_avg]
            std_rmses = [np.std([r['heads'][hname]['rmse']
                                 for r in all_ablation_results['baseline_full']])]
            for gk in group_keys:
                rmses = [r['heads'][hname]['rmse']
                         for r in all_ablation_results.get(gk, [])]
                avg_rmses.append(np.mean(rmses) if rmses else 0)
                std_rmses.append(np.std(rmses) if rmses else 0)

            en_labels = [GROUP_LABELS.get(lb, lb) for lb in labels]

            fig, ax = plt.subplots(figsize=(max(8, len(labels) * 2.2), 5.5))
            colors = ['#3498db'] + \
                     ['#e74c3c' if e > base_avg else
                      '#27ae60' if e < base_avg else '#95a5a6'
                      for e in avg_rmses[1:]]

            bars = ax.bar(range(len(en_labels)), avg_rmses, color=colors,
                          edgecolor='k', alpha=0.85, width=0.6, zorder=3)

            for i, (avg, std) in enumerate(zip(avg_rmses, std_rmses)):
                ax.plot([i, i], [avg - std, avg + std], color='#555555',
                        linewidth=1.0, alpha=0.4, zorder=4)

            ax.axhline(base_avg, color='#555555', linestyle='--', alpha=0.5,
                       linewidth=1.0, zorder=1)
            ax.set_xticks(range(len(en_labels)))
            ax.set_xticklabels(en_labels, fontsize=11, ha='center')
            ax.set_ylabel(f'Avg RMSE [{self.meta["unit"]}]')

            hlabel = HEAD_LABELS.get(hname, hname)
            ax.set_title(f'Channel Group Ablation -- {hlabel}', pad=15)
            _apply_sci_yaxis(ax)
            _academic_style(ax)

            for i, v in enumerate(avg_rmses):
                pct = (v - base_avg) / base_avg * 100 if base_avg != 0 else 0
                pct_str = f'\n({pct:+.1f}%)' if i > 0 else ''
                ax.text(i, v + max(avg_rmses) * 0.015, f'{v:.5f}{pct_str}',
                        ha='center', fontsize=10, fontweight='bold')

            plt.tight_layout()
            _save(fig, self.out / f"ablation_groups_{hname}.png")
            print(f"  ablation_groups_{hname}.png")

    # =========================================================================
    # 8. 跨目标消融热力图 (ΔRMSE)
    # =========================================================================
    @staticmethod
    def plot_ablation_heatmap(cross_target_data, aux_channel_names, output_path):
        targets = list(cross_target_data.keys())
        target_labels = [HEAD_LABELS.get(t, t) for t in targets]
        channel_labels = [CHANNEL_LABELS.get(ch, ch) for ch in aux_channel_names]
        n_t, n_c = len(targets), len(aux_channel_names)
        matrix = np.zeros((n_t, n_c))
        for i, t in enumerate(targets):
            for j, ch in enumerate(aux_channel_names):
                matrix[i, j] = cross_target_data[t].get(ch, 0.0)

        vmax = max(np.max(np.abs(matrix)) * 1.1, 0.001)

        fig, ax = plt.subplots(figsize=(max(10, n_c*1.2), max(4, n_t*0.9+2)))
        im = ax.imshow(matrix, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
        ax.set_xticks(range(n_c))
        ax.set_xticklabels(channel_labels, rotation=35, ha='right', fontsize=12)
        ax.set_yticks(range(n_t))
        ax.set_yticklabels(target_labels, fontsize=12)

        for i in range(n_t):
            for j in range(n_c):
                v = matrix[i, j]
                color = 'white' if abs(v) > vmax * 0.6 else 'black'
                ax.text(j, i, f'{v:+.5f}', ha='center', va='center',
                        fontsize=9, color=color, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(r'$\Delta$RMSE (positive = channel is useful)', fontsize=12)
        ax.set_title('Input Channel Importance Across Targets', fontsize=14, pad=12)
        plt.tight_layout()
        _save(fig, output_path)
        print(f"  {output_path}")

        csv_path = str(output_path).replace('.png', '.csv')
        rows = [['target/head'] + aux_channel_names]
        for i, t in enumerate(targets):
            rows.append([t] + [f'{matrix[i,j]:.8f}' for j in range(n_c)])
        pd.DataFrame(rows[1:], columns=rows[0]).to_csv(csv_path, index=False)

    # =========================================================================
    # 9. 架构消融 ΔRMSE 条状图
    # =========================================================================
    def plot_arch_ablation_bar(self, all_results: Dict[str, List[Dict]]):
        if 'full_model' not in all_results:
            return

        head_names = self.meta['head_names']
        baseline_rmses = {}
        for hname in head_names:
            baseline_rmses[hname] = np.mean([
                r['heads'][hname]['rmse'] for r in all_results['full_model']])

        variants = [k for k in all_results if k != 'full_model']
        if not variants:
            return

        delta_matrix, std_matrix = {}, {}
        for var in variants:
            delta_matrix[var], std_matrix[var] = {}, {}
            for hname in head_names:
                rmses = [r['heads'][hname]['rmse'] for r in all_results[var]]
                delta_matrix[var][hname] = np.mean(rmses) - baseline_rmses[hname]
                std_matrix[var][hname] = np.std(rmses)

        n_v, n_h = len(variants), len(head_names)
        en_variants = [ARCH_VARIANT_LABELS.get(v, v) for v in variants]

        fig, ax = plt.subplots(figsize=(max(9, 2.2 * n_v + 1), 5.8))
        head_colors = ['#2c3e50', '#e74c3c', '#3498db', '#2ecc71'][:n_h]
        bar_w = 0.65 / n_h
        x_pos = np.arange(n_v)

        for hi, hname in enumerate(head_names):
            offsets = x_pos + (hi - (n_h - 1) / 2) * bar_w
            vals = [delta_matrix[v][hname] for v in variants]
            stds = [std_matrix[v][hname] for v in variants]
            hlabel = HEAD_LABELS.get(hname, hname)
            ax.bar(offsets, vals, bar_w * 0.88,
                   label=hlabel, color=head_colors[hi],
                   edgecolor='k', alpha=0.85, zorder=3)

            for off, v, s in zip(offsets, vals, stds):
                ax.plot([off, off], [v - s, v + s], color='#555555',
                        linewidth=0.8, alpha=0.35, zorder=4)

            for off, v in zip(offsets, vals):
                pct = v / baseline_rmses[hname] * 100 if baseline_rmses[hname] != 0 else 0
                y_off = max(abs(v) * 0.1, 5e-6) * (1 if v >= 0 else -1)
                va = 'bottom' if v >= 0 else 'top'
                ax.text(off, v + y_off, f'{pct:+.1f}%', ha='center', va=va,
                        fontsize=8, fontweight='bold',
                        color='#c0392b' if v > 0 else '#1e8449')

        ax.axhline(0, color='#333333', linewidth=0.6, zorder=1)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(en_variants, fontsize=11, ha='center')
        ax.set_ylabel(f'$\\Delta$RMSE [{self.meta["unit"]}] (vs full model)')
        ax.set_title('Architecture Ablation: Component Impact on RMSE', pad=15)
        _academic_style(ax)

        if n_h > 1:
            ax.legend(fontsize=10, loc='best', framealpha=0.9, edgecolor='#cccccc')

        plt.tight_layout()
        _save(fig, self.out / 'arch_ablation_bar.png')
        print(f"  arch_ablation_bar.png")

        rows_csv = []
        for var in variants:
            row = {'variant': ARCH_VARIANT_LABELS.get(var, var)}
            for hname in head_names:
                row[f'{hname}_delta_rmse'] = f'{delta_matrix[var][hname]:.8f}'
                pct = delta_matrix[var][hname] / baseline_rmses[hname] * 100 if baseline_rmses[hname] != 0 else 0
                row[f'{hname}_delta_pct'] = f'{pct:+.1f}%'
            rows_csv.append(row)
        pd.DataFrame(rows_csv).to_csv(self.out / 'arch_ablation_bar.csv', index=False)

    # =========================================================================
    # 10. 架构消融各配置预测曲线对比
    # =========================================================================
    def plot_arch_ablation_curves(self, all_results: Dict[str, List[Dict]]):
        if 'full_model' not in all_results:
            return

        full_results = all_results['full_model']
        variants = [k for k in all_results if k != 'full_model']
        t = self.display_time

        for hname in self.meta['head_names']:
            tag = f"arch_curves_{hname}"
            subdir = self.out / tag
            subdir.mkdir(parents=True, exist_ok=True)
            print(f"\n[{tag}]")

            ylim = self._compute_ylim(full_results, hname)
            full_by_cfg = {f"{r['material'][:3]}_{r['cohesion']}": r for r in full_results}

            for var_name in variants:
                var_results = all_results[var_name]
                var_by_cfg = {f"{r['material'][:3]}_{r['cohesion']}": r for r in var_results}

                configs = sorted(full_by_cfg.keys())
                n_cfgs = len(configs)
                nc = min(4, n_cfgs)
                nr = (n_cfgs + nc - 1) // nc

                fig, axes = plt.subplots(nr, nc, figsize=(5.2 * nc, 3.8 * nr),
                                         squeeze=False)

                var_color = ARCH_VARIANT_COLORS.get(var_name, '#e74c3c')
                var_label = ARCH_VARIANT_LABELS.get(var_name, var_name)

                for idx, cfg in enumerate(configs):
                    row, col = divmod(idx, nc)
                    ax = axes[row][col]

                    if cfg not in full_by_cfg or cfg not in var_by_cfg:
                        ax.set_visible(False)
                        continue

                    h_full = full_by_cfg[cfg]['heads'][hname]
                    h_var = var_by_cfg[cfg]['heads'][hname]

                    arr = self._sc(h_full.get('target'))
                    if arr is not None and len(arr) == len(t):
                        ax.plot(t, arr, color=COLOR_F1, linewidth=1.8, linestyle='-')
                    arr = self._sc(h_full.get('input'))
                    if arr is not None and len(arr) == len(t):
                        ax.plot(t, arr, color=COLOR_F2, linewidth=1.2, linestyle='--', alpha=0.7)
                    arr = self._sc(h_full.get('pred'))
                    if arr is not None and len(arr) == len(t):
                        ax.plot(t, arr, color=COLOR_ML, linewidth=1.5, linestyle='-')
                    arr = self._sc(h_var.get('pred'))
                    if arr is not None and len(arr) == len(t):
                        ax.plot(t, arr, color=var_color, linewidth=1.5, linestyle='-.')

                    ax.set_xlim(0, self.display_xmax)
                    ax.set_xticks(self.x_ticks)
                    ax.set_ylim(ylim)
                    ax.grid(True, alpha=0.3, linestyle=':')
                    ax.set_xlabel('Time [s]', fontsize=9)
                    if col == 0:
                        ax.set_ylabel(self.meta['ylabel'], fontsize=9)

                    rmse_f = h_full['rmse']
                    rmse_v = h_var['rmse']
                    delta = rmse_v - rmse_f
                    ax.set_title(
                        f"{cfg}  Full={rmse_f:.5f}  "
                        f"{var_label}={rmse_v:.5f}  "
                        f"(\u0394={delta:+.5f})", fontsize=8)

                    if idx == 0:
                        ax.legend(['F1', 'F2', 'Full', var_label],
                                  fontsize=6.5, loc='lower right', framealpha=0.7)

                for idx in range(n_cfgs, nr * nc):
                    row, col = divmod(idx, nc)
                    axes[row][col].set_visible(False)

                hlabel = HEAD_LABELS.get(hname, hname)
                fig.suptitle(f'Architecture Ablation: {var_label} -- {hlabel}',
                             fontsize=13, fontweight='bold', y=1.02)
                plt.tight_layout()
                _save(fig, subdir / f"curves_{var_name}.png")
                print(f"  curves_{var_name}.png")

    # =========================================================================
    # 11. 架构消融参数量 vs RMSE 散点图
    # =========================================================================
    def plot_arch_ablation_params(self, all_results: Dict[str, List[Dict]],
                                  param_counts: Dict[str, int]):
        if not param_counts:
            return

        head_names = self.meta['head_names']
        n_h = len(head_names)

        fig, axes = plt.subplots(1, n_h, figsize=(5.8 * n_h, 5), squeeze=False)
        axes = axes[0]

        for hi, hname in enumerate(head_names):
            ax = axes[hi]
            for exp_name, results in all_results.items():
                if exp_name not in param_counts:
                    continue
                n_p = param_counts[exp_name]
                rmses = [r['heads'][hname]['rmse'] for r in results]
                avg_rmse = np.mean(rmses)
                std_rmse = np.std(rmses)
                color = ARCH_VARIANT_COLORS.get(exp_name, '#7f8c8d')
                label = ARCH_VARIANT_LABELS.get(exp_name, exp_name)

                if exp_name == 'full_model':
                    ax.scatter(n_p, avg_rmse, s=160, c=color, marker='*',
                               zorder=5, edgecolors='black', linewidths=0.8, label=label)
                else:
                    ax.scatter(n_p, avg_rmse, s=80, c=color, marker='o',
                               zorder=4, edgecolors='white', linewidths=0.5, label=label)
                ax.errorbar(n_p, avg_rmse, yerr=std_rmse, fmt='none',
                            ecolor=color, alpha=0.4, capsize=3, linewidth=0.8)

            hlabel = HEAD_LABELS.get(hname, hname)
            ax.set_xlabel('Model Parameters')
            ax.set_ylabel(f'Avg RMSE [{self.meta["unit"]}]')
            ax.set_title(hlabel, fontweight='bold')
            ax.legend(fontsize=8, loc='best', framealpha=0.9)
            _academic_style(ax)
            ax.ticklabel_format(axis='x', style='sci', scilimits=(3, 5))

        fig.suptitle('Architecture Efficiency: Parameters vs. RMSE',
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        _save(fig, self.out / 'arch_ablation_params.png')
        print(f"  arch_ablation_params.png")

    # =========================================================================
    # 12. 架构消融 CSV 汇总导出
    # =========================================================================
    def _export_arch_ablation_csv(self, all_results: Dict[str, List[Dict]],
                                   param_counts: Dict[str, int]):
        head_names = self.meta['head_names']
        baseline_rmses = {}
        for hname in head_names:
            baseline_rmses[hname] = np.mean([
                r['heads'][hname]['rmse'] for r in all_results.get('full_model', [])])

        rows = []
        for exp_name in ['full_model'] + [k for k in all_results if k != 'full_model']:
            if exp_name not in all_results:
                continue
            results = all_results[exp_name]
            row = {'variant': ARCH_VARIANT_LABELS.get(exp_name, exp_name),
                   'params': param_counts.get(exp_name, '')}
            for hname in head_names:
                rmses = [r['heads'][hname]['rmse'] for r in results]
                maes = [r['heads'][hname]['mae'] for r in results]
                smapes = [r['heads'][hname]['smape'] for r in results]
                avg = np.mean(rmses)
                delta = avg - baseline_rmses.get(hname, avg)
                row[f'{hname}_avg_rmse'] = f'{avg:.6f}'
                row[f'{hname}_delta_rmse'] = (f'{delta:+.6f}' if exp_name != 'full_model' else '-')
                row[f'{hname}_avg_mae'] = f'{np.mean(maes):.6f}'
                row[f'{hname}_avg_smape'] = f'{np.mean(smapes):.4f}'
            rows.append(row)
        pd.DataFrame(rows).to_csv(self.out / 'arch_ablation_summary.csv', index=False)
        print(f"  arch_ablation_summary.csv")

    # =========================================================================
    # 架构消融主入口
    # =========================================================================
    def run_arch_ablation(self, all_results: Dict[str, List[Dict]],
                          param_counts: Optional[Dict[str, int]] = None):
        print("=" * 70)
        print(f"V14 Architecture Ablation Visualizer -- {self.meta['label']}")
        print("=" * 70)

        self.plot_arch_ablation_bar(all_results)
        self.plot_arch_ablation_curves(all_results)
        if param_counts:
            self.plot_arch_ablation_params(all_results, param_counts)
        self._export_arch_ablation_csv(all_results, param_counts or {})

        print(f"\nArch Ablation Done! -> {self.out}")

    # =========================================================================
    # 主入口
    # =========================================================================
    def run(self, input_dir: str):
        print("=" * 70)
        print(f"V14 Visualizer -- {self.meta['label']}")
        print(f"Display: 0 -> {self.display_xmax:.1f}s (steady state)")
        print("=" * 70)

        results = self.load(input_dir)
        if not results:
            return

        self.plot_predictions(results)
        self.plot_by_material(results)
        self.plot_hierarchy(results)
        self.plot_summary(results)
        self.export_csv(results)

        print(f"\nDone! -> {self.out}")

    def run_ablation(self, all_ablation_results: Dict[str, List[Dict]]):
        print("=" * 70)
        print(f"V14 Ablation Visualizer -- {self.meta['label']}")
        print("=" * 70)

        self.plot_ablation_bar(all_ablation_results)
        self.plot_ablation_groups(all_ablation_results)

        print(f"\nAblation Done! -> {self.out}")


# =============================================================================
# CLI
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description='V14 Unified Visualizer')
    parser.add_argument('--target', required=True, choices=['ke', 'ke_comp', 'gt'])
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--output_dir', default=None)
    parser.add_argument('--nth', type=int, default=15, help='CSV decimation interval')
    args = parser.parse_args()

    out = args.output_dir or str(Path(args.input_dir) / 'viz_v14')
    cfg = VizConfig(csv_nth=args.nth)
    VisualizerV14(args.target, out, cfg).run(args.input_dir)


if __name__ == "__main__":
    main()
