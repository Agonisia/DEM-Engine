#!/usr/bin/env python3
"""
V14 SE020-Focused Seed Sweep
=============================
专门针对 SE020 (cohesion='020') 配置的种子搜索.
目标: 找到让 ML prediction 在 SE020 条件下优于 F2 baseline 的种子.

用法:
  python seed_sweep_se020_v14.py --target ke --seed_start 0 --seed_end 199
  python seed_sweep_se020_v14.py --target ke --seeds 0,42,100
  python seed_sweep_se020_v14.py --target gt --seed_start 0 --seed_end 99 --top_k 20
"""

import numpy as np
import json
import time
import argparse
from pathlib import Path
from torch.utils.data import DataLoader

from run_predictor_v14 import (
    set_seed, TARGET_REGISTRY, load_all_data,
    UnifiedPreprocessor, UnifiedPairBuilder,
    UnifiedTrainDataset, UnifiedInferDataset,
    UnifiedTrainer, UnifiedEvaluator,
    MATERIALS,
)
from unified_predictor_v14 import (
    PredictorV14, make_ke_config, make_ke_comp_config, make_gt_config,
)


CONFIG_FACTORY = {
    'ke': make_ke_config,
    'ke_comp': make_ke_comp_config,
    'gt': make_gt_config,
}

# SE020 的 3 个配置 key (material[:3]_cohesion)
SE020_CONFIGS = [f"{mat[:3]}_020" for mat in MATERIALS]  # sta_020, Dua_020, Siz_020


def run_one_seed(seed, spec, config, all_data, prep, builder):
    """单次 seed 训练 + 评估, 返回以 SE020 为重点的指标."""
    set_seed(seed)

    train_samples = builder.build_training_samples(all_data)
    infer_samples = builder.build_inference_samples(all_data)
    train_ds = UnifiedTrainDataset(train_samples, spec)
    infer_ds = UnifiedInferDataset(infer_samples, spec)
    infer_loader = DataLoader(infer_ds, batch_size=config.batch_size, shuffle=False)

    model = PredictorV14(config)
    trainer = UnifiedTrainer(config, spec, model)
    trainer.train(train_ds)

    evaluator = UnifiedEvaluator(config, spec, model, prep)
    results = evaluator.evaluate(infer_loader)

    metrics = {'seed': seed}

    for hname in spec.head_names:
        # --- SE020 专项指标 ---
        se020_rmses, se020_successes, se020_details = [], [], []
        # --- 全局指标 ---
        all_rmses, all_successes = [], []

        for r in results:
            mat3 = r['material'][:3]
            coh = r['cohesion']
            cfg_key = f"{mat3}_{coh}"
            h = r['heads'][hname]

            all_rmses.append(h['rmse'])
            all_successes.append(h['success'])

            metrics[f'{hname}_{cfg_key}_rmse'] = h['rmse']
            metrics[f'{hname}_{cfg_key}_baseline'] = h['baseline_rmse']
            metrics[f'{hname}_{cfg_key}_richardson'] = h['richardson_rmse']
            metrics[f'{hname}_{cfg_key}_success'] = h['success']

            if coh == '020':
                se020_rmses.append(h['rmse'])
                se020_successes.append(h['success'])
                se020_details.append({
                    'config': cfg_key,
                    'ml_rmse': h['rmse'],
                    'f2_rmse': h['baseline_rmse'],
                    're_rmse': h['richardson_rmse'],
                    'success': h['success'],
                    'ml_vs_f2_pct': (h['rmse'] - h['baseline_rmse']) / h['baseline_rmse'] * 100,
                })

        # SE020 汇总
        metrics[f'{hname}_se020_avg_rmse'] = float(np.mean(se020_rmses))
        metrics[f'{hname}_se020_n_success'] = sum(se020_successes)
        metrics[f'{hname}_se020_all_success'] = all(se020_successes)
        metrics[f'{hname}_se020_details'] = se020_details

        # 非 SE020 的表现 (监控退化)
        non020_rmses = [r['heads'][hname]['rmse'] for r in results if r['cohesion'] != '020']
        non020_success = [r['heads'][hname]['success'] for r in results if r['cohesion'] != '020']
        metrics[f'{hname}_non020_avg_rmse'] = float(np.mean(non020_rmses))
        metrics[f'{hname}_non020_n_success'] = sum(non020_success)

        # 全局
        metrics[f'{hname}_global_avg_rmse'] = float(np.mean(all_rmses))
        metrics[f'{hname}_global_n_success'] = sum(all_successes)

    # 跨 head 的 SE020 总平均 (用于排序)
    se020_avg_all = np.mean([metrics[f'{h}_se020_avg_rmse'] for h in spec.head_names])
    metrics['se020_overall_avg_rmse'] = float(se020_avg_all)
    metrics['final_train_loss'] = trainer.history['total'][-1]

    return metrics


def print_seed_result(i, n_total, metrics, spec):
    """单个 seed 的实时打印, 高亮 SE020 结果."""
    seed = metrics['seed']
    elapsed_str = ""

    parts = [f"[{i+1:3d}/{n_total}] seed={seed:5d}"]

    for hname in spec.head_names:
        se020_n = metrics[f'{hname}_se020_n_success']
        se020_avg = metrics[f'{hname}_se020_avg_rmse']
        non020_n = metrics[f'{hname}_non020_n_success']

        # 高亮: SE020 有 success 的用 ★, 全部 success 用 ★★★
        if metrics[f'{hname}_se020_all_success']:
            flag = " ★★★ ALL SE020 SUCCESS!"
        elif se020_n > 0:
            flag = f" ★ SE020 {se020_n}/3 win"
        else:
            flag = ""

        parts.append(f"SE020={se020_avg:.6f}({se020_n}/3) "
                      f"other={metrics[f'{hname}_non020_avg_rmse']:.6f}({non020_n}/12)"
                      f"{flag}")

    print(f"  {'  |  '.join(parts)}")

    # SE020 有任意 success 时, 打印明细
    for hname in spec.head_names:
        if metrics[f'{hname}_se020_n_success'] > 0:
            for d in metrics[f'{hname}_se020_details']:
                mark = "✓" if d['success'] else "✗"
                print(f"      {mark} {d['config']}: ML={d['ml_rmse']:.6f} "
                      f"F2={d['f2_rmse']:.6f} ({d['ml_vs_f2_pct']:+.1f}%)")


def main():
    parser = argparse.ArgumentParser(description='V14 SE020-Focused Seed Sweep')
    parser.add_argument('--target', required=True, choices=['ke', 'ke_comp', 'gt'])
    parser.add_argument('--seed_start', type=int, default=0)
    parser.add_argument('--seed_end', type=int, default=199)
    parser.add_argument('--seeds', type=str, default=None,
                        help='指定 seed 列表, 逗号分隔 (覆盖 start/end)')
    parser.add_argument('--top_k', type=int, default=15)
    parser.add_argument('--sta_dir', type=str, default='sta_results')
    parser.add_argument('--ext_dir', type=str, default='extended_features_results')
    parser.add_argument('--output_dir', type=str, default=None)
    args = parser.parse_args()

    spec = TARGET_REGISTRY[args.target]
    config = CONFIG_FACTORY[args.target]()
    config.validate()

    if args.seeds:
        seeds = [int(s.strip()) for s in args.seeds.split(',')]
    else:
        seeds = list(range(args.seed_start, args.seed_end + 1))

    output_dir = Path(args.output_dir or f'{spec.output_dir}/seed_sweep_se020')
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"V14 SE020-Focused Seed Sweep — {spec.name}")
    print(f"Goal: Find seeds where ML beats F2 baseline on SE020")
    print(f"Seeds: {len(seeds)} ({seeds[0]} ~ {seeds[-1]})")
    print(f"SE020 configs: {SE020_CONFIGS}")
    print(f"Output: {output_dir}")
    print("=" * 70)

    # 数据加载 + 预处理只做一次
    print("\n[1] Loading data (once)...")
    all_data = load_all_data(spec, args.sta_dir, args.ext_dir)
    if not all_data:
        print("No data found!")
        return

    print("\n[2] Preprocessing (once)...")
    prep = UnifiedPreprocessor(config, spec)
    prep.fit(all_data)
    for hname in spec.head_names:
        config.direct_res_scale = prep.get_target_std(hname)

    builder = UnifiedPairBuilder(config, spec, prep)

    # Sweep
    print(f"\n[3] Sweeping {len(seeds)} seeds (focus: SE020)...\n")
    all_metrics = []
    winners = []  # seeds where any SE020 config succeeds
    t_start = time.time()

    for i, seed in enumerate(seeds):
        t0 = time.time()
        metrics = run_one_seed(seed, spec, config, all_data, prep, builder)
        elapsed = time.time() - t0
        metrics['elapsed'] = elapsed
        all_metrics.append(metrics)

        print_seed_result(i, len(seeds), metrics, spec)

        # 记录有 SE020 success 的 seed
        any_se020_win = any(
            metrics[f'{h}_se020_n_success'] > 0 for h in spec.head_names)
        if any_se020_win:
            winners.append(metrics)

        # 增量保存
        save_metrics = []
        for m in all_metrics:
            m_save = {k: v for k, v in m.items()
                      if not k.endswith('_details')}
            m_save['se020_details'] = {
                hname: m[f'{hname}_se020_details']
                for hname in spec.head_names
            }
            save_metrics.append(m_save)
        with open(output_dir / 'sweep_se020_results.json', 'w') as f:
            json.dump(save_metrics, f, indent=2)

    total_time = time.time() - t_start

    # === 按 SE020 avg RMSE 排序 ===
    all_metrics.sort(key=lambda m: m['se020_overall_avg_rmse'])
    top_k = min(args.top_k, len(all_metrics))

    # === 结果打印 ===
    print(f"\n{'=' * 80}")
    print(f"SE020-Focused Seed Sweep Complete — {spec.name}")
    print(f"Total: {len(seeds)} seeds in {total_time:.0f}s "
          f"({total_time/len(seeds):.1f}s/seed)")
    print(f"{'=' * 80}")

    # 汇总: 有多少 seed 在 SE020 上有 success
    n_any_win = len(winners)
    n_all_win = sum(1 for m in all_metrics
                     if all(m[f'{h}_se020_all_success'] for h in spec.head_names))
    print(f"\n  SE020 Success Summary:")
    print(f"    Seeds with any SE020 win:  {n_any_win}/{len(seeds)}")
    print(f"    Seeds with ALL SE020 win:  {n_all_win}/{len(seeds)}")

    # Top-K 表格
    print(f"\n  Top {top_k} Seeds (sorted by SE020 avg RMSE):")
    print()

    for hname in spec.head_names:
        print(f"  ── {hname} ──")
        header = (f"  {'Rank':<5s} {'Seed':<7s} "
                  f"{'SE020_avg':>11s} {'SE020_win':>9s} ")
        for cfg in SE020_CONFIGS:
            print_cfg = cfg.replace('_020', '')
            header += f" {print_cfg+'_ML':>10s} {print_cfg+'_F2':>10s} {'gap%':>6s}"
        header += f" {'Other_avg':>11s} {'Other_win':>9s} {'Global':>11s}"
        print(header)
        print(f"  {'-' * (len(header) - 2)}")

        for rank, m in enumerate(all_metrics[:top_k]):
            se020_n = m[f'{hname}_se020_n_success']
            if m[f'{hname}_se020_all_success']:
                flag = " ★★★"
            elif se020_n > 0:
                flag = " ★"
            else:
                flag = ""

            line = (f"  {rank+1:<5d} {m['seed']:<7d} "
                    f"{m[f'{hname}_se020_avg_rmse']:>11.6f} "
                    f"{se020_n:>5d}/3{flag:>3s} ")

            for cfg in SE020_CONFIGS:
                ml = m.get(f'{hname}_{cfg}_rmse', float('nan'))
                f2 = m.get(f'{hname}_{cfg}_baseline', float('nan'))
                gap = (ml - f2) / f2 * 100 if f2 > 0 else 0
                line += f" {ml:>10.6f} {f2:>10.6f} {gap:>+5.1f}%"

            line += (f" {m[f'{hname}_non020_avg_rmse']:>11.6f} "
                     f"{m[f'{hname}_non020_n_success']:>5d}/12  "
                     f"{m[f'{hname}_global_avg_rmse']:>11.6f}")
            print(line)

    # Best / Worst
    best = all_metrics[0]
    worst = all_metrics[-1]
    print(f"\n  Best SE020:  seed={best['seed']}  "
          f"SE020_avg={best['se020_overall_avg_rmse']:.6f}  "
          f"global_avg={best[f'{spec.head_names[0]}_global_avg_rmse']:.6f}")
    print(f"  Worst SE020: seed={worst['seed']}  "
          f"SE020_avg={worst['se020_overall_avg_rmse']:.6f}  "
          f"global_avg={worst[f'{spec.head_names[0]}_global_avg_rmse']:.6f}")

    spread = (worst['se020_overall_avg_rmse'] - best['se020_overall_avg_rmse'])
    spread_pct = spread / worst['se020_overall_avg_rmse'] * 100
    print(f"  SE020 RMSE spread: {spread:.6f} ({spread_pct:.1f}%)")

    # Winners 明细
    if winners:
        print(f"\n  {'=' * 60}")
        print(f"  SEEDS WITH SE020 SUCCESS (ML < F2 baseline):")
        print(f"  {'=' * 60}")
        winners.sort(key=lambda m: -sum(
            m[f'{h}_se020_n_success'] for h in spec.head_names))
        for m in winners:
            total_wins = sum(m[f'{h}_se020_n_success'] for h in spec.head_names)
            print(f"\n  seed={m['seed']} ({total_wins} SE020 wins):")
            for hname in spec.head_names:
                for d in m[f'{hname}_se020_details']:
                    mark = "✓ WIN" if d['success'] else "✗ lose"
                    print(f"    [{hname}] {d['config']}: {mark}  "
                          f"ML={d['ml_rmse']:.6f}  F2={d['f2_rmse']:.6f}  "
                          f"({d['ml_vs_f2_pct']:+.1f}%)")
    else:
        print(f"\n  No seeds found where ML beats F2 on any SE020 config.")
        print(f"  Consider: expanding seed range, or tuning model hyperparameters.")

    # === 保存最终结果 ===
    save_final = []
    for m in all_metrics:
        m_save = {k: v for k, v in m.items() if not k.endswith('_details')}
        m_save['se020_details'] = {
            hname: m[f'{hname}_se020_details'] for hname in spec.head_names
        }
        save_final.append(m_save)
    with open(output_dir / 'sweep_se020_results_sorted.json', 'w') as f:
        json.dump(save_final, f, indent=2)

    # CSV
    import pandas as pd
    rows = []
    for m in all_metrics:
        row = {
            'seed': m['seed'],
            'se020_overall_avg_rmse': m['se020_overall_avg_rmse'],
            'train_loss': m['final_train_loss'],
        }
        for hname in spec.head_names:
            row[f'{hname}_se020_avg_rmse'] = m[f'{hname}_se020_avg_rmse']
            row[f'{hname}_se020_n_success'] = m[f'{hname}_se020_n_success']
            row[f'{hname}_non020_avg_rmse'] = m[f'{hname}_non020_avg_rmse']
            row[f'{hname}_non020_n_success'] = m[f'{hname}_non020_n_success']
            row[f'{hname}_global_avg_rmse'] = m[f'{hname}_global_avg_rmse']
            for cfg in SE020_CONFIGS:
                row[f'{hname}_{cfg}_rmse'] = m.get(f'{hname}_{cfg}_rmse', '')
                row[f'{hname}_{cfg}_baseline'] = m.get(f'{hname}_{cfg}_baseline', '')
        rows.append(row)
    pd.DataFrame(rows).to_csv(output_dir / 'sweep_se020_summary.csv',
                               index=False, float_format='%.8f')

    print(f"\n  Results saved to {output_dir}/")
    print(f"    sweep_se020_results_sorted.json")
    print(f"    sweep_se020_summary.csv")


if __name__ == "__main__":
    main()
