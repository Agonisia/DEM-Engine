#!/usr/bin/env python3
"""
=============================================================================
 Dataset 2 Generator: Material Type (Standard vs SizeDiff)
 Optimized: reuse existing Standard data, only extract SizeDiff from raw
=============================================================================
"""

import os
import glob
import re
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm

# ===========================================================================
# Configuration
# ===========================================================================
EXISTING_STANDARD_CSV = "Assignment2_PDF_Dataset.csv"  # Already extracted
SIZEDIFF_BASE_PATH = "/media/huyuze/Fanxiang1/DEME_Data/backup251115"
SIZEDIFF_PATTERN = r'SUPMixerOutput_SizeDiff_f(\d+)se(\d+)$'

GLOBAL_BINS = 50
V_RANGE = (0, 2.0)
BIN_EDGES = np.linspace(V_RANGE[0], V_RANGE[1], GLOBAL_BINS + 1)

FEATURE_NOISE_RATIO = 0.10
LABEL_FLIP_RATIO = 0.10
RANDOM_SEED = 42

OUTPUT_FINAL = "Assignment2_MaterialType_Dataset.csv"


# ===========================================================================
# Feature Extraction (SizeDiff only)
# ===========================================================================
def extract_features_from_frame(file_path):
    """Read a single DEM output CSV and compute the velocity PDF histogram."""
    try:
        df = pd.read_csv(file_path, usecols=['v_x', 'v_y', 'v_z'])
        v_total = np.sqrt(df['v_x']**2 + df['v_y']**2 + df['v_z']**2)
        hist, _ = np.histogram(v_total, bins=BIN_EDGES, density=True)
        return hist
    except Exception as e:
        print(f"    [ERROR] {file_path}: {e}")
        return None


def extract_sizediff(base_path):
    """Extract SizeDiff samples from raw DEM files."""
    base_path = Path(base_path)
    samples = []

    if not base_path.exists():
        print(f"  [ERROR] Path does not exist: {base_path}")
        return samples

    print(f"  Scanning: {base_path}")

    for folder_path in sorted(base_path.iterdir()):
        if not folder_path.is_dir():
            continue

        match = re.match(SIZEDIFF_PATTERN, folder_path.name)
        if not match:
            continue

        scale_factor = int(match.group(1))
        surface_energy = match.group(2)

        # Skip se005 — not used in this assignment
        if surface_energy == "005":
            continue

        csv_files = sorted(glob.glob(str(folder_path / 'mixer_output_*.csv')))
        if not csv_files:
            continue

        # Extract frame numbers, take last 50% (steady state)
        frames = []
        for f in csv_files:
            m = re.search(r'_(\d+)\.csv', f)
            if m:
                frames.append((f, int(m.group(1))))
        frames.sort(key=lambda x: x[1])

        cutoff = int(len(frames) * 0.5)
        steady_frames = frames[cutoff:]

        print(f"    {folder_path.name} -> {len(steady_frames)} frames")

        for file_path, frame_num in tqdm(steady_frames,
                                          desc="      Extracting", leave=False):
            features = extract_features_from_frame(file_path)
            if features is not None:
                sample = {
                    'material_type': 'SizeDiff',
                    'surface_energy': surface_energy,
                    'scale_factor': scale_factor,
                    'frame': frame_num,
                    'source_folder': folder_path.name,
                }
                for i, val in enumerate(features):
                    sample[f'bin_{i}'] = val
                samples.append(sample)

    return samples


# ===========================================================================
# Main
# ===========================================================================
def main():
    print("=" * 60)
    print(" Dataset 2: Standard vs SizeDiff")
    print("=" * 60)

    feat_cols = [f'bin_{i}' for i in range(GLOBAL_BINS)]

    # --- Step 1: Load existing Standard data ---
    print(f"\n[1/4] Loading existing Standard data: {EXISTING_STANDARD_CSV}")
    df_std = pd.read_csv(EXISTING_STANDARD_CSV)
    # Keep only the columns we need
    keep_cols = ['material_type', 'surface_energy', 'scale_factor',
                 'frame', 'source_folder'] + feat_cols
    # Make sure material_type is set
    df_std['material_type'] = 'Standard'
    df_std = df_std[[c for c in keep_cols if c in df_std.columns]]
    print(f"  Standard samples: {len(df_std)}")

    # --- Step 2: Extract SizeDiff ---
    print(f"\n[2/4] Extracting SizeDiff from raw DEM files...")
    sizediff_samples = extract_sizediff(SIZEDIFF_BASE_PATH)

    if not sizediff_samples:
        print("\n  [ERROR] No SizeDiff samples found! Check path and pattern.")
        return

    df_sd = pd.DataFrame(sizediff_samples)
    print(f"\n  SizeDiff samples: {len(df_sd)}")

    # --- Step 3: Combine and balance ---
    print(f"\n[3/4] Combining and balancing...")
    df_all = pd.concat([df_std, df_sd], ignore_index=True)
    print(f"  Combined total: {len(df_all)}")
    print(f"  Distribution before balancing:")
    print(df_all['material_type'].value_counts().to_string())

    # Downsample to the smaller class
    min_count = df_all['material_type'].value_counts().min()
    rng = np.random.RandomState(RANDOM_SEED)

    balanced = []
    for cls in df_all['material_type'].unique():
        subset = df_all[df_all['material_type'] == cls]
        if len(subset) > min_count:
            subset = subset.sample(n=min_count, random_state=RANDOM_SEED)
        balanced.append(subset)
        print(f"    {cls}: -> {len(subset)} samples")

    df = pd.concat(balanced, ignore_index=True).sample(
        frac=1, random_state=RANDOM_SEED
    ).reset_index(drop=True)

    # --- Step 4: Inject noise ---
    print(f"\n[4/4] Injecting noise...")

    # Feature noise
    X = df[feat_cols].values.copy()
    for j in range(X.shape[1]):
        col_std = np.std(X[:, j])
        if col_std > 0:
            X[:, j] += rng.normal(0, col_std * FEATURE_NOISE_RATIO, size=X.shape[0])
    df[feat_cols] = np.clip(X, 0, None)

    # Label noise
    labels = df['material_type'].astype(str).values.copy()
    unique_labels = sorted(set(labels))
    n_flip = int(len(labels) * LABEL_FLIP_RATIO)
    flip_idx = rng.choice(len(labels), size=n_flip, replace=False)

    flipped = 0
    for idx in flip_idx:
        candidates = [l for l in unique_labels if l != labels[idx]]
        if candidates:
            labels[idx] = rng.choice(candidates)
            flipped += 1
    df['material_type'] = labels

    # --- Save ---
    df.to_csv(OUTPUT_FINAL, index=False)

    print(f"\n  Final dataset: {len(df)} samples")
    print(f"  Classes: {unique_labels}")
    print(f"  Distribution:")
    print(df['material_type'].value_counts().to_string())
    print(f"  Flipped: {flipped} ({flipped/len(df)*100:.1f}%)")
    print(f"  Bayes ceiling: ~{(1 - flipped/len(df))*100:.0f}%")
    print(f"\n  Saved: {OUTPUT_FINAL}")
    print("=" * 60)


if __name__ == "__main__":
    main()