#!/usr/bin/env python3
"""
=============================================================================
 Mid-Term Dataset Generator: Cross-Material Surface Energy Classification
 INFO6105 Spring 2026 | Peize Li
=============================================================================

Purpose:
    Generate a NEW dataset for the mid-term project that is DISTINCT from
    Assignment 2.  This dataset combines ALL THREE material types
    (Standard, SizeDiff, DualDensity) for 5-class surface energy
    classification — a harder, more realistic inverse problem than the
    single-material dataset used in Assignment 2.

What's different from Assignment 2:
    - Assignment 2 Dataset 1: Standard-only, SE 5-class
    - Assignment 2 Dataset 2: Standard vs SizeDiff, material binary
    - THIS dataset: Standard + SizeDiff + DualDensity, SE 5-class
      → The model must learn to predict surface energy REGARDLESS of
        material composition, which is the real industrial challenge.

Usage:
    1. First, make sure Assignment2_PDF_Dataset.csv contains all 3 material
       types.  If it only has Standard, re-run sta_xx.py on a directory
       that includes all three material folders.
    2. Run this script:  python generate_midterm_dataset.py
    3. Output: Midterm_CrossMaterial_SE_Dataset.csv

Noise injection:
    - Feature noise: 15% additive Gaussian (per-column std-scaled)
    - Label flip: 20% uniform random across 5 SE classes
    - Theoretical Bayes ceiling: ~84%
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ===================== CONFIGURATION =====================
# Input: the raw PDF dataset extracted by sta_xx.py
# This should contain columns: material_type, surface_energy, bin_0..bin_49
INPUT_CSV = "Midterm_Raw_PDF_Dataset.csv"

# Output: the new mid-term dataset
OUTPUT_CSV = "Midterm_CrossMaterial_SE_Dataset.csv"

# Noise parameters (same severity as Assignment 2 for comparable analysis)
FEATURE_NOISE_RATIO = 0.15   # 15% of per-column std
LABEL_FLIP_RATIO = 0.20      # 20% uniform random label flip
RANDOM_SEED = 2026            # Different seed from Assignment 2 (which used default)

# =========================================================

def main():
    rng = np.random.RandomState(RANDOM_SEED)

    print("=" * 60)
    print(" Mid-Term Dataset Generator")
    print(" Cross-Material Surface Energy Classification")
    print("=" * 60)

    # --- Step 1: Load raw data ---
    print(f"\n[1/5] Loading raw data: {INPUT_CSV}")
    if not Path(INPUT_CSV).exists():
        print(f"  ERROR: {INPUT_CSV} not found!")
        print(f"  Please run sta_xx_multi.py first to extract features from DEM data.")
        return

    df = pd.read_csv(INPUT_CSV)
    feat_cols = [c for c in df.columns if c.startswith('bin_')]
    print(f"  Total samples: {len(df)}")
    print(f"  Features: {len(feat_cols)} velocity bins")

    # --- Step 2: Check material types ---
    print(f"\n[2/5] Checking material composition...")
    if 'material_type' not in df.columns:
        print("  WARNING: 'material_type' column not found.")
        print("  Assuming all samples are Standard. This won't work for the midterm!")
        print("  Please re-run sta_xx.py on a directory with all 3 material types.")
        return

    mat_counts = df['material_type'].value_counts()
    print(f"  Material distribution:")
    for mat, cnt in mat_counts.items():
        print(f"    {mat}: {cnt} samples")

    required_types = {'Standard', 'SizeDiff', 'DualDensity'}
    found_types = set(df['material_type'].unique())
    missing = required_types - found_types
    if missing:
        print(f"\n  ERROR: Missing material types: {missing}")
        print(f"  Please re-run sta_xx.py on a directory that includes all 3 types.")
        return

    print(f"  All 3 material types found.")

    # --- Step 3: Check and clean surface energy labels ---
    print(f"\n[3/5] Checking surface energy labels...")
    se_counts = df['surface_energy'].value_counts().sort_index()
    print(f"  SE distribution (before noise):")
    for se, cnt in se_counts.items():
        print(f"    SE-{se}: {cnt} samples")

    # Filter out any SE values we don't want (e.g., se005 if present)
    # Keep only the 5 standard levels
    valid_se = ['000', '010', '020', '040']
    # Handle both string and numeric SE values
    df['surface_energy'] = df['surface_energy'].astype(str).str.strip()
    
    # Check what format the SE values are in
    se_unique = df['surface_energy'].unique()
    print(f"  Unique SE values found: {se_unique}")
    
    # Try to map numeric values to string codes if needed
    se_map = {'0': '000', '10': '010', '20': '020', '40': '040', '80': '080',
              '0.0': '000', '10.0': '010', '20.0': '020', '40.0': '040', '80.0': '080'}
    df['surface_energy'] = df['surface_energy'].map(lambda x: se_map.get(x, x))

    before_filter = len(df)
    df = df[df['surface_energy'].isin(valid_se)].reset_index(drop=True)
    if len(df) < before_filter:
        print(f"  Filtered out {before_filter - len(df)} samples with non-standard SE values")

    print(f"  Final sample count: {len(df)}")

    # --- Step 4: Inject noise ---
    print(f"\n[4/5] Injecting noise (seed={RANDOM_SEED})...")

    # 4a. Feature noise: per-column Gaussian scaled by column std
    print(f"  Feature noise: {FEATURE_NOISE_RATIO*100:.0f}% of per-column std...")
    X = df[feat_cols].values.copy()
    for j in range(X.shape[1]):
        col_std = np.std(X[:, j])
        if col_std > 0:
            X[:, j] += rng.normal(0, col_std * FEATURE_NOISE_RATIO, size=X.shape[0])
    # Clip to non-negative (PDF values can't be negative)
    df[feat_cols] = np.clip(X, 0, None)

    # 4b. Label noise: uniform random flip
    labels = df['surface_energy'].values.copy()
    all_labels = sorted(df['surface_energy'].unique())
    n_flip = int(len(labels) * LABEL_FLIP_RATIO)
    flip_idx = rng.choice(len(labels), size=n_flip, replace=False)

    actually_flipped = 0
    for idx in flip_idx:
        # Pick a DIFFERENT label (guaranteed flip, not maybe-same)
        candidates = [l for l in all_labels if l != labels[idx]]
        if candidates:
            labels[idx] = rng.choice(candidates)
            actually_flipped += 1

    df['surface_energy'] = labels

    effective_flip = actually_flipped / len(df)
    bayes_ceiling = 1.0 - effective_flip
    print(f"  Label noise: {actually_flipped}/{len(df)} flipped "
          f"({effective_flip*100:.1f}%)")
    print(f"  Theoretical Bayes ceiling: ~{bayes_ceiling*100:.0f}%")

    # --- Step 5: Shuffle and save ---
    print(f"\n[5/5] Saving dataset...")
    df = df.sample(frac=1, random_state=RANDOM_SEED).reset_index(drop=True)
    df.to_csv(OUTPUT_CSV, index=False)

    # Final summary
    print(f"\n{'='*60}")
    print(f" Dataset generated: {OUTPUT_CSV}")
    print(f"{'='*60}")
    print(f"  Total samples:      {len(df)}")
    print(f"  Features:           {len(feat_cols)} velocity PDF bins")
    print(f"  Classes:            5 (SE-{', SE-'.join(all_labels)})")
    print(f"  Material types:     {sorted(df['material_type'].unique())}")
    print(f"  Feature noise:      {FEATURE_NOISE_RATIO*100:.0f}%")
    print(f"  Label flip rate:    {effective_flip*100:.1f}%")
    print(f"  Bayes ceiling:      ~{bayes_ceiling*100:.0f}%")
    print(f"\n  Per-material sample counts:")
    for mat in sorted(df['material_type'].unique()):
        cnt = len(df[df['material_type'] == mat])
        print(f"    {mat}: {cnt}")
    print(f"\n  Per-SE sample counts:")
    for se in all_labels:
        cnt = len(df[df['surface_energy'] == se])
        print(f"    SE-{se}: {cnt}")
    print(f"\n  Next step: run  run_midterm_explainability.py")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()