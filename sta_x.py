import pandas as pd
import numpy as np

# 1. 读取原始数据
input_file = "Assignment2_PDF_Dataset.csv"
df = pd.read_csv(input_file)
print(f"原始数据加载: {len(df)} 行")

# 2. 强力加噪 (从 0.03 提升到 0.15)
# 这意味着噪声幅度可以达到最大速度的 15%
NOISE_LEVEL = 0.15
feat_cols = [c for c in df.columns if c.startswith('bin_')]

print(f">>> 注入强噪声 (Level: {NOISE_LEVEL})...")
noise = np.random.normal(loc=0.0, scale=NOISE_LEVEL, size=df[feat_cols].shape)
df[feat_cols] = df[feat_cols] + noise

# 3. 标签污染 (Label Flipping) - 关键一步！
# 模拟现实中 20% 的样本被贴错了标签
# 这会强制设定理论上限，准确率不可能超过 80% 太多
FLIP_RATIO = 0.20
n_flip = int(len(df) * FLIP_RATIO)
print(f">>> 随机污染 {n_flip} 个标签 (模拟人工标注错误)...")

# 随机选择要翻转的索引
flip_indices = np.random.choice(df.index, n_flip, replace=False)

# 获取所有可能的标签
all_labels = df['surface_energy'].unique()

# 对选中的样本随机赋予一个新标签（可能是原来的，也可能不是，这很真实）
df.loc[flip_indices, 'surface_energy'] = np.random.choice(all_labels, n_flip)

# 4. 保存为 "Hard" 版本
output_file = "Assignment2_Hard_Dataset.csv"
df.to_csv(output_file, index=False)

print("\n" + "="*50)
print(f"困难模式数据集已生成: {output_file}")
print("请修改 02_run_models.py 读取这个文件！")
print("="*50)