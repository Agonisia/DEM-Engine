#!/usr/bin/env python3
"""
脚本 3: ANN 深度剖析 (ANN In-Depth Analysis) - 改进版
功能: 针对含噪 DEM 数据集，全面对比不同 MLP 架构和超参数的表现。
改进: 
  1. 增加 Training Accuracy 输出（用于对比 train/test gap 体现 overfitting）
  2. 增加所有实验的 Learning Curve 对比图（合并到一张图）
  3. 保留原有的 Loss Curves 和 Confusion Matrix
约束: 仅使用 matplotlib 和 sklearn，无额外依赖。
"""

import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.model_selection import train_test_split, learning_curve
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, ConfusionMatrixDisplay

# --- 配置 ---
DATASET_PATH = "Assignment2_Hard_Dataset.csv"
OUTPUT_DIR = "ann_analysis_results"
Path(OUTPUT_DIR).mkdir(exist_ok=True)

plt.style.use('default')


def load_and_prep_data(filepath):
    print(">>> 正在加载数据...")
    df = pd.read_csv(filepath)
    feat_cols = [c for c in df.columns if c.startswith('bin_')]
    X = df[feat_cols].values

    le = LabelEncoder()
    y = le.fit_transform(df['surface_energy'].astype(str))
    class_names = [f"SE-{c}" for c in le.classes_]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    return X_train_scaled, X_test_scaled, y_train, y_test, class_names


def plot_loss_curves(models_dict):
    """绘制多个模型的 Loss 曲线对比"""
    print("    正在绘制 Loss Curves...")
    plt.figure(figsize=(10, 6))

    for name, model in models_dict.items():
        if hasattr(model, 'loss_curve_'):
            plt.plot(model.loss_curve_, label=name, linewidth=2)

    plt.title("Training Loss Curve Comparison", fontsize=14)
    plt.xlabel("Iterations (Epochs)", fontsize=12)
    plt.ylabel("Loss", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/ANN_Loss_Curves.png", dpi=150)
    plt.close()


def plot_all_learning_curves(experiments_dict, X_train, y_train, cv=3):
    """
    绘制所有实验的 Learning Curve 合并到一张图（cross-validation score vs training size）
    这是作业要求的关键图：performance as a function of training size
    """
    print("    正在绘制合并 Learning Curves（所有实验）...")
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()
    
    # 也准备一个 overlay 图
    fig_overlay, ax_overlay = plt.subplots(figsize=(10, 6))
    
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']
    
    for idx, (name, mlp) in enumerate(experiments_dict.items()):
        print(f"      计算 {name} 的 learning curve...")
        train_sizes, train_scores, test_scores = learning_curve(
            mlp, X_train, y_train, cv=cv, n_jobs=-1,
            train_sizes=np.linspace(0.1, 1.0, 5),
            scoring='accuracy'
        )

        train_mean = np.mean(train_scores, axis=1)
        train_std = np.std(train_scores, axis=1)
        test_mean = np.mean(test_scores, axis=1)
        test_std = np.std(test_scores, axis=1)

        # --- 子图 ---
        ax = axes[idx]
        ax.fill_between(train_sizes, train_mean - train_std,
                        train_mean + train_std, alpha=0.1, color="r")
        ax.fill_between(train_sizes, test_mean - test_std,
                        test_mean + test_std, alpha=0.1, color="g")
        ax.plot(train_sizes, train_mean, 'o-', color="r", label="Training score")
        ax.plot(train_sizes, test_mean, 'o-', color="g", label="CV score")
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("Training examples", fontsize=9)
        ax.set_ylabel("Accuracy", fontsize=9)
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.set_ylim([0.5, 1.05])

        # --- Overlay 图 (仅 CV score) ---
        ax_overlay.plot(train_sizes, test_mean, 'o-', color=colors[idx],
                        label=name, linewidth=2)
        ax_overlay.fill_between(train_sizes, test_mean - test_std,
                                test_mean + test_std, alpha=0.1, color=colors[idx])

    # 隐藏多余的子图
    for j in range(len(experiments_dict), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Learning Curves: All ANN Experiments", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{OUTPUT_DIR}/ANN_Learning_Curves_Individual.png", dpi=150)
    plt.close(fig)

    # Overlay 图
    ax_overlay.set_title("Learning Curves Comparison (Cross-Validation Score)", fontsize=14)
    ax_overlay.set_xlabel("Training examples", fontsize=12)
    ax_overlay.set_ylabel("Cross-Validation Accuracy", fontsize=12)
    ax_overlay.legend(loc="lower right", fontsize=9)
    ax_overlay.grid(True, linestyle='--', alpha=0.7)
    ax_overlay.set_ylim([0.5, 0.95])
    fig_overlay.tight_layout()
    fig_overlay.savefig(f"{OUTPUT_DIR}/ANN_Learning_Curves_Comparison.png", dpi=150)
    plt.close(fig_overlay)


def plot_single_learning_curve(estimator, title, X, y, cv=3):
    """绘制单个模型的学习曲线（保留兼容性）"""
    print(f"    正在绘制学习曲线: {title} ...")
    plt.figure(figsize=(8, 6))
    plt.title(f"Learning Curve: {title}")
    plt.xlabel("Training examples")
    plt.ylabel("Accuracy Score")

    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=-1,
        train_sizes=np.linspace(0.1, 1.0, 5)
    )

    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1, color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g", label="Cross-validation score")
    plt.legend(loc="best")

    safe_title = title.replace(' ', '_').replace(',', '').replace('=', '').replace(':', '')
    plt.savefig(f"{OUTPUT_DIR}/LC_{safe_title}.png", dpi=150)
    plt.close()


def run_ann_experiments():
    X_train, X_test, y_train, y_test, class_names = load_and_prep_data(DATASET_PATH)

    print("\n" + "=" * 60)
    print("开始执行 ANN 深度对比实验 (改进版)")
    print("=" * 60)

    # 定义实验组
    experiments = {
        "Exp 1: Shallow & Wide (ReLU)": MLPClassifier(
            hidden_layer_sizes=(128,), activation='relu', max_iter=1000, random_state=42),
        "Exp 2: Deep & Narrow (ReLU)": MLPClassifier(
            hidden_layer_sizes=(64, 32, 16), activation='relu', max_iter=1000, random_state=42),
        "Exp 3: Tanh Activation": MLPClassifier(
            hidden_layer_sizes=(64, 32), activation='tanh', max_iter=1000, random_state=42),
        "Exp 4: High Regularization (Alpha=0.1)": MLPClassifier(
            hidden_layer_sizes=(64, 32), alpha=0.1, max_iter=1000, random_state=42),
        "Exp 5: Early Stopping": MLPClassifier(
            hidden_layer_sizes=(64, 32), early_stopping=True, max_iter=1000, random_state=42)
    }

    results = []
    trained_models = {}
    best_model_name = None
    best_accuracy = 0

    for name, mlp in experiments.items():
        print(f"\n训练模型: {name}")
        start_time = time.time()

        mlp.fit(X_train, y_train)

        train_time = time.time() - start_time
        
        # ===== 改进: 同时记录 train accuracy 和 test accuracy =====
        train_acc = mlp.score(X_train, y_train)
        test_acc = mlp.score(X_test, y_test)
        iterations = mlp.n_iter_

        print(f"   Train Accuracy: {train_acc:.4f} | Test Accuracy: {test_acc:.4f} | "
              f"Iterations: {iterations} | Time: {train_time:.2f}s")

        trained_models[name] = mlp
        results.append({
            "Experiment": name,
            "Train_Accuracy": train_acc,
            "Test_Accuracy": test_acc,
            "Iterations": iterations,
            "Train_Time_sec": round(train_time, 2)
        })

        if test_acc > best_accuracy:
            best_accuracy = test_acc
            best_model_name = name

    # --- 绘制所有模型的 Loss 曲线对比 ---
    plot_loss_curves(trained_models)

    # --- 绘制最佳模型的单独 Learning Curve ---
    print(f"\n绘制最佳模型 ({best_model_name}) 的单独学习曲线...")
    plot_single_learning_curve(trained_models[best_model_name], 
                               best_model_name, X_train, y_train)

    # --- 改进: 绘制所有模型的合并 Learning Curves ---
    print("\n绘制所有模型的合并 Learning Curves...")
    plot_all_learning_curves(experiments, X_train, y_train)

    # --- 绘制最佳模型的混淆矩阵 ---
    print(f"\n绘制最佳模型 ({best_model_name}) 的混淆矩阵...")
    best_mlp = trained_models[best_model_name]
    fig, ax = plt.subplots(figsize=(8, 6))
    disp = ConfusionMatrixDisplay.from_estimator(
        best_mlp, X_test, y_test,
        display_labels=class_names,
        cmap=plt.cm.Blues, ax=ax
    )
    plt.title(f"Confusion Matrix: {best_model_name}")
    plt.savefig(f"{OUTPUT_DIR}/ANN_Best_Confusion_Matrix.png", dpi=150)
    plt.close()

    # --- 汇总输出 ---
    print("\n" + "=" * 60)
    print("实验结果汇总 (Summary)")
    print("=" * 60)
    res_df = pd.DataFrame(results).sort_values(by="Test_Accuracy", ascending=False)
    print(res_df.to_string(index=False))
    res_df.to_csv(f"{OUTPUT_DIR}/ann_experiments_metrics.csv", index=False)

    # --- 打印 LaTeX 表格 (方便直接复制到报告) ---
    print("\n" + "=" * 60)
    print("LaTeX Table (可直接复制到报告):")
    print("=" * 60)
    for _, row in res_df.iterrows():
        exp_short = row['Experiment'].replace('(ReLU)', '').replace('(Alpha=0.1)', r'($\\alpha=0.1$)').strip()
        print(f"  {exp_short} & {row['Train_Accuracy']*100:.2f}\\% & "
              f"{row['Test_Accuracy']*100:.2f}\\% & {row['Iterations']} & "
              f"{row['Train_Time_sec']} \\\\")

    print(f"\n图表和数据已保存至: {OUTPUT_DIR}/ 文件夹")
    print("\n生成的文件清单:")
    print(f"  1. {OUTPUT_DIR}/ANN_Loss_Curves.png          - Loss曲线对比")
    print(f"  2. {OUTPUT_DIR}/ANN_Best_Confusion_Matrix.png - 最佳模型混淆矩阵")
    print(f"  3. {OUTPUT_DIR}/ANN_Learning_Curves_Comparison.png - 所有模型CV曲线对比")
    print(f"  4. {OUTPUT_DIR}/ANN_Learning_Curves_Individual.png - 各模型单独学习曲线")
    print(f"  5. {OUTPUT_DIR}/LC_Exp_5_Early_Stopping.png   - 最佳模型单独学习曲线")
    print(f"  6. {OUTPUT_DIR}/ann_experiments_metrics.csv    - 结果数据表")


if __name__ == "__main__":
    run_ann_experiments()