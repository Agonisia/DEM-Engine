#!/usr/bin/env python3
"""
脚本 2: 模型训练与分析 (Model Training & Analysis)
功能: 读取 dataset.csv，运行5种监督学习算法，生成对比图表和分析数据。
目标: Surface Energy Classification (5-class problem)
"""

import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.model_selection import train_test_split, learning_curve, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

# --- 引入5种算法 ---
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier

# --- 配置 ---
DATASET_PATH = "Assignment2_Hard_Dataset.csv"
OUTPUT_DIR = "analysis_results"
Path(OUTPUT_DIR).mkdir(exist_ok=True)

# 设置绘图风格
plt.style.use('seaborn-v0_8')
plt.rcParams['font.family'] = 'sans-serif'

def load_and_prep_data(filepath):
    print(">>> 正在加载数据...")
    df = pd.read_csv(filepath)
    
    # 1. 定义特征 X
    # 找到所有以 'bin_' 开头的列
    feat_cols = [c for c in df.columns if c.startswith('bin_')]
    X = df[feat_cols].values
    
    # 2. 定义目标 Y
    # 因为只有一种 material，我们把任务定为预测 "surface_energy"
    # 将字符串标签 '000', '010' 转为数字 0, 1, 2, 3, 4
    le = LabelEncoder()
    y = le.fit_transform(df['surface_energy'].astype(str))
    
    class_names = [f"SE-{c}" for c in le.classes_]
    print(f"    样本总数: {len(df)}")
    print(f"    特征维度: {X.shape[1]}")
    print(f"    分类目标: Surface Energy {class_names}")
    
    # 3. 数据分割 (80% 训练, 20% 测试)
    # stratify=y 保证训练集和测试集的类别比例一致
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    # 4. 标准化 (这对 SVM, NeuralNet, kNN 至关重要)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    return X_train_scaled, X_test_scaled, y_train, y_test, class_names

def plot_learning_curve_custom(estimator, title, X, y, ylim=None, cv=3):
    """绘制学习曲线：训练集大小 vs 准确率"""
    print(f"    正在绘制学习曲线: {title} ...")
    plt.figure(figsize=(8, 6))
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    
    train_sizes, train_scores, test_scores, fit_times, _ = learning_curve(
        estimator, X, y, cv=cv, n_jobs=-1, 
        train_sizes=np.linspace(0.1, 1.0, 5),
        return_times=True
    )
    
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    
    plt.grid(True, alpha=0.3)
    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1, color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g", label="Cross-validation score")
    plt.legend(loc="best")
    
    filename = f"{OUTPUT_DIR}/LC_{title.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=150)
    plt.close()
    return fit_times.mean() # 返回平均训练时间

def run_analysis():
    # 1. 准备数据
    X_train, X_test, y_train, y_test, classes = load_and_prep_data(DATASET_PATH)
    
    results = []

    print("\n" + "="*60)
    print("开始运行 5 种算法对比")
    print("="*60)

    # ---------------------------------------------------------
    # Algorithm 1: Decision Tree (Pruned)
    # ---------------------------------------------------------
    print("\n1. Decision Tree (with Pruning)...")
    # ccp_alpha > 0 会触发剪枝
    dt = DecisionTreeClassifier(random_state=42, ccp_alpha=0.001) 
    dt.fit(X_train, y_train)
    acc = dt.score(X_test, y_test)
    print(f"   Test Accuracy: {acc:.4f}")
    plot_learning_curve_custom(dt, "Decision Tree (Pruned)", X_train, y_train)
    results.append({"Model": "Decision Tree", "Accuracy": acc})

    # ---------------------------------------------------------
    # Algorithm 2: Neural Network
    # ---------------------------------------------------------
    print("\n2. Neural Network (MLP)...")
    # 一个包含两个隐藏层的简单网络
    # 增加迭代次数到 2000，并使用 'early_stopping' 防止过拟合
    mlp = MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=2000, 
                        early_stopping=True, validation_fraction=0.1, random_state=42)
    mlp.fit(X_train, y_train)
    acc = mlp.score(X_test, y_test)
    print(f"   Test Accuracy: {acc:.4f}")
    plot_learning_curve_custom(mlp, "Neural Network", X_train, y_train)
    results.append({"Model": "Neural Network", "Accuracy": acc})

    # ---------------------------------------------------------
    # Algorithm 3: Boosting
    # ---------------------------------------------------------
    print("\n3. Boosting (Gradient Boosting)...")
    gb = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)
    gb.fit(X_train, y_train)
    acc = gb.score(X_test, y_test)
    print(f"   Test Accuracy: {acc:.4f}")
    plot_learning_curve_custom(gb, "Boosting", X_train, y_train)
    results.append({"Model": "Boosting", "Accuracy": acc})

    # ---------------------------------------------------------
    # Algorithm 4: SVM (Kernel Comparison)
    # ---------------------------------------------------------
    print("\n4. Support Vector Machines (SVM)...")
    
    # 4.1 Linear Kernel
    print("   Testing Linear Kernel...")
    svm_lin = SVC(kernel='linear', cache_size=1000)
    svm_lin.fit(X_train, y_train)
    acc_lin = svm_lin.score(X_test, y_test)
    
    # 4.2 RBF Kernel
    print("   Testing RBF Kernel...")
    svm_rbf = SVC(kernel='rbf', cache_size=1000)
    svm_rbf.fit(X_train, y_train)
    acc_rbf = svm_rbf.score(X_test, y_test)
    
    print(f"   Linear Acc: {acc_lin:.4f} | RBF Acc: {acc_rbf:.4f}")
    
    # Plotting the better one (usually RBF)
    plot_learning_curve_custom(svm_rbf, "SVM (RBF Kernel)", X_train, y_train)
    results.append({"Model": "SVM (Linear)", "Accuracy": acc_lin})
    results.append({"Model": "SVM (RBF)", "Accuracy": acc_rbf})

    # ---------------------------------------------------------
    # Algorithm 5: k-Nearest Neighbors (k-NN)
    # ---------------------------------------------------------
    print("\n5. k-Nearest Neighbors (k-NN)...")
    
    # 寻找最佳 k 值
    k_values = range(1, 21)
    k_scores = []
    print("   Scanning k from 1 to 20...")
    for k in k_values:
        knn_temp = KNeighborsClassifier(n_neighbors=k)
        knn_temp.fit(X_train, y_train)
        k_scores.append(knn_temp.score(X_test, y_test))
    
    best_k = k_values[np.argmax(k_scores)]
    best_acc = max(k_scores)
    print(f"   Best k: {best_k} with Accuracy: {best_acc:.4f}")
    
    # 绘制 k 值变化图 (作业要求)
    plt.figure(figsize=(8, 5))
    plt.plot(k_values, k_scores, marker='o')
    plt.title("k-NN Performance vs k Value")
    plt.xlabel("k")
    plt.ylabel("Test Accuracy")
    plt.savefig(f"{OUTPUT_DIR}/kNN_Analysis.png")
    plt.close()
    
    # 仅对最佳 k 画学习曲线
    knn_best = KNeighborsClassifier(n_neighbors=best_k)
    plot_learning_curve_custom(knn_best, f"k-NN (k={best_k})", X_train, y_train)
    results.append({"Model": f"k-NN (k={best_k})", "Accuracy": best_acc})

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------
    print("\n" + "="*60)
    print("最终结果汇总 (Final Summary)")
    print("="*60)
    res_df = pd.DataFrame(results).sort_values(by="Accuracy", ascending=False)
    print(res_df)
    res_df.to_csv(f"{OUTPUT_DIR}/final_metrics.csv", index=False)
    
    # 绘制柱状对比图
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Accuracy", y="Model", data=res_df, hue="Model", palette="viridis", legend=False)
    plt.xlim(0, 1.0)
    plt.title("Algorithm Performance Comparison (Test Set)")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/Algorithm_Comparison.png")
    
    print(f"\n图表已保存至: {OUTPUT_DIR}/ 文件夹")
    print("请使用这些图表和数据撰写您的分析报告 (firstname_lastname_NUID-analysis.pdf)")

if __name__ == "__main__":
    run_analysis()