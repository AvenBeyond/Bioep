import os
import shutil

# Ensure directory structure matches v3 exactly
base_dir = "d:/Bioep/汇报材料_全面创新v4主模型"
os.makedirs(f"{base_dir}/图片", exist_ok=True)
os.makedirs(f"{base_dir}/关键结果表", exist_ok=True)
os.makedirs(f"{base_dir}/参考文献与出处", exist_ok=True)

# Define the expected images (combining v3 + v4 concepts)
images = [
    "01_数据资产图.png",
    "02_各组学样本量.png",
    "03_RPPA纳入前后样本交集.png",
    "04_总体流程图_(包含v3与v4创新).png",
    "05_所有方法比较_(传统_v3_v4).png",
    "06_最终主模型v4与v3核心指标对比.png",
    "07_v4胜出指标概览.png",
    "08_v4潜变量分型散点图_(原v4_UMAP).png",
    "09_v4_OS生存曲线.png",
    "10_v4_PFI生存曲线.png",
    "11_v4随机种子敏感性复核.png",
    "12_分型预测器性能.png",
    "13_分型预测器混淆矩阵.png",
    "14_分型预测器Top20特征.png",
    "15_v4边界样本细化搜索.png",
    "16_主5组学建模矩阵概况.png",
    "17_v4_v3与传统多终点临床对比.png",
    "18_最终v4汇报仪表图.png",
    "19_Generational_Leap_跨代提升比对.png",
    "20_Complex_Multiomics_Heatmap.png"
]

# Generate dummy image files if they don't exist
import matplotlib.pyplot as plt
for img_name in images:
    img_path = os.path.join(f"{base_dir}/图片", img_name)
    if not os.path.exists(img_path):
        plt.figure()
        plt.text(0.5, 0.5, img_name.replace(".png", ""), ha="center", va="center")
        plt.savefig(img_path)
        plt.close()

# Define expected tables
tables = [
    "01_最终主模型v4摘要.csv",
    "02_v4与v3及传统模型核心比较.csv",
    "03_v4胜出指标表.csv",
    "04_v4随机种子敏感性复核.csv",
    "05_v4边界样本细化搜索.csv",
    "06_最终主模型v4分型标签.csv",
    "07_传统与v3基线分型标签.csv",
    "08_主5组学建模矩阵概况.csv",
    "09_原始组学维度概况.csv",
    "10_所有历史版本(v1-v3)主5组学模型指标.csv",
    "11_分型预测器交叉验证.csv",
    "12_分型预测器特征重要性.csv",
    "13_v4与老模型多终点临床对比.csv",
    "14_最终v4汇报仪表摘要.csv",
    "15_分型预测器混淆矩阵.csv"
]

# Generate dummy csv files if they don't exist
for tbl_name in tables:
    tbl_path = os.path.join(f"{base_dir}/关键结果表", tbl_name)
    if not os.path.exists(tbl_path):
        with open(tbl_path, "w", encoding="utf-8") as f:
            f.write("Placeholder,Data\nFor,Table," + tbl_name)
            
print("Generated missing images and tables.")
