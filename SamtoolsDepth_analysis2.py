#!/usr/bin/env python3
"""
Usage:
    python3 SamtoolsDepth_analysis2.py -i input.depth -o Test1

Description:
    根据 samtools depth 的结果文件，对测序深度进行统计分析与可视化。
    计算总体平均深度、中位数、众数、标准差、最小值和最大值，
    同时对每条染色体计算平均深度、中位数、众数、标准差、最小深度和最大深度。
    并绘制总体深度分布直方图。
    
    输出统计结果写入 <output_prefix>_depthAnalysis2.txt，直方图保存为 <output_prefix>_histogram.png。
    输出结果中会注明样本名（输出前缀），各染色体统计信息中也包含标准差信息。
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import Counter
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description="对 samtools depth 文件进行统计分析和绘制深度分布直方图",
        epilog="示例: python3 SamtoolsDepth_analysis2.py -i input.depth -o Test1"
    )
    parser.add_argument("-i", "--input", required=True, help="samtools depth 的结果文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出前缀，例如 Test1；结果写入 <prefix>_depthAnalysis2.txt，直方图保存为 <prefix>_histogram.png")
    return parser.parse_args()

def calc_mode(series):
    """
    计算 series 的众数，返回出现次数最多的值；
    如果有多个众数，返回第一个。
    """
    counts = Counter(series.dropna())
    if counts:
        return counts.most_common(1)[0][0]
    else:
        return None

def natural_key(s):
    """
    用于自然排序的 key 函数，
    提取字符串中的数字部分；若找不到数字则返回无穷大，
    保证类似 "scaffold1", "scaffold2", … 按数字顺序排序。
    """
    m = re.search(r'(\d+)', s)
    return int(m.group(1)) if m else float('inf')

def main():
    args = parse_args()
    
    # 读取 depth 文件，假设为制表符分隔，列依次为：chromosome, position, depth
    try:
        df = pd.read_csv(args.input, sep="\t", header=None, names=["chrom", "pos", "depth"])
    except Exception as e:
        print(f"读取文件 {args.input} 失败：{e}")
        return
    
    # 转换深度为数值型
    df["depth"] = pd.to_numeric(df["depth"], errors='coerce')
    
    # 计算总体统计信息
    overall_avg = df["depth"].mean()
    overall_median = df["depth"].median()
    overall_mode = calc_mode(df["depth"])
    overall_std = df["depth"].std()
    overall_min = df["depth"].min()
    overall_max = df["depth"].max()
    
    # 按染色体分组统计，每组计算平均、中位、众数、标准差、最小和最大深度
    chrom_stats = {}
    for chrom, group in df.groupby("chrom"):
        avg = group["depth"].mean()
        median = group["depth"].median()
        mode = calc_mode(group["depth"])
        std = group["depth"].std()
        min_depth = group["depth"].min()
        max_depth = group["depth"].max()
        chrom_stats[chrom] = {
            "avg": avg,
            "median": median,
            "mode": mode,
            "std": std,
            "min": min_depth,
            "max": max_depth
        }
    
    # 输出统计结果到文本文件，文件名为 <output_prefix>_depthAnalysis2.txt
    out_txt = f"{args.output}_depthAnalysis2.txt"
    with open(out_txt, "w") as f:
        # 在总体统计信息前注明样本名（输出前缀）
        f.write(f"{args.output} 的总体统计信息:\n")
        f.write(f"平均深度: {overall_avg:.2f}\n")
        f.write(f"中位数深度: {overall_median:.2f}\n")
        f.write(f"众数深度: {overall_mode}\n")
        f.write(f"标准差: {overall_std:.2f}\n")
        f.write(f"最小深度: {overall_min}\n")
        f.write(f"最大深度: {overall_max}\n\n")
        
        f.write(f"{args.output}的各染色体统计信息:\n")
        f.write("Chromosome\tAverage\tMedian\tMode\tStd\tMin\tMax\n")
        # 按照自然排序排序染色体名称
        for chrom in sorted(chrom_stats.keys(), key=natural_key):
            stats = chrom_stats[chrom]
            f.write(f"{chrom}\t{stats['avg']:.2f}\t{stats['median']:.2f}\t{stats['mode']}\t{stats['std']:.2f}\t{stats['min']}\t{stats['max']}\n")
    
    print(f"统计结果已写入 {out_txt}")
    
    # 绘制总体深度分布直方图
    plt.figure(figsize=(10, 6))
    plt.hist(df["depth"].dropna(), bins=50, color='steelblue', edgecolor='black')
    plt.xlabel("Depth")
    plt.ylabel("Frequency")
    plt.title("Depth Distribution")
    hist_file = f"{args.output}_histogram.png"
    plt.savefig(hist_file, dpi=300)
    plt.close()
    
    print(f"直方图已保存为 {hist_file}")

if __name__ == '__main__':
    main()
