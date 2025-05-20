#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage: Bam_analysis1.py [options] --i FILE

该脚本基于samtools depth结果绘制图形，用于评估BAM文件的覆盖情况。
可选功能：
  --zxt   绘制折线图（显示指定染色体或所有染色体上深度沿位置的变化）
  --zft   绘制直方图（显示深度值的频数分布）
  --ljfbt 绘制累计分布图（显示深度值的累积分布）
  --chr   指定绘制的染色体：传入 "all" 绘制所有染色体；传入数字选取文件中第 n 条染色体（例如：--chr 1）
  --window 指定染色体上的位置范围，例如 --window 10 1000
  --o     输出文件前缀，生成的图片将以该前缀保存为pdf格式

示例：
  python Bam_analysis1.py --i depth.txt --zft --chr all --window 10 1000 --o result
"""

import argparse
import sys
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="该脚本基于samtools depth结果绘制图形，用于评估BAM文件的覆盖情况。",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--i", "--input", required=True, type=str,
                        help="输入的depth结果文件，要求为3列：染色体、位置、深度")
    parser.add_argument("--zxt", action="store_true", default=False,
                        help="开启折线图功能（显示深度沿位置的变化）")
    parser.add_argument("--zft", action="store_true", default=False,
                        help="开启直方图功能（显示深度值的频数分布）")
    parser.add_argument("--ljfbt", action="store_true", default=False,
                        help="开启累计分布图功能（显示深度值的累积分布）")
    parser.add_argument("--chr", type=str, default="all",
                        help="指定绘制的染色体：传入 'all' 绘制所有染色体；传入数字选取文件中第 n 条染色体（例如：--chr 1）")
    parser.add_argument("--window", type=int, nargs=2, metavar=("START", "END"),
                        help="指定染色体上的位置范围，例如 --window 10 1000")
    parser.add_argument("--o", "--output", type=str, default="output",
                        help="输出文件前缀，生成的图片将以该前缀保存为pdf格式")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.i):
        sys.exit("错误：输入文件不存在！")
    
    # 读取 depth 结果文件，假设以空白字符或制表符分隔
    try:
        depth_df = pd.read_csv(args.i, sep=r"\s+", header=None)
    except Exception as e:
        sys.exit(f"读取文件时出错: {e}")
    
    if depth_df.shape[1] < 3:
        sys.exit("错误：输入文件格式错误，至少需要3列（染色体、位置、深度）")
    
    # 设置列名
    depth_df.columns = ["chrom", "pos", "depth"]
    
    # 获取文件中所有染色体（按照出现顺序，不重复）
    unique_chr = depth_df["chrom"].drop_duplicates().tolist()
    num_chr = len(unique_chr)
    
    # 处理 --chr 参数
    chr_arg = args.chr.strip().lower()
    if chr_arg == "all":
        selected_chr = unique_chr
    else:
        try:
            chr_index = int(chr_arg)
        except ValueError:
            sys.exit("错误：--chr 参数输入不合法，请输入 'all' 或正整数。")
        if chr_index < 1 or chr_index > num_chr:
            sys.exit(f"错误：超出染色体条数范围：文件中只有 {num_chr} 条染色体，您选择了第 {chr_index} 条。")
        selected_chr = [ unique_chr[chr_index - 1] ]
    
    # 筛选数据：只保留所选染色体的数据
    plot_data = depth_df[ depth_df["chrom"].isin(selected_chr) ]
    
    # 处理 --window 参数
    window = args.window
    if window:
        start, end = window
        if start > end:
            sys.exit("错误：窗口的起始位置不能大于结束位置。")
        # 筛选数据
        plot_data = plot_data[(plot_data["pos"] >= start) & (plot_data["pos"] <= end)]
        # 检查是否有数据
        if plot_data.empty:
            sys.exit("错误：指定窗口范围内没有数据。")
    else:
        start, end = plot_data["pos"].min(), plot_data["pos"].max()
    
    # 根据参数绘图
    output_prefix = args.o

    # 1. 绘制折线图（--zxt）：按位置展示深度变化
    if args.zxt:
        # 如果有多个染色体，则为每个染色体绘制子图
        if len(selected_chr) > 1:
            fig, axes = plt.subplots(len(selected_chr), 1, figsize=(10, 5*len(selected_chr)), sharex=False)
            if len(selected_chr) == 1:
                axes = [axes]
            for ax, chrom in zip(axes, selected_chr):
                subdata = plot_data[ plot_data["chrom"] == chrom ]
                ax.plot(subdata["pos"], subdata["depth"], color="steelblue")
                ax.set_title(f"Chromosome {chrom} Depth Line Plot")
                ax.set_xlabel("Position")
                ax.set_ylabel("Depth")
                if window:
                    ax.set_xlim(start, end)
            plt.tight_layout()
        else:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(plot_data["pos"], plot_data["depth"], color="steelblue")
            ax.set_title(f"Chromosome {selected_chr[0]} Depth Line Plot")
            ax.set_xlabel("Position")
            ax.set_ylabel("Depth")
            if window:
                ax.set_xlim(start, end)
            plt.tight_layout()
        pdf_file = f"{output_prefix}_line.pdf"
        fig.savefig(pdf_file)
        plt.close(fig)
        print(f"折线图已保存到 {pdf_file}")
    
    # 2. 绘制直方图（--zft）：显示深度值的频数分布
    if args.zft:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(plot_data["depth"], bins=50, color="skyblue", edgecolor="black")
        ax.set_title(f"Depth Histogram for Chromosome {', '.join(selected_chr)}")
        ax.set_xlabel("Depth")
        ax.set_ylabel("Frequency")
        plt.tight_layout()
        pdf_file = f"{output_prefix}_hist.pdf"
        fig.savefig(pdf_file)
        plt.close(fig)
        print(f"直方图已保存到 {pdf_file}")
    
    # 3. 绘制累计分布图（--ljfbt）：显示深度累积分布
    if args.ljfbt:
        # 计算每个深度值的频数和累积比例
        cum_df = plot_data.groupby("depth").size().reset_index(name="freq")
        cum_df = cum_df.sort_values("depth")
        cum_df["cumfreq"] = cum_df["freq"].cumsum()
        total = cum_df["freq"].sum()
        cum_df["cumperc"] = cum_df["cumfreq"] / total
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(cum_df["depth"], cum_df["cumperc"], color="darkgreen")
        ax.set_title(f"Depth Cumulative Distribution for Chromosome {', '.join(selected_chr)}")
        ax.set_xlabel("Depth")
        ax.set_ylabel("Cumulative Proportion")
        plt.tight_layout()
        pdf_file = f"{output_prefix}_cum.pdf"
        fig.savefig(pdf_file)
        plt.close(fig)
        print(f"累计分布图已保存到 {pdf_file}")
    
    if not (args.zxt or args.zft or args.ljfbt):
        print("警告：未选择任何绘图参数 (--zxt, --zft, --ljfbt)，脚本未生成任何图形。")

if __name__ == "__main__":
    main()
