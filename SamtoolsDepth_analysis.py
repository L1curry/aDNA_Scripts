#!/usr/bin/env python3
"""
Usage:
    python3 SamtoolsDepth_analysis.py -depth depth.txt -chrlength chrlength.txt -o output_prefix

Description:
    根据samtools计算的depth文件以及参考基因组中每条染色体的长度，
    提取每条染色体上未测序到的区间，并计算覆盖率（Coverage）和未覆盖区域占比（Missing(%)）。
    Coverage = (总长 - 缺失长度)/总长 * 100，Missing(%) = 100 - Coverage。
    输出结果为文本文件，每行包含染色体、缺失区间、Coverage(%) 和 Missing(%) 四列。
"""

import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="提取染色体上未测序区间，并计算覆盖率和未覆盖区域占比",
        epilog="示例: python3 SamtoolsDepth_analysis.py -depth depth.txt -chrlength chrlength.txt -o output_prefix"
    )
    parser.add_argument("-depth", required=True, help="输入的depth文件（samtools depth 结果），格式：chr<TAB>pos<TAB>depth")
    parser.add_argument("-chrlength", required=True, help="染色体长度文件，每行两列：chr<TAB>length")
    parser.add_argument("-o", required=True, help="输出文件前缀，结果将写入 <prefix>_depthAnalysis.txt")
    return parser.parse_args()

def read_chrlength(chrlength_file):
    """
    读取染色体长度文件，返回字典 {chrom: length}
    """
    chr_lengths = {}
    with open(chrlength_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            chrom, length = parts[0], parts[1]
            try:
                chr_lengths[chrom] = int(length)
            except ValueError:
                sys.stderr.write(f"警告：无法转换长度为整数：{line}\n")
    return chr_lengths

def read_depth(depth_file):
    """
    读取depth文件，返回字典，键为染色体名，值为已测序位置的列表（整数），列表为升序
    假设depth文件格式：chr<TAB>pos<TAB>depth
    """
    depth_dict = {}
    with open(depth_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if chrom not in depth_dict:
                depth_dict[chrom] = []
            depth_dict[chrom].append(pos)
    for chrom in depth_dict:
        depth_dict[chrom].sort()
    return depth_dict

def get_missing_intervals(chr_length, positions):
    """
    给定某条染色体长度和已测序位置列表（升序），返回未测序的连续区间列表，
    每个区间为 (start, end)（1-based闭区间）。
    若positions为空，则返回 [(1, chr_length)]
    """
    missing = []
    if not positions:
        missing.append((1, chr_length))
        return missing
    if positions[0] > 1:
        missing.append((1, positions[0]-1))
    prev = positions[0]
    for pos in positions[1:]:
        if pos - prev > 1:
            missing.append((prev+1, pos-1))
        prev = pos
    if positions[-1] < chr_length:
        missing.append((positions[-1]+1, chr_length))
    return missing

def main():
    args = parse_args()
    # 读取染色体长度文件
    chr_lengths = read_chrlength(args.chrlength)
    # 读取depth文件，获取每条染色体已测序的位置列表
    depth_dict = read_depth(args.depth)
    
    output_file = f"{args.o}_depthAnalysis.txt"
    
    with open(output_file, 'w') as out:
        # 输出增加一列 Missing(%)
        out.write("Chromosome\tMissing_Intervals\tCoverage(%)\tMissing(%)\n")
        for chrom, length in chr_lengths.items():
            positions = depth_dict.get(chrom, [])
            missing_intervals = get_missing_intervals(length, positions)
            missing_length = sum(end - start + 1 for start, end in missing_intervals)
            coverage = (length - missing_length) / length * 100
            missing_percent = 100 - coverage
            if missing_intervals:
                intervals_str = ",".join(f"{s}-{e}" for s, e in missing_intervals)
            else:
                intervals_str = "None"
            out.write(f"{chrom}\t{intervals_str}\t{coverage:.2f}\t{missing_percent:.2f}\n")
    
    print(f"处理完成，结果写入 {output_file}")

if __name__ == '__main__':
    main()
