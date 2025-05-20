#!/usr/bin/env python3
"""
使用示例:
    python3 Bam_Maptype.py -bam input.bam -o output_prefix -qual 20 -proc 4

描述:
    统计BAM文件中各染色体上的比对记录情况，包括：
      1. 总序列数量
      2. 重复序列数量（标记为PCR duplicate的序列）
      3. 低质量比对序列数量（mapping quality小于指定阈值）
      4. 重复且低质量的序列数量
      5. 保留的序列数量（既非重复且质量达到要求）
    脚本采用多进程加速，每个进程负责一个染色体的数据处理，
    输出结果中第一行为运行命令的注释信息，后续按染色体（按数字顺序排序）分组输出统计结果。
"""

import argparse
import pysam
import sys
import concurrent.futures
import os
import re

def parse_args():
    # 使用 argparse 解析命令行参数
    parser = argparse.ArgumentParser(
        description="统计BAM文件中重复序列和低质量比对序列的数量（按染色体分组，多进程加速）",
        epilog="示例: python3  Bam_Maptype.py.py -bam input.bam -o output_prefix -qual 20 -proc 4"
    )
    parser.add_argument("-bam", required=True, help="输入的BAM文件路径")
    parser.add_argument("-o", required=True, help="输出文件前缀，结果将写入 <prefix>_analysis.txt")
    parser.add_argument("-qual", type=int, default=0, help="低质量阈值（mapping quality），默认为0（不筛选）")
    parser.add_argument("-proc", type=int, default=os.cpu_count(), help="进程数，默认为系统CPU核心数")
    return parser.parse_args()

def process_chrom(chrom, bam_path, qual):
    """
    处理单个染色体的统计信息
    参数:
      chrom: 染色体名称
      bam_path: BAM文件路径
      qual: 低质量阈值
    返回:
      (chrom, stats) 其中 stats 为字典，包含统计信息
    """
    # 初始化统计字典
    stats = {"total": 0, "dup": 0, "lowqual": 0, "dup_lowqual": 0, "retained": 0}
    try:
        # 每个进程独立打开BAM文件
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        sys.stderr.write(f"错误：无法打开BAM文件 {bam_path} 对染色体 {chrom}: {e}\n")
        return chrom, stats

    # 遍历当前染色体的比对记录
    for read in bam.fetch(chrom):
        if read.is_unmapped:  # 跳过未比对的序列
            continue
        stats["total"] += 1
        # 判断低质量：若 mapping quality 小于指定阈值，则视为低质量（如果阈值大于0）
        is_lowqual = (read.mapping_quality < qual) if qual > 0 else False
        # 判断重复：使用 read.is_duplicate 判断PCR duplicate标志
        is_dup = read.is_duplicate

        if is_dup:
            stats["dup"] += 1
        if is_lowqual:
            stats["lowqual"] += 1
        if is_dup and is_lowqual:
            stats["dup_lowqual"] += 1
        if (not is_dup) and (not is_lowqual):
            stats["retained"] += 1
    bam.close()
    return chrom, stats

def natural_key(s):
    """
    用于自然排序的key函数，提取字符串中的数字部分
    如果找不到数字，则返回一个较大值，保证排序时排在后面
    """
    m = re.search(r'(\d+)', s)
    return int(m.group(1)) if m else float('inf')

def main():
    args = parse_args()

    # 尝试打开BAM文件，获取所有参考染色体名称
    try:
        bam = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        sys.stderr.write(f"错误：无法打开BAM文件 {args.bam}: {e}\n")
        sys.exit(1)
    chroms = bam.references  # 获取所有染色体名称（元组）
    bam.close()

    # 使用多进程处理每个染色体，指定进程数为 -proc 参数
    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.proc) as executor:
        future_to_chrom = {executor.submit(process_chrom, chrom, args.bam, args.qual): chrom for chrom in chroms}
        for future in concurrent.futures.as_completed(future_to_chrom):
            chrom, stats = future.result()
            results[chrom] = stats

    # 构造输出文件路径，输出文件为 <prefix>_mapReadsType.txt
    output_file = f"{args.o}_mapReadsType.txt"
    with open(output_file, "w") as fout:
        # 输出注释行，记录运行命令
        fout.write("## 运行命令: " + " ".join(sys.argv) + "\n")
        # 写入表头
        fout.write("Chromosome\tTotal\tDuplicates\tLowQual\tDup_LowQual\tRetained\n")
        # 按照自然排序排序染色体名称，如 scaffold1, scaffold2, scaffold3, ...
        sorted_chroms = sorted(results.keys(), key=natural_key)
        for chrom in sorted_chroms:
            s = results[chrom]
            fout.write(f"{chrom}\t{s['total']}\t{s['dup']}\t{s['lowqual']}\t{s['dup_lowqual']}\t{s['retained']}\n")
    print(f"统计完成，结果写入 {output_file}")

if __name__ == "__main__":
    main()
