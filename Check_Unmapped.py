#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本功能：
  1. 从输入 BAM 中筛选 unmapped reads（即第3列为“*”）。
  2. 对这些 unmapped reads 根据其序列（query_sequence）去重：
     - 首次出现的序列写入 <prefix>_pass.bam
     - 重复出现的序列写入 <prefix>_fail.bam（保留所有重复条目）
  3. 跳过所有已比对的 reads（mapped），不写入任何输出文件。

用法：
  python3 Check_Unmapped.py \
    -bam input.bam \
    -o output_prefix

输出：
  output_prefix_pass.bam   # 保留的唯一 unmapped reads  
  output_prefix_fail.bam   # 重复的 unmapped reads
"""

import sys
import time
import argparse
import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="过滤 BAM：保留唯一的 unmapped reads，重复的写入 fail 文件"
    )
    parser.add_argument(
        "-bam", "--input", required=True,
        help="输入 BAM 文件路径"
    )
    parser.add_argument(
        "-o", "--prefix", required=True,
        help="输出文件名前缀"
    )
    return parser.parse_args()

def main():
    start_time = time.time()
    args = parse_args()

    in_bam_path = args.input
    prefix = args.prefix
    out_pass = f"{prefix}_pass.bam"
    out_fail = f"{prefix}_fail.bam"

    print(f"[INFO] 打开输入 BAM：{in_bam_path}")
    bam_in = pysam.AlignmentFile(in_bam_path, "rb")
    bam_pass = pysam.AlignmentFile(out_pass, "wb", header=bam_in.header)
    bam_fail = pysam.AlignmentFile(out_fail, "wb", header=bam_in.header)

    seen_seqs = set()
    total = 0
    kept = 0
    cnt_dup = 0

    print("[INFO] 开始逐条处理 reads ...")
    for rec in bam_in.fetch(until_eof=True):
        total += 1

        # 跳过所有已比对的 reads（mapped）
        if not rec.is_unmapped:
            continue

        seq = rec.query_sequence
        if seq in seen_seqs:
            # 序列重复 -> 写入 fail
            bam_fail.write(rec)
            cnt_dup += 1
        else:
            # 首次出现 -> 写入 pass
            seen_seqs.add(seq)
            bam_pass.write(rec)
            kept += 1

        # 每处理 1,000,000 条输出一次进度
        if total % 1_000_000 == 0:
            print(f"[INFO] 已处理 {total:,} 条 reads，保留 {kept:,} 条，重复 {cnt_dup:,} 条")

    # 关闭文件句柄
    bam_in.close()
    bam_pass.close()
    bam_fail.close()

    elapsed = time.time() - start_time
    print(f"[INFO] 处理完毕，用时 {elapsed:.2f} 秒")
    print(f"[INFO] 总读取:       {total:,} 条")
    print(f"[INFO] 保留 (pass):  {kept:,} 条 -> {out_pass}")
    print(f"[INFO] 重复 (fail):  {cnt_dup:,} 条 -> {out_fail}")

if __name__ == "__main__":
    main()
