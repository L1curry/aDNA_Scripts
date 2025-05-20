#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage: Bam_extractMD_analysis.py -i input.txt -o output_prefix [-ld NEAR_LENGTH]

This script analyzes an input text file where each line is formatted like:
Read LH00708:99:22M32KLT4:5:2279:34323:9618 (length 41) MD: 41

It extracts the read length (after "length") and the MD string (after "MD:"), then computes statistics for each read length:
  - Number of reads and percentage.
  - Base mismatch events, insertion events and new 'indel0' events based on the MD field.
    * If the MD string consists only of digits, the read is considered fully matched.
    * A single letter (A/T/C/G) indicates a mismatch (count 1 event).
    * A deletion block starting with '^' (e.g. "^T") is counted as an insertion event (1 event per base in the block).
    * Additionally, if a number "0" is encountered and its preceding token was a single letter (i.e. a mismatch event), then an extra event is added with type 'indel0'.
    * If unexpected characters are encountered, the script prints the offending line and stops.
If the optional parameter -ld is provided (e.g. 10), the first and last ld bases of each read are defined as the proximal region, and additional statistics for mismatches/insertions in the proximal region are computed.
Two output files are produced:
  1. A text file (output_prefix + "_stats.txt") with the detailed statistics.
  2. Two bar charts in PDF format:
     - One chart shows the overall distribution of events by type (e.g. mis_A, ins_T, indel0, etc.).
     - A second chart shows the distribution of total error counts per read (errors are aggregated, where a deletion block or an indel0 event is counted as 1 error).
     
Example:
  python Bam_extractMD_analysis.py -i input.txt -o result -ld 10
"""

import argparse
import sys
import re
from collections import defaultdict, Counter

import matplotlib.pyplot as plt
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze read length and MD field from an input text file to compute mismatch/insertion statistics.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", required=True, type=str,
                        help="Input text file; each line formatted like:\nRead LH00708:99:22M32KLT4:5:2279:34323:9618 (length 41) MD: 41")
    parser.add_argument("-o", required=True, type=str,
                        help="Output file prefix. Statistics will be saved to {prefix}_stats.txt; bar charts will be saved as PDF.")
    parser.add_argument("-ld", type=int, default=None,
                        help="Optional. If specified (e.g. 10), defines the proximal region as the first and last ld bases.")
    return parser.parse_args()

def parse_md(md_str, read_length):
    """
    解析 MD 字符串，返回事件列表。
    每个事件为元组 (pos, event_type, base)
      pos: 在read上的位置（1-indexed）
      event_type: "mismatch" 表示错配，"insertion" 表示删除块产生的插入事件，
                  "indel0" 表示新规则下的事件（当数字0前面不是数字而是字母时计1）
      base: 错配或插入的碱基（对于删除块中的每个碱基均单独计数）
    说明：删除块（以'^'开头）每个字母都生成一个事件，但不消耗read碱基；数字则仅用来移动位置。
    新规则：若匹配到的数字为 "0"，且其前一个匹配的类型为单个字母（即 mismatch），则额外添加一个 indel0 事件。
    """
    events = []
    current_pos = 0
    # 定义正则表达式，匹配数字、删除块（以'^'开头）、单个字母（A/T/C/G）
    pattern = re.compile(r'(\d+)|(\^[ATCG]+)|([ATCG])')
    prev_type = None  # 记录上一个匹配的类型： "number", "deletion", "mismatch"
    prev_text = None  # 记录上一个匹配的文本
    for match in pattern.finditer(md_str):
        if match.group(1):  # 数字部分
            token = match.group(1)
            # 如果数字为 "0" 且上一个token类型为 "mismatch"，则新增一个 indel0 事件
            if token == "0" and prev_type == "mismatch":
                # 添加 indel0 事件，位置不变（current_pos 不改变）
                events.append((current_pos+1, 'indel0', '0'))
                # 不改变 current_pos
            else:
                num = int(token)
                current_pos += num
            current_token_type = "number"
        elif match.group(2):  # 删除块，如 ^T 或 ^TG
            deletion = match.group(2)
            # 对于删除块中的每个字母都计为一个 insertion 事件
            for base in deletion[1:]:
                if base not in "ATCG":
                    raise ValueError(f"Unexpected character {base} in deletion block: {deletion}")
                events.append((current_pos+1, 'insertion', base))
            # 删除块不消耗read碱基
            current_token_type = "deletion"
        elif match.group(3):  # 单个字母 mismatch
            base = match.group(3)
            events.append((current_pos+1, 'mismatch', base))
            current_pos += 1
            current_token_type = "mismatch"
        else:
            raise ValueError("Unknown parsing error in MD string.")
        # 更新 prev_type 和 prev_text
        prev_type = current_token_type
        prev_text = match.group(0)
    if current_pos > read_length:
        raise ValueError(f"Parsed position {current_pos} exceeds read length {read_length}.")
    return events

def compute_error_count(md_str):
    """
    计算 MD 字符串中的聚合错误数。
    - 数字不计入错误。
    - 删除块（例如 '^T' 或 '^TG'）计为1错误（不论长度）。
    - 单个字母（A/T/C/G）计为1错误。
    - 对于新规则出现的数字"0"（如果前面是字母），也计为1错误。
    说明：由于 aggregated error count 的计算不需要区分详细类型，因此这里与之前保持一致，
    我们仍使用正则表达式忽略数字。
    """
    count = 0
    pattern = re.compile(r'(\d+)|(\^[ATCG]+)|([ATCG])')
    for match in pattern.finditer(md_str):
        if match.group(1):
            token = match.group(1)
            # 对于"0"不计入错误，因为新规则已在详细统计中处理；聚合错误不改变
            continue
        elif match.group(2):
            count += 1
        elif match.group(3):
            count += 1
    return count

def main():
    args = parse_args()
    
    input_file = args.i
    output_prefix = args.o
    near_len = args.ld  # 如果提供，则为近端区长度
    
    total_reads = 0
    # 按read长度统计信息
    length_stats = defaultdict(lambda: {
        'count': 0,
        'mismatches': Counter(),         # 详细错配
        'insertions': Counter(),          # 详细插入（删除块产生的事件）
        'indel0': Counter(),              # 新规则下的事件（数字0前面为字母）
        'near_end_mismatches': Counter(), # 近端区错配
        'near_end_insertions': Counter(), # 近端区插入
        'near_end_indel0': Counter()      # 近端区新规则事件
    })
    # 全局详细事件统计（用于第一张柱状图）
    overall_mismatches = Counter()
    overall_insertions = Counter()
    overall_indel0 = Counter()
    # 全局聚合错误计数（用于第二张柱状图）
    overall_err_counts = Counter()
    
    try:
        infile = open(input_file, "r")
    except Exception as e:
        sys.exit(f"Error: Unable to open input file: {e}")
    
    for line in infile:
        line = line.strip()
        if not line:
            continue
        total_reads += 1
        
        # 提取 read 长度，例如 "(length 41)"
        len_match = re.search(r'\(length\s+(\d+)\)', line)
        if not len_match:
            sys.exit(f"Error: Cannot find read length in line: {line}")
        try:
            read_length = int(len_match.group(1))
        except ValueError:
            sys.exit(f"Error: Read length is not an integer in line: {line}")
        
        # 提取 MD 字段：查找 "MD:" 后的非空字符
        md_match = re.search(r'MD:\s*(\S+)', line)
        if not md_match:
            sys.exit(f"Error: Cannot find MD field in line: {line}")
        md_str = md_match.group(1)
        
        # 更新该 read 长度下的计数
        length_stats[read_length]['count'] += 1
        
        # 如果 MD 字符串仅由数字组成，则认为该 read 完全匹配，错误数为 0
        if re.fullmatch(r'\d+', md_str):
            overall_err_counts[0] += 1
            continue
        
        # 解析 MD 字符串获取详细事件列表；若出错则停止
        try:
            events = parse_md(md_str, read_length)
        except ValueError as ve:
            sys.exit(f"Error while parsing MD field in line: {line}\n{ve}")
        
        # 计算该 read 的聚合错误数（不区分详细类型）
        err_count = compute_error_count(md_str)
        overall_err_counts[err_count] += 1
        
        # 根据解析得到的事件更新统计信息
        for pos, event_type, base in events:
            if event_type == 'mismatch':
                length_stats[read_length]['mismatches'][base] += 1
                overall_mismatches[base] += 1
                if near_len is not None:
                    if pos <= near_len or pos > (read_length - near_len):
                        length_stats[read_length]['near_end_mismatches'][base] += 1
            elif event_type == 'insertion':
                length_stats[read_length]['insertions'][base] += 1
                overall_insertions[base] += 1
                if near_len is not None:
                    if pos <= near_len or pos > (read_length - near_len):
                        length_stats[read_length]['near_end_insertions'][base] += 1
            elif event_type == 'indel0':
                length_stats[read_length]['indel0'][base] += 1
                overall_indel0[base] += 1
                if near_len is not None:
                    if pos <= near_len or pos > (read_length - near_len):
                        length_stats[read_length]['near_end_indel0'][base] += 1
            else:
                sys.exit(f"Error: Unknown event type {event_type} in line: {line}")
    
    infile.close()
    
    # 输出统计结果到文本文件
    stats_output = f"{output_prefix}_stats.txt"
    try:
        outf = open(stats_output, "w")
    except Exception as e:
        sys.exit(f"Error: Unable to create output stats file: {e}")
    
    # 输出表头（使用英文）
    header = "Read_Length\tCount\tPercent\tMismatch_Details\tInsertion_Details\tIndel0_Details"
    if near_len is not None:
        header += "\tProximal_Mismatches\tProximal_Insertions\tProximal_Indel0"
    outf.write(header + "\n")
    
    total_count = sum(item['count'] for item in length_stats.values())
    for rlen in sorted(length_stats.keys()):
        count = length_stats[rlen]['count']
        percent = (count / total_count * 100) if total_count > 0 else 0
        mismatch_detail = ", ".join([f"{b}:{c}" for b, c in length_stats[rlen]['mismatches'].items()]) if length_stats[rlen]['mismatches'] else "None"
        insertion_detail = ", ".join([f"{b}:{c}" for b, c in length_stats[rlen]['insertions'].items()]) if length_stats[rlen]['insertions'] else "None"
        indel0_detail = ", ".join([f"{b}:{c}" for b, c in length_stats[rlen]['indel0'].items()]) if length_stats[rlen]['indel0'] else "None"
        line_out = f"{rlen}\t{count}\t{percent:.2f}%\t{mismatch_detail}\t{insertion_detail}\t{indel0_detail}"
        if near_len is not None:
            near_mm = ", ".join([f"{b}:{c}" for b, c in length_stats[rlen]['near_end_mismatches'].items()]) if length_stats[rlen]['near_end_mismatches'] else "None"
            near_ins = ", ".join([f"{b}:{c}" for b, c in length_stats[rlen]['near_end_insertions'].items()]) if length_stats[rlen]['near_end_insertions'] else "None"
            near_indel0 = ", ".join([f"{b}:{c}" for b, c in length_stats[rlen]['near_end_indel0'].items()]) if length_stats[rlen]['near_end_indel0'] else "None"
            line_out += f"\t{near_mm}\t{near_ins}\t{near_indel0}"
        outf.write(line_out + "\n")
    outf.close()
    print(f"Statistics saved to {stats_output}")
    
    # -----------------------
    # 第一张柱状图：总体详细事件分布
    # 合并错配、插入和 indel0 的统计，标签前缀分别为 mis_、ins_、indel0_
    combined = Counter()
    for base, cnt in overall_mismatches.items():
        combined[f"mis_{base}"] = cnt
    for base, cnt in overall_insertions.items():
        combined[f"ins_{base}"] = cnt
    for base, cnt in overall_indel0.items():
        combined[f"indel0_{base}"] = cnt
    if combined:
        event_types = list(combined.keys())
        counts = [combined[et] for et in event_types]
        plt.figure(figsize=(10,6))
        bars = plt.bar(event_types, counts, color='skyblue', edgecolor='black')
        plt.xlabel("Event Type")
        plt.ylabel("Count")
        plt.title("Overall Mismatch/Insertion/Indel0 Distribution")
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, height, f'{int(height)}', ha='center', va='bottom')
        plt.tight_layout()
        bar_output = f"{output_prefix}_bar.pdf"
        plt.savefig(bar_output)
        plt.close()
        print(f"Detailed event bar chart saved to {bar_output}")
    else:
        print("No mismatch/insertion/indel0 events detected; no detailed event bar chart generated.")
    
    # -----------------------
    # 第二张柱状图：每个 read 聚合错误计数的分布
    if overall_err_counts:
        err_counts = sorted(overall_err_counts.keys())
        frequencies = [overall_err_counts[ec] for ec in err_counts]
        plt.figure(figsize=(10,6))
        bars = plt.bar([str(ec) for ec in err_counts], frequencies, color='lightgreen', edgecolor='black')
        plt.xlabel("Error Count per Read")
        plt.ylabel("Number of Reads")
        plt.title("Distribution of Aggregated Error Counts per Read")
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, height, f'{int(height)}', ha='center', va='bottom')
        plt.tight_layout()
        err_bar_output = f"{output_prefix}_err_bar.pdf"
        plt.savefig(err_bar_output)
        plt.close()
        print(f"Aggregated error count bar chart saved to {err_bar_output}")
    else:
        print("No error events detected; no aggregated error count bar chart generated.")

if __name__ == "__main__":
    main()
