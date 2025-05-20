#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage: Bam_extractMD.py -i input.bam -o output_prefix -length <LOWER> <UPPER>

该脚本从 BAM 文件中提取 MD 标签。用户通过 -length 参数指定读长筛选条件：
  - 若参数为 "min 50"：提取长度 <= 50 bp 的记录；
  - 若参数为 "51 100"：提取长度 > 50 bp 且 <= 100 bp 的记录；
  - 若参数为 "101 max"：提取长度 >= 101 bp 的记录；
  - 若参数为 "min max"：不考虑读长，提取所有 reads 的 MD 标签。
  
示例：
  python Bam_extractMD.py -i sample.bam -o result -length min 50
  python Bam_extractMD.py -i sample.bam -o result -length min max
"""

import argparse
import sys
import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="该脚本根据指定的读长区间提取 BAM 文件中比对记录的 MD 标签。",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", required=True, type=str,
                        help="输入的 BAM 文件。")
    parser.add_argument("-o", required=True, type=str,
                        help="输出文件前缀。提取的 MD 标签将保存在 {prefix}_MD.txt 中。")
    parser.add_argument("-length", nargs=2, required=True, metavar=('LOWER', 'UPPER'),
                        help=("指定读长区间，必须提供两个值：\n"
                              "  'min 50'   表示查找长度 <= 50 bp 的记录\n"
                              "  '51 100'   表示查找长度 > 50 bp 且 <= 100 bp 的记录\n"
                              "  '101 max'  表示查找长度 >= 101 bp 的记录\n"
                              "  'min max'  表示不限制读长，提取所有 reads"))
    return parser.parse_args()

def main():
    args = parse_args()
    
    input_bam = args.i
    output_prefix = args.o
    lower_val, upper_val = args.length  # 例如：lower_val = "min", upper_val = "50" 或 "max"
    
    # 当 -length 参数为 "min max" 时，不做读长筛选，提取所有 reads 的 MD 标签
    if lower_val.lower() == "min" and upper_val.lower() == "max":
        def length_condition(l):
            return True
    # 如果只提供了 "min" 开头，则要求第二个值为整数，筛选读长 <= upper_limit
    elif lower_val.lower() == "min":
        try:
            upper_limit = int(upper_val)
        except ValueError:
            sys.exit("错误：当第一个值为 'min' 时，第二个值必须为整数或 'max'。")
        def length_condition(l):
            return l <= upper_limit
    # 如果只提供了 "max" 结尾，则要求第一个值为整数，筛选读长 >= lower_limit
    elif upper_val.lower() == "max":
        try:
            lower_limit = int(lower_val)
        except ValueError:
            sys.exit("错误：当第二个值为 'max' 时，第一个值必须为整数或 'min'。")
        def length_condition(l):
            return l >= lower_limit
    else:
        try:
            lower_limit = int(lower_val)
            upper_limit = int(upper_val)
        except ValueError:
            sys.exit("错误：-length 参数的两个值必须为 'min'/'max' 或整数。")
        def length_condition(l):
            return (l > (lower_limit - 1)) and (l <= upper_limit)
    
    # 打开 BAM 文件（只读模式）
    try:
        bamfile = pysam.AlignmentFile(input_bam, "rb")
    except Exception as e:
        sys.exit(f"错误：无法打开 BAM 文件: {e}")
    
    # 打开输出文件，保存提取的 MD 标签
    output_file = f"{output_prefix}_MD.txt"
    try:
        outfile = open(output_file, "w")
    except Exception as e:
        sys.exit(f"错误：无法创建输出文件: {e}")
    
    # 遍历 BAM 文件中所有比对记录（直接迭代，无需索引）
    for read in bamfile:
        qlen = read.query_length
        if qlen is None:
            continue  # 若无法获取序列长度，则跳过
        if length_condition(qlen):
            try:
                md_tag = read.get_tag("MD")
                outfile.write(f"Read {read.query_name} (length {qlen}) MD: {md_tag}\n")
            except KeyError:
                outfile.write(f"Read {read.query_name} (length {qlen}) 无 MD 标签\n")
    bamfile.close()
    outfile.close()
    print(f"MD标签提取完成，结果保存到 {output_file}")

if __name__ == "__main__":
    main()

