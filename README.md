# 此项目脚本用于 aDNA 群体遗传分析
----------------------------------------------------------------------------------
# Bam_extractMD.py

## 简介
`Bam_extractMD.py` 是一个用于从 BAM 文件中批量提取 MD 标签的小工具。通过自定义 read 长度区间，快速筛选并获取所需的 MD 信息，方便下游分析。

## 功能
- 遍历 BAM 文件中的每条比对记录
- 根据指定的 read 长度区间进行筛选
- 提取并输出 MD 标签（无 MD 标签时给出提示）
- 无需索引即可直接处理 BAM 文件

## 所需 Python 包
- Python 3.6+
- pysam

安装命令：
pip install pysam


## 示例命令

```
# 提取所有长度 ≤ 50 bp 的 MD 标签
python Bam_extractMD.py -i sample.bam -o result -length min 50

# 提取长度在 51-100 bp 之间的 MD 标签
python Bam_extractMD.py -i sample.bam -o result -length 51 100

# 提取长度 ≥ 101 bp 的 MD 标签
python Bam_extractMD.py -i sample.bam -o result -length 101 max

# 不限长度，提取所有 MD 标签
python Bam_extractMD.py -i sample.bam -o result -length min max
```



# Bam_analysis1.py

## 简介
`Bam_analysis1.py` 是一个基于 samtools depth 结果绘制覆盖深度图形的脚本，用于评估 BAM 文件在染色体或指定区域上的测序深度分布情况。

## 功能
- 根据 samtools depth 输出的三列数据（染色体、位置、深度）绘制不同类型的图形：
  - 折线图（--zxt）：展示深度随基因组位置的变化
  - 直方图（--zft）：展示深度值的频数分布
  - 累计分布图（--ljfbt）：展示深度值的累积百分比分布
- 支持按染色体筛选：传入 "all" 绘制所有染色体，或传入整数选择第 n 条染色体
- 支持窗口范围参数（--window START END），只绘制指定区间的深度数据
- 自动处理输入文件格式校验，错误提示友好

## 所需 Python 包
- Python 3.6+
- pandas
- numpy
- matplotlib

安装命令：

pip install pandas numpy matplotlib

## 示例命令
```bash
# 对 depth.txt 在所有染色体上绘制深度直方图，输出前缀为 result
python Bam_analysis1.py --i depth.txt --zft --chr all --o result

# 在染色体 1 的位置 10 到 1000 范围内绘制深度折线图
python Bam_analysis1.py --i depth.txt --zxt --chr 1 --window 10 1000 --o result

# 对所有染色体生成累计分布图
python Bam_analysis1.py --i depth.txt --ljfbt --chr all --o result

# 同时生成折线图和直方图
python Bam_analysis1.py --i depth.txt --zxt --zft --chr all --o result
```



#  Bam_Maptype.py

## 简介
` Bam_Maptype.py` 是一个基于多进程加速的脚本，用于按染色体统计 BAM 文件中不同类型的比对记录情况，包括重复序列、低质量序列及其组合，并输出综合统计结果。

## 功能
- 统计各染色体上：
  - 总序列数量（Total）
  - PCR 重复序列数量（Duplicates）
  - 低质量比对序列数量（LowQual，mapping quality < 指定阈值）
  - 同时重复且低质量的序列数量（Dup_LowQual）
  - 保留序列数量（Retained，既非重复且质量达标）
- 支持自定义低质量阈值（`-qual`）
- 支持自定义进程数（`-proc`），默认使用系统 CPU 核心数
- 按染色体自然排序输出结果，脚本第一行为运行命令注释，后续按染色体逐行输出统计

## 所需 Python 包
- Python 3.6+
- pysam

安装命令：
pip install pysam

## 示例命令
```bash
python3  Bam_Maptype.py -bam input.bam -o output_prefix -qual 20 -proc 4
```



# Bam_extractMD_analysis.py

## 简介
`Bam_extractMD_analysis.py` 用于解析输入格式为：Read ... (length N) MD: \<MD字符串> 的文本文件，提取每条 reads 的长度和 MD 字段，并基于 MD 信息统计错配（mismatch）、插入（insertion）和新规则下的 indel0 事件。同时可选计算近端（首尾）Ld 基础上的统计。

## 功能
- 从每行中提取 read 长度和 MD 字符串
- 统计各 read 长度下：
  - reads 数量及百分比
  - 错配事件（mismatch）
  - 插入事件（insertion，基于 MD 中删除块）
  - indel0 事件（新规则下“0”前为字母时）
- 可选近端区域 (`-ld`) 内错配/插入/indel0 的统计
- 输出：
  1. `{prefix}_stats.txt` 详细统计表格
  2. `{prefix}_bar.pdf` 事件类型分布柱状图
  3. `{prefix}_err_bar.pdf` per-read 错误聚合数分布柱状图

## 所需 Python 包
- Python 3.6+
- matplotlib
- numpy

安装命令：

pip install matplotlib numpy


## 示例命令

```bash
# 基本统计并绘图
python3 Bam_extractMD_analysis.py -i input.txt -o result

# 启用近端区域统计（首尾10个碱基）
python3 Bam_extractMD_analysis.py -i input.txt -o result -ld 10
```



# SamtoolsDepth_analysis.py

## 简介
`SamtoolsDepth_analysis.py` 用于基于 samtools depth 结果文件和参考基因组染色体长度文件，识别各染色体上未测序的连续区间，并计算覆盖率（Coverage）和缺失比例（Missing(%)）。

## 功能
- 读取 samtools depth 文件（格式：`chr<TAB>pos<TAB>depth`）
- 读取染色体长度文件（格式：`chr<TAB>length`）
- 对每条染色体：
  - 提取未测序区间（连续未出现的位点）
  - 计算覆盖率：`(总长 - 缺失长度)/总长 * 100`
  - 计算缺失比例：`100 - Coverage`
- 输出包含染色体、缺失区间、Coverage(%) 和 Missing(%) 四列的结果文件

## 所需 Python 包
- Python 3.6+（无需额外第三方包）

## 示例命令
```bash
python3 SamtoolsDepth_analysis.py -depth depth.txt -chrlength chrlength.txt -o result_prefix
````



# SamtoolsDepth_analysis2.py

## 简介
`SamtoolsDepth_analysis2.py` 用于对 samtools depth 生成的深度文件进行统计分析与可视化，计算总体和各染色体的平均深度、中位数、众数、标准差、最小值及最大值，并绘制总体深度分布直方图。

## 功能
- 读取 samtools depth 文件（格式：`chrom<TAB>pos<TAB>depth`）
- 计算总体深度统计指标：平均值、 中位数、众数、标准差、最小值、最大值
- 按染色体分组计算相同统计指标，并按自然顺序排序输出
- 绘制总体深度分布直方图并保存为图片
- 输出结果文本和图像文件，文件名前缀由用户指定

## 所需 Python 包
- Python 3.6+
- pandas
- numpy
- matplotlib

安装命令：
pip install pandas numpy matplotlib

## 示例命令

```bash
python3 SamtoolsDepth_analysis2.py -i input.depth -o Test1
```



# Check_Unmapped.py

## 简介
`Check_Unmapped.py` 用于从 BAM 文件中筛选和去重 unmapped reads。脚本将首次出现的唯一序列写入 `_pass.bam`，重复的序列写入 `_fail.bam`，并输出处理统计信息。

## 功能
- 遍历输入 BAM，跳过所有已比对（mapped）的 reads
- 对 unmapped reads 根据 `query_sequence` 去重：
  - 首次出现的序列写入 `<prefix>_pass.bam`
  - 重复出现的序列写入 `<prefix>_fail.bam`
- 按每处理 1,000,000 条 reads 输出进度信息
- 处理完成后打印总计数、保留和重复条目统计及输出文件路径

## 所需 Python 包
- Python 3.6+
- pysam

安装命令：
pip install pysam


## 示例命令

```bash
python3 Check_Unmapped.py -bam input.bam -o output_prefix
```


