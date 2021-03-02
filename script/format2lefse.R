#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).


# 1. 分析前准备：帮助、参数、依赖包和读取文件

# 命令行运行为当前目录；Rstudio手动运行脚本前需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：转换OTU表为LEfSe输入文件
# otutab_rare Script functions: Format OTU table, taxonomy and metadata to LEfSe input file
# Main steps:
# - Read data table otutab.txt
# - Read data table taxonomy.txt
# - Read data table metadata.txt
# - Select column name of metadata
# - Abundance threshold
# - Save LEfSe input in lefse.txt

# 程序使用示例
# USAGE
# Default
# # 显示脚本帮助
# Rscript ./script/otutab_rare.R -h
# # 默认读取result/otutab.txt，按最小值抽平，输出输入文件前缀的标准化表和alpha多样性
# Rscript ./script/otutab_rare.R
# # 完整参数，输出文件名默认为alpha指数类型
# Rscript ./script/otutab_rare.R -i result/otutab.txt \
# -d 10000 \
# -n result/otutab_rare.txt \
# -o alpha/vegan_diversity.txt
# options(warn = -1) # Turn off warning



# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
              help="Input reads count file; such as OTU table [default %default]"),
  make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt",
              help="Input taxonomy file [default %default]"),
  make_option(c("-d", "--design"), type="character", default="result/metadata.tsv",
              help="Experiment design or sample metadata [default %default]"),
  make_option(c("-T", "--threshold"), type="numeric", default=0.1,
              help="Threshold of abundance percentage, such as 0.1% [default %default]"),
  make_option(c("-n", "--group"), type="character", default="Group",
              help="Column name of group [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result/LEfSe.txt",
              help="Output directory and filename [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The input feature table is ", opts$input,  sep = ""))
print(paste("Input taxonomy file: ", opts$taxonomy,  sep = ""))
print(paste("Experiment design or sample metadata: ", opts$design, sep = ""))
print(paste("Threshold of abundance percentage: ", opts$threshold,  sep = ""))
print(paste("Column name of group: ", opts$group,  sep = ""))
print(paste("Output directory and filename: ", opts$output,  sep = ""))


# 1. 读取OTU表
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# 2. 读取物种注释
tax = read.table(opts$taxonomy, header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F)
# 3. 读取样本分组信息
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="")
# OTU丰度筛选阈值，默认0.1%，主要为圈图展示合适数据的OTU
# thre = 0.1
# 样本分组列名
# groupID = "Group"

# 生成各分类级汇总特征表
suppressWarnings(suppressMessages(library(amplicon)))
format2lefse(otutab, tax, metadata, opts$threshold, opts$group, paste0(opts$output,".txt"))
format2lefse2(otutab, tax, metadata, opts$threshold, opts$group, paste0(opts$output,"2.txt"))
# 在当前目录生成文件opts$output
