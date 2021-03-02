#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).


#----1. 参数 Parameters#----

# 命令行运行为当前目录；Rstudio手动运行脚本前需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1.1 功能描述 Function description#----

# 程序功能：转换OTU表为STAMP输入文件
# Functions of format2stamp script: Format OTU table and taxonomy to STAMP input file
# Main steps:
# - Read data table otutab.txt
# - Read data table taxonomy.txt
# - Abundance threshold
# - Save STAMP input of each taxonomy levles in output directory

# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}


#----1.2 参数缺省值 Default values#----

## 设置输入输出文件和参数
# 修改下面`default=`后面的文件和参数。
# 输入文件"-i", "--input"，默认result/otutab.txt; 特征表(OTU/ASV)
# 实验设计"-t", "--taxonomy"，默认result/taxonomy.txt，包括特征的界门纲目科属种注释
# 分组列名"-T", "---threshold"，默认0.1%，结果只保留此丰度以上的特征
# 分组列名"-o", "--output"，默认为输出目录及文件名前缀；
# 解析参数-h显示帮助信息

option_list <- list(
  make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
              help="Input reads count OTU/ASV table [default %default]"),
  make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt",
              help="Input taxonomy in 7 levels [default %default]"),
  make_option(c("-T", "--threshold"), type="numeric", default=0.1,
              help="Threshold of relative abundance, such as 0.1% [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result/stamp/tax",
              help="Output directory and prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The input feature table is ", opts$input,  sep = ""))
print(paste("Input taxonomy file: ", opts$taxonomy,  sep = ""))
print(paste("Threshold of abundance percentage: ", opts$threshold,  sep = ""))
print(paste("Output directory and filename prefix: ", opts$output,  sep = ""))

#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))

#----2. 读取文件 Read files#----

#----2.1 特征表 OTU表#----
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
#----2.2 物种组成矩阵Taxonomy matrix#----
tax = read.table(opts$taxonomy, header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F)
# OTU丰度筛选阈值，默认0.1%，主要为圈图展示合适数据的OTU
# thre = 0.1

#----3. 分析保存 Analysis and saving#----

# 生成各分类级汇总特征表
format2stamp(otutab = otutab, taxonomy = tax, thre = opts$threshold, prefix = opts$output)
# 在输出目录生成tax_1-8共7个级别+OTU过滤文件
