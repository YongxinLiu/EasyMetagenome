#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. 
# A practical guide to amplicon and metagenomic analysis of microbiome data. 
# Protein Cell 41, 1-16, doi:10.1007/s13238-020-00724-8 (2020).
# Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
# NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. 
# Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：Beta多样性分析——限制性主坐标分析
# Functions: Beta diversity——Constrained PCoA

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为距离矩阵(如bray_curtis.txt)+分组信息(metadata.tsv)
#
# 输入文件"-i", "--input"，result/beta/bray_curtis.txt; beta多样性距离矩阵文件，有多种距离可选，常用bray_curtis、jaccard、unifrac和unifrac_binary
#
# 实验设计"-d", "--design"，默认`metadata.tsv`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.tsv中的Group列作为分组信息，可修改为任意列名；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小
#
# 分组列名"-o", "--output"，输出目录，默认同输入+.cpcoa.pdf；
#


#----1.2 参数缺少值 Default values#----
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/beta/bray_curtis.txt",
                help="Alpha rarefaction richness [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=59,
                help="Figure heidth in mm [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output==""){opts$output=paste0(opts$input,".cpcoa.pdf")}
opts$output = gsub(".txt","",opts$output)
suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))


#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

#----2.1 距离矩阵Distance matrix#----
distance_mat = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")



#----3. 绘图保存 Plotting and saving#----

#----3.1 绘图 Plotting#----
# 输入矩阵矩阵、元数据和分组列，返回ggplot2对象
p = beta_cpcoa_dis(distance_mat, metadata, groupID = opts$group)

#---3.2 保存 Saving#----
# 大家可以修改图片名称和位置，长宽单位为毫米
ggsave(opts$output, p, width = opts$width, height = opts$height, units = "mm")
