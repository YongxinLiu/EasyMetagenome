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

# 1.1 程序功能描述和主要步骤

# 程序功能：Alpha多样性柱状图+统计
# Functions: Alpha boxplot

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为原始alpha多样性(vegan.txt)+分组信息(metadata.tsv)
#
# 输入文件"-i", "--input"，result/alpha/vegan.txt; alpha多样性表格
#
# 实验设计"-d", "--design"，默认`metadata.tsv`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.tsv中的Group列作为分组信息，可修改为任意列名；
#
# 分组列名"-o", "--output"，默认为输出目录，图片文件名为alpha_boxplot_+多样性指数名+.pdf；统计文本位于代码运行目录中alpha_boxplot_TukeyHSD.txt；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


# 1.2 解析命令行
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
    make_option(c("-i", "--input"), type="character", default="result/alpha/vegan.txt",
                help="Alpha diversity matrix [default %default]"),
    make_option(c("-a", "--alpha_index"), type="character", default="richness",
                help="Group name [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.tsv",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/alpha/",
                help="Output pdf directory, with prefix alpha_boxplot_; Stat in alpha_boxplot_TukeyHSD.txt [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=59,
                help="Figure heigth in mm [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  suppressWarnings(dir.create(opts$output))
}

if(opts$output==""){opts$output="./"}
suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


#---- 2. 依赖关系检查、安装和加载 ----

suppressWarnings(suppressMessages(library(amplicon)))

# 3. 读取输入文件

# 读取OTU表
alpha_div = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")

# 读取实验设计
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

# 4. 绘图

p = alpha_barplot(alpha_div, metadata, index = opts$alpha_index, groupID = opts$group)

# 5. 保存图片Saving figure
# 保存图片，大家可以修改图片名称和位置，长宽单位为毫米
ggsave(paste0(opts$output,"/alpha_barplot_",opts$alpha_index,".pdf"), p, width = opts$width, height = opts$height, units = "mm")
