#!/usr/bin/env Rscript

# Copyright 2016-2020, Tao Wen,  Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# 1. EasyAmplicon(TBD)
# 2. Jingying Zhang, Yong-Xin Liu, et.al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：筛选OTU表中N行，根据物种注释绘制网络气泡图
# Functions: Plot maptree based on OTU table and taxonomy file

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为特征表(otutab.txt)+物种注释(taxonomy.txt)
#
# 输入文件"-i", "--input"，otutab.txt; 特征表
#
# 物种注释"-t", "--taxonomy"，taxonomy.txt; 物种注释
#
# 图片宽"-w", "--width"，默认183 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认118 mm，根据图像布局可适当增大或缩小


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
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="OTU table [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt",
                help="Taxonomy in 8 fields [default %default]"),
    make_option(c("-N", "--topN"), type="numeric", default=200,
                help="Top N OTU or ASV [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/tax/tax_maptree.pdf",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=183,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=118,
                help="Figure heidth in mm [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}


# 2. 依赖关系检查、安装和加载

suppressWarnings(suppressMessages(library(amplicon)))


# 3. 读取输入文件

# 读取OTU表
otutab = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")
otutab = as.matrix(otutab)
# 读取物种组成
taxonomy = read.table(opts$taxonomy, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
taxonomy = as.matrix(taxonomy)

# 数据转换为maptree格式
mapdata = format2maptree(otutab, taxonomy, opts$topN)
#按照平均丰度修改大小和按照门水平上颜色
mapadd = tax_maptree(mapdata)
(p = mapadd[[1]])

# Saving figure
# 保存图片，大家可以修改图片名称和位置，长宽单位为毫米
ggsave(opts$output, p, width = opts$width, height = opts$height, units = "mm")
