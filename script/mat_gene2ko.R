#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：按丰度筛选OTU，再汇总为各分类级文件，作为STAMP和lefse输入
# Functions: Boxplot show alpha richness among groups
# Main steps: 
# - Reads data table input.txt
# - Calculate pvalue and save in output.txt
# - Draw boxplot and save in output.pdf

# 程序使用示例
# USAGE
# Default
# # 显示脚本帮助
# Rscript ~/github/Metagenome/denovo1/script/mat_gene2ko.R -h
# Rscript ~/github/Metagenome/denovo1/script/mat_gene2ko.R -i genetab.txt -o kotab.txt

options(warn = -1) # Turn off warning


# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="temp/24eggnog/gene_anno.count",
                help="OTU [default %default]"),
    make_option(c("-n", "--norm"), type="numeric", default=1000000,
                help="OTU [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/24eggnog/annotab",
                help="output directory or prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  # 显示输入输出确认是否正确
  print(paste("Input gene table: ", opts$input,  sep = ""))
  print(paste("Normalization to : ", opts$norm, sep = ""))
  print(paste("Output KO table: ", opts$output, sep = ""))
}
# 合并每个级别分类学表

# 0. 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 1. 读取gene表
gene = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)

# 2. 按KO列求合
KO = gene %>% group_by(KO) %>% summarise_all(sum)

# 整理为数值矩阵
# 把第一列变为列名，纯数值才能操作
df = as.data.frame(KO)
row.names(df) = df$KO
df = df[,-1]
# colSums(df)
# 四余五入取整，ceiling()  向上取整, floor() 向下取整
df = floor(df+0.5)
colSums(df)

# 保存KO表
write.table("ID\t", file=paste(opts$output, ".count", sep=""), append = F, sep="\t", eol = "", quote=F, row.names=F, col.names=F)
write.table(df, file=paste(opts$output, ".count", sep = ""), append = T, sep="\t", quote=F, row.names=T, col.names=T)

# 3. 标准化为1、100、1000000
norm = round(t(t(df)/colSums(df,na=T)) * opts$norm, 2)
colSums(norm)

write.table("ID\t", file=paste(opts$output, ".tpm", sep=""), append = F, sep="\t", eol = "", quote=F, row.names=F, col.names=F)
write.table(norm, file=paste(opts$output, ".tpm", sep = ""), append = T, sep="\t", quote=F, row.names=T, col.names=T)

