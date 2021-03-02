#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：OTU按实验设计分组，计算菌值
# Functions: Calculate mean OTU abundance for tree
# Main steps: 
# - Reads data table input.txt


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
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="Feature table [default %default]"),
    make_option(c("-T", "--thre"), type="numeric", default=0,
                help="Threshold of abundance, such as 0.001 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.tsv",
                help="Experiment design or sample metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Column name of group [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/otutab_mean.txt",
                help="output directory and prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))

  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$input, "_", opts$group, sep = "")}

  # 显示输入输出确认是否正确
  print(paste("Feature table: ", opts$input,  sep = ""))
  print(paste("Metadata: ", opts$design,  sep = ""))
  print(paste("Group name: ", opts$group,  sep = ""))
  print(paste("Abundance threshold: ", opts$thre,  sep = ""))
  print(paste("Output filename: ", opts$output, sep = ""))
}

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

# 1. 读取OTU表
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# head(otutab)

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F) 

# 将选定的分组列统一命名为group
design$group=design[,opts$group]

# 标准化，并筛选高丰度菌均值最小百万分之一0.0001%
norm = t(otutab)/colSums(otutab,na=T)*100
# rowSums(norm)
idx = colMeans(norm) > opts$thre
HA = norm[,idx]
# dim(HA)
# rowSums(HA)

# 按group合并
merge=cbind(HA, design[,c("group"),drop=F])
HA_group_mean = merge %>% group_by(group) %>% summarise_all(mean)
HA_t = as.data.frame(cbind(c("All", round(colMeans(HA), digits = 6)),t(HA_group_mean)), stringsAsFactors = F)
rownames(HA_t)[1] = "OTUID"
write.table(HA_t, file=paste(opts$output, "", sep = ""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=F)

