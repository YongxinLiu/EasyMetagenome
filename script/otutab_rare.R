#!/usr/bin/env Rscript

# Copyright 2016-2021 Yong-Xin Liu <yxliu@genetics.ac.cn / metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. (2021). A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12, doi: https://doi.org/10.1007/s13238-020-00724-8


# 1. 分析前准备：帮助、参数、依赖包和读取文件

# 命令行运行为当前目录；Rstudio手动运行脚本前需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：OTU表抽平，计算Alpha多样性
# otutab_rare Script functions: Rarefaction and alpha diversity
# Main steps: 
# - Read data table taxonomy_count.txt
# - Select lowest taxonomy
# - Rarefaction to same depth
# - Alpha diversity calculation
# - Save rare table and alpha diversity

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
              help="Input reads count file; such as OTU table, kraken2 taxonomy counts table [default %default]"),
  make_option(c("-d", "--depth"), type="numeric", default=0,
              help="Rarefaction depth; default 0 to auto using minimum [default %default]"),
  make_option(c("-s", "--seed"), type="numeric", default=1,
              help="Random sample seed, can set any integer [default %default]"),  make_option(c("-n", "--normalize"), type="character", default="result/otutab_rare.txt",
              help="Normalized filename [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result/alpha/vegan.txt",
              help="Output alpha diversity filename [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The input feature table is ", opts$input,  sep = ""))
# print(paste("Normalized filename: ", opts$normalize,  sep = ""))
# print(paste("Output alpha diversity: ", opts$output, sep = ""))

# suppressWarnings(dir.create("alpha/"))

# 1.3 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("vegan")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


# 1.4 读取输入文件

# 默认的quote会跳过2/3的数据行减少，产生NA，改为空
species = read.table(opts$input, header=T, sep="\t", quote = "", row.names=1, comment.char="") 

# 2. 计算过程

## 2.1 抽样至最少值或指定值Rarefaction
print(paste0("Samples size are:"))
colSums(species)
min = min(colSums(species))
if (opts$depth==0){
  opts$depth=min}
print(paste("Rarefaction depth ", opts$depth, ". If depth set 0 will using sample minimum size ", min, sep = ""))
# 抽样数为零则使用最小值，否则使用指定抽样数

# print(paste0("Rarefaction depth is ", min))
# vegan::rrarefy抽平至最小值或指定值
print(paste("Random sample number: ", opts$seed,  sep = ""))
set.seed(opts$seed)
otu = vegan::rrarefy(t(species), opts$depth)
# print(paste0("All sample rarefaction as following"))
# rowSums(otu)
# 筛选大于抽样数量的样本
idx = rowSums(otu) >= opts$depth
# 保存丢弃的样本列表
print(paste0("The discard samples: ",rownames(otu[!idx,])))
suppressWarnings(write.table(rownames(otu[!idx,]), file=paste(opts$normalize,".discard",sep=""), append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F))

# 筛选大于抽样数量的样本的特征表
otu = otu[idx,]


## 2.3 Alpha diversity
# vegan::estimateR计算obs, chao1和ACE指数
estimateR = t(estimateR(otu))[,c(1,2,4)]
colnames(estimateR) = c("richness", "chao1", "ACE")
# vegan::diversity计算多样性指数shannon, simpson和invsimpson
shannon = diversity(otu, index = "shannon")
simpson = diversity(otu, index = "simpson")
invsimpson = diversity(otu, index = "invsimpson")
# 合并6种指数
alpha_div = cbind(estimateR, shannon, simpson, invsimpson)
print(paste0("Calculate six alpha diversities by estimateR and diversity"))
head(alpha_div, n=1)

# 3. 结果输出
# 3.1 保存抽平的物种表
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("#OTUID\t", file=paste(opts$normalize,sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(t(otu), file=paste(opts$normalize,sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

# 3.2 保存alpha多样性指数
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("SampleID\t", file=paste(opts$output,sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(alpha_div, file=paste(opts$output,sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

print(paste("Name of rarefaction file ", opts$normalize,  sep = ""))
print(paste("Output alpha diversity filename ", opts$output, sep = ""))
print(paste("The discard samples file ",  opts$normalize, ".discard", sep = ""))
