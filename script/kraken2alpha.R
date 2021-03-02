#!/usr/bin/env Rscript

# Copyright 2016-2021 Yong-Xin Liu <yxliu@genetics.ac.cn / metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. (2021). A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12, doi: https://doi.org/10.1007/s13238-020-00724-8


# 1. 分析前准备：帮助、参数、依赖包和读取文件

# 命令行运行为当前目录；Rstudio手动运行脚本前需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：Kraken2物种表抽平，计算Alpha多样性
# kraken2alpha Script functions: Rarefaction and alpha diversity
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
# Rscript ./script/kraken2alpha.R -h
# # 默认读取result/kraken2/mpa，按最小值抽平，输出输入文件前缀的标准化表和alpha多样性
# Rscript ./script/kraken2alpha.R
# # 完整参数，输出文件名默认为alpha指数类型
# -i kraken2物种分类汇总表
# -d 采样深度，默认为0则自动获取样本最小值
# -n 抽样标准化表文件名，空为输入文件名
# -o 输出文件名，
# Rscript ./script/kraken2alpha.R -i result/kraken2/taxonomy_count.txt \
# -d 10000 \
# -n normalized_count.txt \
# -o alpha_diversity.txt
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
    make_option(c("-i", "--input"), type="character", default="result/kraken2/tax_count.mpa",
                help="Input kraken2 reads count result [default %default]"),
    make_option(c("-d", "--depth"), type="numeric", default=0,
                help="Rarefaction depth, default minimum [default %default]"),
    make_option(c("-s", "--species"), type="character", default="",
                help="Name of species file; default input.txt [default %default]"),
    make_option(c("-n", "--normalize"), type="character", default="",
                help="Name of normalize file; default input.norm [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output alpha diversity; default input.alpha [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  prefix = gsub(".mpa$", "", opts$input, perl = T)
  if (opts$species==""){
    opts$species=paste(prefix, ".txt", sep = "")}
  if (opts$normalize==""){
    opts$normalize=paste(prefix, ".norm", sep = "")}
  if (opts$output==""){
    opts$output=paste(prefix, ".alpha", sep = "")}

  # 显示输入输出确认是否正确
  print(paste("Input file: ", opts$input,  sep = ""))
  print(paste("Rarefaction depth: ", opts$depth, ". If 0, means using sample minimum size.", sep = ""))
  print(paste("Name of normalize file ", opts$normalize,  sep = ""))
  print(paste("Output alpha diversity filename ", opts$output, sep = ""))
}


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

# 读取karken2文件
# 默认的quote会跳过2/3的数据行减少，产生NA，改为空
taxonomy = read.table(opts$input, header=T, sep="\t", quote = "", row.names=1, comment.char="") 
print(paste0("All taxonomy annotations are ", dim(taxonomy)[1], " lines!"))
# 去除NA，否则无法计算
taxonomy = na.omit(taxonomy)
# 显示样本总数据，有冗余 
colSums(taxonomy)

# 2. 计算过程

## 2.1 选择注释到种的非零、非NA行
# 2021年起kraken2 2.1.1起s__变为s_，即双下划线为单下划线
idx = grep("s_", rownames(taxonomy))
# 筛选注释到种的行,21180行变为17030行
species = taxonomy[idx, ]
# 筛选非零和缺失行,减少为1391行
idx = rowSums(species)>0
species = na.omit(species[idx,])
print(paste0("Detected non-zero species are ", dim(species)[1], "."))

## 2.2 抽样至最少值或指定值Rarefaction
print(paste0("Samples size are:"))
colSums(species)
summary(colSums(species))
min = min(colSums(species))
# 抽样数为零则使用最小值，否则使用指定抽样数
if (opts$depth==0){
  opts$depth=min}
print(paste0("Rarefaction depth is ", min))
# vegan::rrarefy抽平至最小值或指定值
set.seed(315)
otu = vegan::rrarefy(t(species), opts$depth)
# print(paste0("All sample rarefaction as following"))
# rowSums(otu)

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
write.table("Taxonomy\t", file=opts$normalize, append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(t(otu), file=opts$normalize, append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

# 3.2 保存alpha多样性指数
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("Alpha_diversity\t", file=opts$output, append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(alpha_div, file=opts$output, append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

# 3.3 保存种水平的物种表
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("Taxonomy\t", file=opts$species, append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(species, file=opts$species, append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
