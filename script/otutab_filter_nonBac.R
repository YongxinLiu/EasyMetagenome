#!/usr/bin/env Rscript

# Copyright 2016-2020, Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# 1. EasyAmplicon(TBD)
# 2. Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录


# 1.1 程序功能描述和主要步骤

# 程序功能：根据sintax注释结果过滤16S扩增子，去除线粒体、叶绿体和非细菌序列
# Functions: Remove chloroplast, mitocondria, and non-bacteria

options(warn = -1) # Turn off warning
## 设置输入输出文件和参数
# 修改下面`default=`后面的文件和参数。
# 输入文件"-i", "--input"，输入文件OTU表
# 物种注释"-t", "--taxonomy"，sintax物种注释文件位置；
# 输出文件"-o", "--output"，输出文件筛选OTU表
# 分组列名"-s", "--stat"，筛选物种注释和统计文件；
# 过滤物种注释"-d", "--discard"，筛选去除的物种。

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
        make_option(c("-i", "--input"), type="character", default="result/raw/otutab.txt",
                    help="OTU table [default %default]"),
        make_option(c("-t", "--taxonomy"), type="character", default="result/raw/otus.sintax",
                    help="sintax taxonomy [default %default]"),
        make_option(c("-s", "--stat"), type="character", default="result/raw/otutab_nonBac.stat",
                    help="Filter stat result [default %default]"),
        make_option(c("-d", "--discard"), type="character", default="result/raw/otus.sintax.discard",
                    help="Filter stat result [default %default]"),
        make_option(c("-o", "--output"), type="character", default="result/otutab.txt",
                    help="Filtered OTU table [default %default]")
    )
    opts = parse_args(OptionParser(option_list=option_list))
#    suppressWarnings(dir.create(opts$output))
}


# # 设置输入文件OTU表和物种注释文件位置
# opts$input = "result/raw/otutab.txt"
# opts$taxonomy = "result/raw/otus.sintax"
# 
# # 设置输出文件筛选OTU表、及筛选物种注释和统计文件位置
# opts$output = "result/otutab.txt"
# opts$stat = "result/raw/otutab_nonBac.txt"

# 读取OTU表和sintax物种注释
otutab = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")
sintax = read.table(opts$taxonomy, header=F, row.names=1, sep="\t", comment.char="")

print(paste0("Input feature table is ", opts$input))
print(paste0("Input sintax taxonomy table is ", opts$taxonomy))

# OTU表汇总
total_reads = colSums(otutab)

# 筛选非细菌、古菌序列
# grep返回位置，grepl返回逻辑值
idx = grepl("Bacteria|Archaea", sintax$V4, perl = T)
nonspecific = sintax[!idx,]
nonspecific_reads = colSums(otutab[rownames(nonspecific),])
sintax = sintax[idx,]

# 筛选叶绿体
idx = grepl("Chloroplast", sintax$V2, perl = T)
chloroplast = sintax[idx,]
chloroplast_reads = colSums(otutab[rownames(chloroplast),])
sintax = sintax[!idx,]

# 筛选线粒体
idx = grepl("Mitochondria", sintax$V2, perl = T)
mitochondria = sintax[idx,]
mitochondria_reads = colSums(otutab[rownames(mitochondria),])
sintax = sintax[!idx,]

# 筛选OTU并按丰度由大到小排序
# 如果直接使用sintax的行名，如果与outtab数量不对等，如过多将出现NA行错误，改为交叉筛选
# otutab = otutab[rownames(sintax),]
idx = rownames(otutab) %in% rownames(sintax)
otutab = otutab[idx,]
idx = order(rowSums(otutab), decreasing = T)
otutab = otutab[idx,]
filtered_reads = colSums(otutab)

# 输出异常OTU编号和总量统计
write.table(rbind(nonspecific, chloroplast, mitochondria), file=opts$discard, append = F, sep="\t", quote=F, row.names=T, col.names=F, eol = "\n")
df = as.data.frame(cbind(total_reads, nonspecific_reads, chloroplast_reads, mitochondria_reads, filtered_reads))
suppressWarnings(write.table(paste("SampleID\t",  sep=""), file=opts$stat, append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = ""))
suppressWarnings(write.table(df, file=paste(opts$stat, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

# 筛选排序后的OTU表
write.table(paste("#OTUID\t",  sep=""), file=paste(opts$output, sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = "")
suppressWarnings(write.table(otutab, file=paste(opts$output, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

# 输出样本基本信息
print(paste0("Summary of samples size in final feature table: "))
summary(filtered_reads)

print(paste0("Onput feature table is ", opts$output))
print(paste0("Detail and statistics in ", opts$stat))
