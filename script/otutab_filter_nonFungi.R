#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 设置输入文件OTU表和物种注释文件位置
otutab_dir = "result/raw/otutab.txt"
sintax_dir = "result/raw/otus.sintax"

# 设置输出文件筛选OTU表、及筛选物种注释和统计文件位置
output_dir = "result/otutab.txt"
stat_dir = "result/raw/otutab_nonBac.txt"

# 读取OTU表和sintax物种注释
otutab = read.table(otutab_dir, header=T, row.names=1, sep="\t", comment.char="")
sintax = read.table(sintax_dir, header=F, row.names=1, sep="\t", comment.char="")

print(paste0("Input feature table is ", otutab_dir))
print(paste0("Input sintax taxonomy table is ", sintax_dir))

# OTU表汇总
total_reads = colSums(otutab)

# 筛选非细菌、古菌序列
# grep返回位置，grepl返回逻辑值
idx = grepl("Fungi", sintax$V4, perl = T)
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
write.table(rbind(nonspecific, chloroplast, mitochondria), file=stat_dir, append = F, sep="\t", quote=F, row.names=T, col.names=F, eol = "\n")
df = as.data.frame(cbind(total_reads, nonspecific_reads, chloroplast_reads, mitochondria_reads, filtered_reads))
suppressWarnings(write.table(paste("\nSampleID\t",  sep=""), file=stat_dir, append = T, sep="\t", quote=F, row.names=F, col.names=F, eol = ""))
suppressWarnings(write.table(df, file=paste(stat_dir, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

# 筛选排序后的OTU表
write.table(paste("#OTUID\t",  sep=""), file=paste(output_dir, sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = "")
suppressWarnings(write.table(otutab, file=paste(output_dir, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

# 输出样本基本信息
print(paste0("Summary of samples size in final feature table: "))
summary(filtered_reads)

print(paste0("Onput feature table is ", output_dir))
print(paste0("Detail and statistics in ", stat_dir))
