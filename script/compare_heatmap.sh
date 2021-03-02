#!/bin/bash

# 模式设置 遇到错误终止流程 Stop when error
set -e

# 默认参数 Default parameter
input=result/compare/A-B.txt
design=result/metadata.tsv
taxonomy=result/taxonomy.txt
group=Group
output=""
order=FALSE
width=4
height=2.5
# 读取差异OTUs文件的非数据行数，一般为7
non_data_line=7
top_tax=10
cluster_cols=TRUE
text_size=7

# 脚本功能描述 Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
  -------------------------------------------------------------------------------
  Filename:    compare_heatmap.sh
Version:     1.7
Date:        2019/12/23
Author:      Yong-Xin Liu
Email:       metagenome@126.com
Website:     https://blog.csdn.net/woodcorpse
Description: Draw heatmap plot by compare result, must have logFC, logCPM and level
Notes:
  -------------------------------------------------------------------------------
  Copyright:   2016-2020 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
  Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).
-------------------------------------------------------------------------------
  Version 1.0 2018/4/9
Draw heatmap plot by compare result, must have logFC, logCPM and level
# All input and output should be in default directory, or give relative or absolute path by -i/-d
  Version 1.7 2019/12/23
Modified to plotting in EasyAmplicon 1.7 in windows

# Input files: result/metadata, result/compare/A-B.txt

# 1. 差异比较OTU
ID	log2FC	log2CPM	PValue	FDR	level	MeanA	MeanB	KO1	KO2	KO3	KO4	KO5	KO6	WT1	WT2	WT3	WT4	WT5	WT6
ASV_1	-0.622	16.009	0.025974025974026	0.116512059369202	Depleted	3.35	5.155	3.34	5.432	2.094	3.585	2.726	2.926	6.237	6.593	4.691	5.414	3.82	4.176
ASV_2	0.758	16.176	0.0151515151515152	0.0951515151515152	Enriched	5.558	3.286	5.777	3.392	6.043	5.816	7.405	4.914	3.418	3.883	2.291	3.017	3.966	3.14

# Output file
1. heatmap plot in pdf and png

OPTIONS:
-d design for each samples, default result/metadata.tsv
-h figure height, default 8
-i OTU table in reads counts, default result/otutab.txt
-l non_data_line number, v2 is 14, v1 is 12
-m statistics method, default edgeR, alternative wilcon
-o output director, default result/tax/
  -p pvaule, default 0.01
-q FDR/qvalue, default 0.05
-s text size, default 7
-w figure width, default 8
-A group name
-B group selected list, empty will not select
-F fold change, default 1.3
-O order of legend, default FALSE alphabet, set TRUE abundance
-? show help of script

Example:
  compare_heatmap.sh -i ${input} -o ${output} -w ${width} -h ${height}

EOF
}


# 参数解析 Analysis parameter
while getopts "d:f:h:i:l:m:n:o:p:q:s:t:w:A:B:C:F:O:" OPTION
do
	case $OPTION in
		d)
			design=$OPTARG
			;;
		h)
			height=$OPTARG
			;;
		i)
			input=$OPTARG
			;;
		l)
			non_data_line=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		s)
			text_size=$OPTARG
			;;
		t)
			taxonomy=$OPTARG
			;;
		w)
			width=$OPTARG
			;;
		A)
			g1=$OPTARG
			;;
		C)
			cluster_cols=$OPTARG
			;;
		O)
			order=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

# 建立脚本目录
mkdir -p script

# 开始写R统计绘图脚本
cat <<END >script/compare_heatmap.R
#!/usr/bin/env Rscript
# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



#----1.1 程序功能描述和主要步骤#----

# 程序功能：差异比较热图展示
# Functions: pheatmap for different compare
# Main steps:
# - Reads data matrix and design
# - Prepare annotation for rows and columns
# - Save heatmap in input folder

# 清空工作环境 Clean enviroment object
rm(list=ls())


#----1.2 安装CRAN来源常用包#----
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("ggplot2","pheatmap","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}



#----2. 读取输入文件#----

# 读取比较列表
input = read.table("${input}", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
# 筛选显著差异部分
input = subset(input, level %in% c("Enriched","Depleted"))

design = read.table("${design}", header=T, row.names=1, sep="\t", comment.char="")
# 统一改实验列为group
design\$group = design\$${g1}

norm = input[,-(1:${non_data_line})]

if (dim(norm)[1]>1){

  idx = rownames(design) %in% colnames(norm)
  design = design[idx,]

  anno_row = data.frame(Level = input\$level, row.names = rownames(input))
  anno_col = data.frame(Group = design\$group, row.names = rownames(design))


  ## 注释文件存在时，添加物种注释，不聚类分组
  if (file.exists("result/taxonomy.txt")){
    taxonomy = read.table("result/taxonomy.txt", sep = "\t", row.names=1, header=T)
    # per = read.delim("result/taxonomy.txt", sep = "\t", row.names=1, header=T)
    # mean = rowMeans(per)
    # per = as.data.frame(mean[order(mean, decreasing = T)])
    # top_tax=head(rownames(per), n=${top_tax})
    #
    # x = input
    #
    # # 将低丰度的门变为Low Abundance
    # x\$tax = x\$Phylum# factor(x\$Phylum, levels=c(as.vector(unique(x\$Phylum)),"Low Abundance"))
    # # 将门中 proteobacteria 替换为纲
    # x[x\$tax %in% "Proteobacteria",]\$tax =  x[x\$tax %in% "Proteobacteria",]\$Class # no level can get value
    #
    # # 判断是否有需要替换为低丰度的类，没有报错的解决
    # if (length(x[!(x\$tax %in% top_tax),]\$tax > 0)){
    #   x[!(x\$tax %in% top_tax),]\$tax = "Low Abundance" # no level can get value
    # }
    #
    # # 颜色还是不能保证一致，因为不同组门数量不同？？
    # x\$tax = factor(x\$tax, levels=sort(c(top_tax,"Low Abundance")))
    input\$Phylum = taxonomy[rownames(input),"Phylum"]
    anno_row = data.frame(Level = input\$level, Taxonomy = input\$Phylum, row.names = rownames(input))
  }

  pheatmap(norm,
           scale = "row",
           cutree_rows=2,cutree_cols = 2,
           cluster_cols = ${cluster_cols},
           annotation_col = anno_col,
           annotation_row = anno_row,
           filename = paste("${output}", ".heatmap.pdf", sep=""),
           width=$width, height=$height,
           annotation_names_row= T,annotation_names_col=T,
           show_rownames=T,show_colnames=T,
           fontsize=$text_size,display_numbers=F)

  # pheatmap(norm,
  #          scale = "row",
  #          cutree_rows=2,cutree_cols = 2,
  #          cluster_cols = ${cluster_cols},
  #          annotation_col = anno_col,
  #          annotation_row = anno_row,
  #          filename = paste("${output}", ".heatmap.png", sep=""),
  #          width=$width, height=$height,
  #          annotation_names_row= T,annotation_names_col=T,
  #          show_rownames=T,show_colnames=T,
  #          fontsize=$text_size,display_numbers=F)
  # # 提示工作完成
  # print(paste("Output in ${output}", "_heatmap.pdf finished.", sep = ""))
}
END

# 执行脚本，脚本运行目录即工作目录(与脚本位置无关)
# output=dirname ${output}
# mkdir -p ${output}
Rscript script/compare_heatmap.R
