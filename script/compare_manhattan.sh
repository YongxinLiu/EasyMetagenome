#!/bin/bash
set -e

# Default parameters
# 输入文件
input=result/compare/KO-WT.txt
taxonomy=result/taxonomy.txt
# 门水平汇总表，用于统一多个差异比较的图例
phylum=result/tax/sum_p.txt
# 图片宽高mm，字号pt
width=183
height=59
text_size=7
# manhattan plot ylab max size
ymax=20
# 调整图例数量
legend=9
level="Phylum"

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    compare_manhattan.sh
Revision:    1.8
Date:        2020/5/24
Author:      Yong-Xin Liu
Email:       yxliu@genetics.ac.cn
Website:     http://bailab.genetics.ac.cn/
Description: This script is used to perform draw manhattan plot
-------------------------------------------------------------------------------
Copyright:   2016-2020 (c) Yong-Xin Liu
License:     GPL
-------------------------------------------------------------------------------
Version 1.0 2018/5/22
Input compared result as input
Version 1.1 2018/6/14
注释门改为门+变形菌纲着色，X轴添加标签
Version 1.7 2019/12/23
输入文件为差异比较结果和物种注释
Version 1.8 2020/5/24
输出文件可自定义

# 输入文件两个
# 1. OTU差异比较结果
ID	log2FC	log2CPM	PValue	FDR	level	MeanA	MeanB	KO1	KO2	KO3	KO4	KO5	KO6	WT1	WT2	WT3	WT4	WT5	WT6
ASV_1	-0.622	16.009	0.025974025974026	0.116512059369202	Depleted	3.35	5.155	3.34	5.432	2.094	3.585	2.726	2.926	6.237	6.593	4.691	5.414	3.82	4.176
ASV_2	0.758	16.176	0.0151515151515152	0.0951515151515152	Enriched	5.558	3.286	5.777	3.392	6.043	5.816	7.405	4.914	3.418	3.883	2.291	3.017	3.966	3.14

# 2. 七级物种注释
OTUID	Kingdom	Phylum	Class	Order	Family	Genus	Species
ASV_1	Bacteria	Actinobacteria	Actinobacteria	Actinomycetales	Thermomonosporaceae	Unassigned	Unassigned
ASV_2	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Comamonadaceae	Pelomonas	Pelomonas_puraquae

# 3. 门水平汇总表
Phylum	KO1	KO2	KO3	KO4	KO5	KO6	OE1	OE2	OE3	OE4	OE5	OE6	WT1	WT2	WT3	WT4	WT5	WT6	All
Proteobacteria	65.2	49.4	60.7	61.6	72.9	67.9	55.8	47.1	50.8	56	58.4	55.8	50.2	59.8	55.3	55.6	56.6	56.8	48.8
Actinobacteria	25.4	40.4	27.9	28.2	17.2	23.8	28.2	34.7	32.7	27.3	30	28.9	34.8	29.4	31.6	33.7	33.5	30.9	26.5

# Output file
1. Manhattan plot: 必须定义输出文件，推荐SampleAvsSampleB.manhattan.pdf，可选png

OPTIONS:
	-v figure height, default 59 mm
	-n show top N taxonomy number, default 9, recommend phylum is 8-13
	-o output director, default same with input
	-s text size, default 7
	-t default result/taxonomy.txt
	-w figure width, default 183
	-h/? show help of script

Example:
    # i差异比较结果,t物种注释,p图例,w宽,v高,s字号,l图例最大值
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_p.txt \
       -w 183 -v 59 -s 7 -l 10 \
       -o result/compare/${compare}.manhattan.p.pdf
    # 上图只有6个门，切换为纲c和-L Class展示细节
    bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/tax/sum_c.txt \
       -w 130 -v 59 -s 7 -l 10 -L Class \
       -o result/compare/${compare}.manhattan.c.pdf
EOF
}


# Analysis parameter
while getopts "i:l:o:p:s:t:v:w:y:L:" OPTION
do
	case $OPTION in
		i)
			input=$OPTARG
			;;
		l)
			legend=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		p)
			phylum=$OPTARG
			;;
		s)
			text_size=$OPTARG
			;;
		t)
			taxonomy=$OPTARG
			;;
		v)
			height=$OPTARG
			;;
		w)
			width=$OPTARG
			;;
		y)
			ymax=$OPTARG
			;;
		L)
			level=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done


cat <<END >script/compare_manhattan.r
# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("ggplot2","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(), panel.grid=element_blank(),
	axis.line.x=element_line(size=.5, colour="black"), axis.line.y=element_line(size=.5, colour="black"),
	axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=${text_size}),
	legend.position="right", legend.background=element_blank(), legend.key=element_blank(), legend.text= element_text(size=${text_size}),
	text=element_text(family="sans", size=${text_size}))

# 实验差异比较结果
x = read.table("${input}", header=T, row.names= 1, sep="\t", stringsAsFactors = F)
# 只提取前14列
x = x[,1:7]
x = na.omit(x)

# 添加物种注释
taxonomy = read.table("${taxonomy}", sep = "\t", row.names=1, header=T, stringsAsFactors = F)
taxonomy = taxonomy[rownames(x),]
x = cbind(x, taxonomy)

# P值求负对数
x\$neglogp = -log10(x\$PValue)

x\$otu=rownames(x)
x = arrange(x, Kingdom, Phylum, Class, Order, Family, Genus, Species)
x\$otu = factor(x\$otu, levels=x\$otu)   # set x order
x\$num = 1:dim(x)[1]

# 读取高丰度门，用于着色
per= read.delim("${phylum}", sep = "\t", row.names=1, header=T)
mean = rowMeans(per)
per = as.data.frame(mean[order(mean, decreasing = T)])
top_tax=head(rownames(per), n=${legend})

# 将低丰度的门变为Low Abundance
x\$tax = x\$${level}

# 将门中 proteobacteria 替换为纲
# x[x\$tax %in% "Proteobacteria",]\$tax = x[x\$tax %in% "Proteobacteria",]\$Class # no level can get value

if (length(unique(x\$tax)) > length(top_tax)){
	x[!(x\$tax %in% top_tax),]\$tax = "Low Abundance" # no level can get value
}

# 调置标签顺序
label = unique(x\$tax)
label = label[!(label %in% "Low Abundance")] # 删除低丰度
x\$tax = factor(x\$tax, levels = c(label, "Low Abundance"))
# 计算标签中位数
temp = x[x\$tax %in% label, c("tax","num")]
mat_mean = aggregate(temp[,-1], by=temp[1], FUN=mean) # mean


# 调整过大的负对数
if (max(x\$neglogp)>${ymax}){
  x[x\$neglogp>${ymax},]\$neglogp = ${ymax}
}
colnames(x)[2]="log2CPM"
# Manhattan plot
FDR = min(x\$neglogp[x\$level!="NotSig"])
p = ggplot(x, aes(x=num, y=neglogp, color=tax, size=log2CPM, shape=level)) +
  geom_point(alpha=.7) +
  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
  scale_shape_manual(values=c(25, 17, 20))+
  scale_size(breaks=c(5, 10, 15)) +
  labs(x="Phylum of Features", y="-log10(P)") +
  labs(title=paste(gsub(".txt", "", basename("${input}")), sep=" ")) +
  main_theme +
  # theme(legend.position="top") +
  scale_x_continuous(breaks=mat_mean\$x, labels=mat_mean\$tax) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) # + ylim(0,${ymax})
  # ylim保持多组y轴范围一致，一般在最终出版时才使用
p
ggsave(file=paste("${output}", sep=""), p, width = ${width}, height = ${height}, useDingbats=F, units = "mm")
# ggsave(file=paste("${input}_manhattan.png", sep=""), p, width = ${width}, height = ${height})

END

Rscript script/compare_manhattan.r
