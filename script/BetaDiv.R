
#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：beta多样性排序 三种方法群落差异检测，群落样本分组两组之间检测
# Functions: beta diversity plot

## 参数解释
## 设置输入输出文件和参数
# 修改下面`default=`后面的文件和参数。
#
# 输入文件为原始otu表格(otutab.txt)+分组信息(metadata.tsv)
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
#
# 排序方法"-m", "--method"；这里可选九中排序方法：DCA, CCA, RDA, NMDS, MDS, PCoA,PCA,LDA,t-sne
#
# 群落显著性检测方法"-p", "--permutation"：可选adonis，anosim，MRPP
#
#距离算法"-x", "--distance"；这里使用phyloseq封装的多种距离：输入距离的编号，这些距离按照顺序排列为：
# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski"  "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial"  "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# "binary" "minkowski""ANY" ；输入对于距离的顺序数据即可使用
#
# 进化树"-t", "--tree"：这里用于指定需要进化树的距离。


options(warn = -1) # Turn off warning


site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="../data/otutab.txt",
                help="otu table [default %default]"),
    make_option(c("-d", "--design"), type="character", default="../data/metadata.tsv",
                help="Design file or metadata [default %default]"),
    make_option(c("-t", "--tree"), type="character", default="",
                help="Community difference detection method [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=170,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=110,
                help="Figure heidth in mm [default %default]"),
    make_option(c("-m", "--method"), type="character", default="DCA",
                help="method could chioce: DCA, CCA, RDA, NMDS, MDS, PCoA,PCA,LDA,t-sne [default %default]"),
    make_option(c("-p", "--permutation"), type="character", default="adonis",
                help="Community difference detection method:adonis，anosim，MRPP [default %default]"),
    make_option(c("-x", "--distance"), type="numeric", default=8,
                help="Community difference detection method [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]")


  )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}


# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf

if(opts$output==""){opts$output="./"}


#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
map = read.delim(opts$design,row.names = 1)
#----2.2 otu table otu#----
otu = read.delim(opts$input,row.names = 1)

# tree = read_tree("./otus.tree")

source("../BetaDiv-wt.R")
source("../pairMicroTest-wt.R")
source("../MicroTest-wt.R")

  otu = phyloseq::otu_table(otu,taxa_are_rows=TRUE)
  map = phyloseq::sample_data(map)
#   ps = phyloseq::merge_phyloseq(otu,map)
#
#   ass = as.character(opts$distance)
#
#
#   dist_methods <- unlist(phyloseq::distanceMethodList)
#
#
# unif <- phyloseq::distance(ps, method=dist_methods[opts$chosemethod])

if (opts$tree=="") {
  tree = NULL
} else{
  tree = phyloseq::read_tree(opts$tree)
}
# tree = read_tree("./otus.tree")


result = BetaDiv(otu, map, tree = tree,
                 dist = opts$distance, group = opts$group, method =opts$method,
                 pvalue.cutoff = 0.05, Micromet = opts$permutation)


#提取图片
p1 = result[[1]]

#提取作图数据
data1 = result[[2]]

# 图形做标签
p2 = result[[3]]

# 提取两两比对结果
data2 = result[[4]]

# 提取整体比较结果
data3 = result[[5]]
dist_methods <- unlist(phyloseq::distanceMethodList)
dist =   dist_methods[opts$distance]

cs = paste(dist,opts$method,opts$permutation,sep = "_")
o=paste0(opts$output,paste("/",cs,"plot.pdf",sep = ""))
# o=paste0(opts$output,"/plot.pdf")
ggsave(o, p1, width = opts$width, height = opts$height, units = "mm")
o=paste0(opts$output,paste("/",cs,"Labelplot.pdf",sep = ""))
ggsave(o, p2, width = opts$width, height = opts$height, units = "mm")


o=paste0(opts$output,"/betaplotData.csv")
write.csv(data1,o)

o=paste0(opts$output,paste("/",cs,"betaplotData.csv",sep = ""))
write.csv(data2,o)

o=paste0(opts$output,paste("/",cs,"total.csv",sep = ""))
write.csv(data3,o)
