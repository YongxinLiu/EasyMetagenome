#!/usr/bin/env Rscript
options(warn = -1) # Turn off warning
# 1.2 解析命令行
# 设置清华源加速下载
#library(optparse)
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, 
                                               quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}

# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="taxonom2.spf",
                help="Metadata of metaphlan4 [default %default]"),
    make_option(c("-n", "--number"), type="numeric", default=492,
                help="Top number threshold [default %default]"),
    make_option(c("-c", "--colname"), type="character", default="Group",
                help="Select group column name to use in metadata.txt, e.g. Group [default %default]"),
    make_option(c("-d", "--design"), type="character", default="metadata.tsv",
                help="Design file or metadata [default %default]"),
    make_option(c("-g", "--group"), type="character", default="all",
                help="Specify the desired group, comma seperate. e.g. A,B,C,D [default %default]"),
    make_option(c("-t", "--type"), type="character", default="heatmap",
                help="specify plot type, e.g. heatmap, bar, [default %default]"),
    make_option(c("-o", "--output"), type="character", default=".",
                help="Output pdf directory. [default %default]"))
  
  opts = parse_args(OptionParser(option_list=option_list))
}
if(opts$output==""){opts$output="."}
suppressWarnings(dir.create(opts$output, showWarnings = F,recursive=T))
# 按丰度筛选，如0.01即代表0.01%，即万分之一
#abundance = opts$abundance
# 按数量筛选，如150即代表最高丰度的150个特征
number = opts$number
if (substr(opts$output,length(opts$output),length(opts$output)) == "/||\\."){
  out <- output}else{out <- paste0(opts$output,"/")}

# 2. 读取输入文件
# 读取OTU表
#all_tax <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
all_tax <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
#all_tax <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy = read.table(opts$input,header = T,row.names = NULL,
                      sep = "\t",quote = "",check.names = F)
#taxonomy = read.table("taxonomy2.spf",header = T,row.names = NULL,
#           sep = "\t",quote = "",check.names = F)
tax_id <- intersect(colnames(taxonomy),all_tax)
#tax_id[length(tax_id)]
# taxonomy
taxonomy <- subset(taxonomy,taxonomy[,tax_id[length(tax_id)]] != "unclassified")
# taxonomy[,tax_id[length(tax_id)]]
for(i in tax_id){
  for (j in c(1:nrow(taxonomy))){
    if (length(grep("__",taxonomy[j,i]))!=0){
      taxonomy[j,i] <- substr(taxonomy[j,i],4,nchar(taxonomy[j,i])) 
    }
    else{taxonomy[j,i] <- taxonomy[j,i]
    }
  }
}
row.names(taxonomy) <- taxonomy[,tax_id[length(tax_id)]]
# 读取实验设计
metadata = read.table(opts$design, sep="\t", header = TRUE, row.names = 1, 
                      stringsAsFactors = F, comment.char = "")
#metadata = read.table("metadata.txt", sep="\t", header = TRUE, row.names = 1, 
#           stringsAsFactors = F, comment.char = "")
#colnames(metadata) <- c("Group")
GROUPNAME <- opts$colname
#GROUPNAME
#GROUPNAME  <- c("Group")
metadata <- subset(metadata,select=GROUPNAME)
group <- opts$group
#group <- "Z,DSS"
if (group == "all"){
  metadata <- metadata
  group <- unique(metadata[,GROUPNAME])
  #group
}else{
  group <- unlist(strsplit(group,","))
  metadata <- subset(metadata,metadata[,GROUPNAME] == group)
}

#head(metadata)
# 3. 过滤
# 标准化并求均值
norm = as.data.frame(t(t(taxonomy[,row.names(metadata)])/colSums(taxonomy[,row.names(metadata)],na=T)*100))
# 丰度由大到小排序
idx = order(rowMeans(norm), decreasing = T)
norm = norm[idx,]
######################################################################################
# 按丰度筛选
#idx_abundance = rowMeans(norm) > abundance
#filtered_taxonomy_abundance = norm[idx_abundance,]
#taxonomy_id_abundance <- taxonomy[idx_abundance,][,c(1:8)]
# 按数量筛选
#idx_number = order(rowMeans(norm), decreasing = T)
#filtered_taxonomy_number = head(norm, number)
#taxonomy_id_number <- taxonomy[row.names(filtered_taxonomy_number),]
######################################################################################
filtered_taxonomy = head(norm, number)
#head(filtered_taxonomy)
#
tax_mean <- data.frame(matrix(ncol=0,nrow=nrow(filtered_taxonomy)))
col_names <- data.frame(matrix(ncol=0,nrow=0))
#head(group)
for (i in group){
  temp <- subset(metadata,metadata[,GROUPNAME] == eval(i))
  temp_tax <- filtered_taxonomy[,row.names(temp)]
  # 标准化并求均值
  #tem_norm = as.data.frame(t(t(temp_tax)/colSums(temp_tax,na=T)*100))
  tem_norm_mean <- round(rowMeans(temp_tax),4)
  min_tem_norm_mean <- min(tem_norm_mean)
  max_tem_norm_mean <- max(tem_norm_mean)
  col_names[i,1] <- eval(i)
  col_names[i,2] <- min_tem_norm_mean
  col_names[i,3] <- max_tem_norm_mean
  tax_mean <- cbind(tem_norm_mean,tax_mean)
}
#col_names
colnames(tax_mean) <- col_names[,1]
colnames(col_names) <- c("Inner2Outer","Range_min","Range_max")
# 过滤
# 保存输出文件
# 过滤的OTU表
spf_file <- paste0(out,"taxonomy_mean.spf")
write.table(paste0(tax_id[length(tax_id)],"\t"), file=spf_file, append = F, sep="\t", quote=F, eol = "", row.names=F, col.names=F)
suppressWarnings(write.table(tax_mean, file=spf_file, append = T, sep="\t", quote=F, row.names=T, col.names=T))

#head(tax_mean)
taxonomy_id <- taxonomy[row.names(tax_mean),]
#head(taxonomy_id)
# 读取筛选后的文件，不设置行名
taxid = taxonomy_id
#head(taxonomy_id)
# 筛选门-属5级+OTUID
#all_tax <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
all_tax <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
#all_tax <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
#all_tax
only_tax_id <- intersect(colnames(taxonomy_id),all_tax)
#only_tax_id
tree = data.frame(taxid[,only_tax_id], stringsAsFactors = F)
#ncol(tree)
#head(tree)
## clarify taxonomy，解决不同级别重名问题，为可识别级别，且与Greengene格式保持一致
tree[,1] = paste("p__",tree[,1],sep = "")
tree[,2] = paste("c__",tree[,2],sep = "")
tree[,3] = paste("o__",tree[,3],sep = "")
# tree[,4] = paste("f__",tree[,4],sep = "")
#head(tree)
if (ncol(tree) == 5){
  head(tree)
  tree[,4] = paste("",tree[,4],sep = "")
}else{if (ncol(tree) == 6){
  tree[,5] = paste("g__",tree[,5],sep = "")
}else{
  tree[,5] = paste("g__",tree[,5],sep = "")
  if (ncol(tree) == 7){
    tree[,6] = paste("s__",tree[,6],sep = "")
  }
}
}
#head(tree)
#tree[,7] = paste("t__",tree[,7],sep = "")
# save tree backbone, 按点分隔格式

# 解决科标签重名问题
idx = tree[,4] %in% "Unassigned"
# 方法1. 重名标签添加数字编号，但结果有太多Unassigned
# tree[idx,4] = paste0(tree[idx,4], 1:length(tree[idx,4]))
# 方法2. 过滤掉科末注释的条目，数量会减少，但图片更美观
tree = tree[!idx,]
tax_mean = tax_mean[!idx,]
#head(tree)
# 简化一些代_的不规则科名
tree[,4] = gsub('_\\w*',"",tree[,4])
tree_file <- paste0(out,"tree1_backbone.txt")
write.table(tree, file=tree_file, sep=".", col.names=F, row.names=F, quote=F)
#head(tree)
# 列出现在有门、纲、目、科、属，用于设置与门对应的背景色
Phylum = unique(tree[,1]) 
Class = unique(tree[,2])
Order = unique(tree[,3])
Family = unique(tree[,4])
Genus = unique(tree[,5])
if (ncol(tree) >5){
  Species = unique(tree[,6])
  if (ncol(tree) >6){
    Strain = unique(tree[,7])
  }
}
# 筛选四大菌门中的科并按门着色
# 修改为目，则将tree的4列改为3列，Family改为Order
#tree
pro = tree[tree[,1]=="p__Proteobacteria",4]
act = tree[tree[,1]=="p__Actinobacteria",4] 
bac = tree[tree[,1]=="p__Bacteroidetes",4]
fir = tree[tree[,1]=="p__Firmicutes",4]
fus = tree[tree[,1]=="p__Fusobacteria",4]
#fir
# 对每个科进行标签、文字旋转、按门注释背景色
# 也可调整为其它级别，如Order, Class或Genus
label_color = data.frame(stringsAsFactors = F)

for (element in Phylum)
for (element in Family)
{
  #element
  anno = data.frame(stringsAsFactors = F)
  anno[1,1] = element
  anno[1,2] = "annotation"
  anno[1,3] = "*"
  # 设置文字旋转90度
  anno[2,1] = element
  anno[2,2] = "annotation_rotation"
  anno[2,3] = "90"
  # 设置背景色，四大门各指定一种色，其它为灰色
  anno[3,1] = element
  anno[3,2] = "annotation_background_color"

  if (element %in% pro)
  {
    anno[3,3] = "#aed6c3"
  } else if (element %in% act)
  {
    anno[3,3] = "#F58D8D"
  } else if (element %in% fir)
  {
    anno[3,3] = "#edb5cf"
  } else if (element %in% bac)
  {
    anno[3,3] = "#bfafcd"
  } else if (element %in% fus)
  {
    anno[3,3] = "#afb09b"
  } else {
    anno[3,3] = "grey"
  }
  label_color = rbind(label_color,anno)
}
tree_file2 <- paste0(out,"tree2_label_color.txt")
write.table(label_color, tree_file2, sep = "\t", quote = F,col.names = F,row.names = F, na="")
#head(label_color)
global <- data.frame(
  "global tree options"=c("ignore_branch_len","total_plotted_degrees", "start_rotation",
                          "clade_separation","branch_bracket_depth","branch_bracket_width",
                          "branch_thickness","branch_color","branch_color_from_ancestor",
                          "annotation_background_width","annotation_background_alpha",
                          "annotation_background_separation","annotation_background_offset",
                          "annotation_legend_font_size","clade_marker_size","clade_marker_shape",
                          "clade_marker_edge_width","clade_marker_edge_color","clade_marker_font_color",
                          "p__Actinobacteria*","p__Proteobacteria*","p__Bacteroidetes*",
                          "p__Firmicutes*","p__Proteobacteria","p__Actinobacteria",
                          "p__Bacteroidetes","p__Firmicutes"),
  #"param" = c(0,328,-87,0.0,0.8,0.5,0.75,"black",1,0.6,0.2,-0.5,0.02,7,16,"o",0.25,"black","k",
  #            rep("clade_marker_color",4),rep("annotation_background_color",4)),
  "param" = c(0,360,-90,
              0.0,0.8,0.5,
              0.75,"black",1,
              0.6,0.2,
              -0.5,0.02,
              7,25,"o",
              0.05,"black","k",
              rep("clade_marker_color",4),rep("annotation_background_color",4)),
  "color" = c(rep("",19),"#EF5656","#2BB065","#47B3DA","#F7A415","#B0FFC0","#FFB0B0",
              "#B0E9FD","#FCDBA2"))
#global
#write.table("#", "graphlan_annotate.txt",sep="\t",quote = F,col.names = F,row.names = F, na="")
annotation_file <- paste0(out,"graphlan_annotate.txt")
write.table(global, annotation_file,sep="\t",quote = F,col.names = F,row.names = F, na="")
write.table(label_color, annotation_file, sep = "\t",append=T, quote = F,col.names = F,row.names = F, na="")
cmd_raw1 <- paste0("graphlan_annotate.py --annot ",annotation_file," ",tree_file,"  ",out,"graphlan.xml")
system(cmd_raw1)
# 绘图，size决定图片大小，越大字越小
cmd_raw <- paste0("graphlan.py ",out,"graphlan.xml ",out,"graphlan1_tree_raw.pdf --size 5")
system(cmd_raw)
#heatmap
#col_names
type <- opts$type
if (type == "heatmap"){
  heat_cfg <- data.frame("ring"=c("ring_internal_separator_thickness","ring_width","ring_color"))
  heat_cfg_param <- data.frame("color"=c(0.5,1))
  #colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
  #colors <- c('blue', 'green', 'red', 'cyan', 'magenta',  'yellow', 'black', 'white')
  #colors <- c('blue', 'green', 'red', 'cyan', 'magenta',  'yellow', 'black', 'white')
  colors <- c('#5ebcc2', '#829d45', '#d175ae', 'cyan', 'magenta',  'yellow', 'black', 'white')
  #system("rm track.heatmap")
  for (i in 1:length(group)){
    temp_ring <- heat_cfg
    temp_ring$width <- rep(eval(i))
    temp_param <- heat_cfg_param
    temp_param[3,1] <- colors[i]
    col_names[i,4] <- colors[i]
    temp <- cbind(temp_ring,temp_param)
    #write.table("#", file="graphlan_annotate.txt", append = T, quote=F, eol = "", row.names=F, col.names=F)
    write.table(temp, file=annotation_file, append = T, sep = "\t", quote=F, row.names=F, col.names=F)
    temp_data <- data.frame(matrix(ncol=0,nrow=nrow(tax_mean)))
    temp_data$"ring_alpha" <- rep("ring_alpha")
    temp_data$"ring_number" <- rep(i)
    temp_data$"alpha" <- tax_mean[,group[i]]
    row.names(temp_data) <- row.names(tax_mean)
    write.table(temp_data,file=annotation_file, append = T, sep = "\t",quote = F,col.names = F,row.names = T, na="")
    col_names
    legend_file <- paste0(out,"Ring_legend_inner-outer.txt")
    write.table(col_names,legend_file,quote = F,col.names = T,row.names = F)
    cmd_heatmap1 <- paste0("graphlan_annotate.py --annot ",annotation_file, " ",tree_file,"  ",out,"graphlan.xml")
    # 绘图，size决定图片大小，越大字越小
    system(cmd_heatmap1)
    cmd_heatmap <- paste0("graphlan.py ",out,"graphlan.xml ",out,"graphlan1_tree_heatmap.pdf --size 6")
    system(cmd_heatmap)
  }
}else{if  (type == "bar"){
  bar_cfg <- data.frame("ring"=c("ring_internal_separator_thickness","ring_width","ring_color"))
  bar_cfg_param <- data.frame("color"=c(0.5,0.5))
  #colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
  colors <- c('blue', 'green', 'red', 'cyan', 'magenta',  'yellow', 'black', 'white')
  #system("rm track.heatmap")
  for (i in 1:length(group)){
    temp_bar <- bar_cfg
    temp_bar$width <- rep(eval(i))
    temp_param <- bar_cfg_param
    temp_param[3,1] <- colors[i]
    col_names[i,4] <- colors[i]
    temp <- cbind(temp_bar,temp_param)
    #write.table("#", file="graphlan_annotate.txt", append = T, quote=F, eol = "", row.names=F, col.names=F)
    write.table(temp, file=annotation_file, append = T, sep = "\t", quote=F, row.names=F, col.names=F)
    temp_data <- data.frame(matrix(ncol=0,nrow=nrow(tax_mean)))
    temp_data$"ring_height" <- rep("ring_height")
    temp_data$"ring_number" <- rep(i)
    temp_data$"height" <- tax_mean[,group[i]]
    row.names(temp_data) <- row.names(tax_mean)
    write.table(temp_data,file=annotation_file, append = T, sep = "\t",quote = F,col.names = F,row.names = T, na="")
    col_names
    legend_file <- paste0(out,"Bar_legend_inner-outer.txt")
    write.table(col_names,legend_file,quote = F,col.names = T,row.names = F)
    cmd_heatmap1 <- paste0("graphlan_annotate.py --annot ",annotation_file, " ",tree_file,"  ",out,"graphlan.xml")
    # 绘图，size决定图片大小，越大字越小
    system(cmd_heatmap1)
    cmd_heatmap <- paste0("graphlan.py ",out,"graphlan.xml ",out,"graphlan1_tree_bar.pdf --size 5")
    system(cmd_heatmap)
  }
}
}


