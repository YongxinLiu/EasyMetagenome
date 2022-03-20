
if (FALSE){
	install.packages("pheatmap", repo="http://cran.us.r-project.org")
}

library(grid)
library(pheatmap)

if(1){
	library(RColorBrewer)
}

find_coordinates = function(n, gaps, m=1:n) {
	if(length(gaps)==0){
		return(list(coord=unit(m/n, "npc"), size=unit(1/n,"npc")))
	}

	if(max(gaps)>n){
		stop("Gaps do not match matrix size")
	}

	size = (1/n)*(unit(1, "npc")-length(gaps)*unit("4", "bigpts"))

	gaps2 = apply(sapply(gaps, function(gap, x){x>gap}, m), 1, sum)
	coord = m * size + (gaps2 * unit("4", "bigpts"))

	return(list(coord=coord, size=size))
}


vjust <- 0
hjust <- 0.5

if(270==270){
	vjust <- 0.5
	hjust <- 0
}else if(270==45){
	vjust <- .5
	hjust <- 1
}else if(270==0){
	vjust <- 1
	hjust <- 0.5
}



draw_colnames_custom <- function (coln, gaps, ...){
	coord = find_coordinates(length(coln),  gaps)
	x = coord$coord - 0.5 * coord$size
	if (270  == 90){
		hjust = 1
		vjust = 0.5
	}
	if (270  == 45){
		hjust = 1
		vjust = 0.5
	}

	res = textGrob(coln, x=x, y=unit(1, "npc")-unit(3, "bigpts"),
		vjust = vjust, hjust=hjust, rot=270, gp=gpar(...))
	return(res)
}


# Overwrite default draw_colnames with your own version
assignInNamespace(x="draw_colnames", value="draw_colnames_custom",
	ns=asNamespace("pheatmap"))

data <- read.table(file="humann2/pathabundance_relab_unstratified.tsv", sep="\t", header=T, row.names=1,
	check.names=F, quote="", comment="")

if ("FALSE" != "FALSE"){
	#data[data==0] <- 1.0000001
	#data[data==1] <- 1.0001
	data <- FALSE(data+1)
}

if (1 == 1){
	legend_breaks = NA
} else if (1 == 2){
	if (Inf == Inf){
		summary_v <- c(t(data))
		legend_breaks <- unique(c(seq(summary_v[1]*0.95,summary_v[2],length=6),
		  seq(summary_v[2],summary_v[3],length=6),
		  seq(summary_v[3],summary_v[5],length=5),
		  seq(summary_v[5],summary_v[6]*1.05,length=5)))
	} else {
		legend_breaks <- unique(c(seq(summary_v[1]*0.95, Inf,
		 length=10), seq(Inf,summary_v[6]*1.05,length=10)))
	}

	if("FALSE" != "FALSE"){
		legend_breaks <- prettyNum(legend_breaks, digits=FALSE)
	}

	print(col)
	print(legend_breaks)
} else {
	legend_breaks <- c()
}


if ("NA" != "NA") {
	annotation_row <- read.table(file="NA", header=T,
		row.names=1, sep="\t", quote="", check.names=F, comment="")
} else {
	annotation_row <- NA
}

if ("group.txt" != "NA") {
	annotation_col <- read.table(file="group.txt", header=T,
		row.names=1, sep="\t", quote="", check.names=F, comment="")
	# Do not remember what this is for?
	#levs <- unique(unlist(lapply(annotation_col, unique)))
	#annotation_col <- data.frame(lapply(annotation_col, factor,
	#	levels=levs), row.names=rownames(annotation_col))
} else {
	annotation_col <- NA
}

data[data>Inf] <- Inf
if ("-Inf" != "-Inf"){
	data[data<-Inf] <- -Inf
}

if ("function" == "function"){
	color_vector <- colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100)
} else if ("function" == "vector"){
	colfunc <- colorRampPalette(colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100), bias=1)
	color_vector <- colfunc(30)
} else {
	color_vector <- colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100)
}

ann_colors = list(NA)

if (ann_colors[1][1] == "NA") {
	ann_colors = NA
}

pheatmap(data, kmean_k=NA, color=color_vector,
scale="none", border_color=NA,
cluster_rows=TRUE, cluster_cols=TRUE,
breaks=legend_breaks, clustering_method="complete",
clustering_distance_rows="correlation",
clustering_distance_cols="correlation",
legend_breaks=legend_breaks, show_rownames=TRUE, show_colnames=TRUE,
main="", annotation_col=annotation_col,
annotation_row=annotation_row,
annotation_colors = ann_colors,
fontsize=14, filename="humann2/pathabundance_relab_unstratified.tsv.pheatmap.pdf", width=20,
height=50)


