#!/bin/bash

#set -x

usage()
{
    cat <<EOF
    ${txtcyn}

    ***CREATED BY Chen Tong (chentong_biology@163.com)***

        Usage:

        $0 options${txtrst}

    ${bldblu}Function${txtrst}:

        This script is used to draw a line or multiple lines using ggplot2.
    You can specify whether or not smooth your line or lines.

    fileformat for -f (suitable for data extracted from one sample, the
                       number of columns is unlimited. Column 'Set' is not necessary)
    ------------------------------------------------------------
        Gene Sample
    g1	h3k27ac
    g2	h3k27ac
    a1	h3k27ac
    a3	h3k27ac
    b4	h3k27ac
    g1	ctcf
    h1	ctcf
    a3	ctcf
    b1	ctcf
    b2	ctcf
    g2	ctcf
    -------------------------------------------------------------

        ${txtbld}OPTIONS${txtrst}:
        -f	Data file (without header line, the first column is the
                      name of genes or otehr things you want to compare, the second
                      column is sample name, tab seperated)
    ${bldred}[NECESSARY]${txtrst}
    -F	If you have the length of each set and the number
    overlaps between or among each set,  pleass give TRUE here.
    ${bldred}[When this is TRUE, the value given to -f would be
              the prefix for the output figure.]${txtrst}
    -t	Title for the picture.
    -a	The name for label1.
    ${bldred}[Necessary when -f is not FALSE,
              one string in your second column,ordered. ]${txtrst}
    -b	The name for label2.
    ${bldred}[Necessary when -f is not FALSE,
              one string in your second column,
              ordered. A parameter to -b is needed for 2-way venn.]${txtrst}
    -c	The name for label3.
    ${bldred}[Necessary when -f is not FALSE,
              one string in your second column,
              ordered. A parameter to -c is needed for 3-way venn.]${txtrst}
    -d	The name for label4.
    ${bldred}[Necessary when -f is not FALSE,
              one string in your second column],
    ordered. A parameter to -d is needed for 4-way venn.]${txtrst}
-g	The name for label5.
${bldred}[Necessary when -f is not FALSE,
          one string in your second column],
ordered. A parameter to -f is needed for 5-way venn.]${txtrst}
-p	The label of output files.
${bldred}[Optional. But when one main file is used to generate
          multiple vennDIagrams, this parameter should be specified.]
${txtrst}
-n	List of numbers for venn plot (used when -F is TRUE).

For two-set venn, the format is "100, 110, 50" represents
(length_a, length_b,
    a_b_overlap).

For three-set venn, the format is "100, 110, 90, 50, 40, 40, 20"
represents (length_a, length_b, length_c,
            a_b_overlap,  b_c_overlap, a_c_overlap, a_b_c_overlap).

For four-set venn, the format is "100, 110, 90, 50, 40, 40, 20"
represents (length_a, length_b, length_c,
            a_b_overlap, a_c_overlap, a_d_overlap, b_c_overlap,
            b_d_overlap, c_d_overlap, abc_overlap, abd_overlap,
            acd_overlap, bcd_overlap, abcd_overlap).

-l	List of labels for venn plot (used when -F is TRUE).

Format: "'a', 'b'" for two-set and "'a', 'b', 'c'" for three-set.

Pay attention to the order.
Both double-quotes and single-quotes are needed.

-C	Color for each area.
[${txtred}Ususlly the number of colors should
    be equal to the number of labels.
    If you manually set colors for 4-way
    venn diagram, the first color will be given to the
    leftmost set, the second will be given to the rightmost
    set, the third will be given to second leftmost and the forth
    will be given to the second rightmost. You may want to change
    the color yourself.
    "'red','green','pink','blue','cyan','green','yellow'" or
    "rgb(255/255,0/255,0/255),rgb(255/255,0/255,255/255),rgb(0/255,0/255,255/255),
		rgb(0/255,255/255,255/255),rgb(0/255,255/255,0/255),rgb(255/255,255/255,0/255)"
    ${txtrst}]
-w	The width of output picture.[${txtred}Default 10${txtrst}]
-u	The height of output picture.[${txtred}Default 10${txtrst}]
-E	The type of output figures.[${txtred}Default pdf, accept
                               eps/ps, tex (pictex), png, jpeg, tiff, bmp, svg and
                               wmf [Only png, eps, png(recommend if you want eps) is available, others are unuseable]${txtrst}]
-r	The resolution of output picture.[${txtred}Default 300 ppi${txtrst}]
-e	Execute or not[${bldred}Default TRUE${txtrst}]
-i	Install depended packages[${bldred}Default FALSE${txtrst}]

Example:
    s-plot vennDiagram -f prefix -F TRUE -n "120, 110, 50"  -l "'a','b'"
s-plot vennDiagram -f file -a h3k27ac -b ctcf
s-plot vennDiagram -f file -a h3k27ac -b ctcf -p k27_ctcf

EOF
}

file=
    numGiven="FALSE"
title=""
label1="CHENTONG"
label2="CHENTONG"
label3="CHENTONG"
label4="CHENTONG"
label5="CHENTONG"
numList=
    labelList=
    execute='TRUE'
ist='FALSE'
uwid=10
vhig=10
res=300
ext='pdf'
line_size=1
#color_v='"cornflowerblue", "green", "yellow", "darkorchid1"'
color_v='"dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"'
prefix=''

while getopts "hf:F:t:a:b:c:d:g:n:l:C:p:w:u:r:e:E:i:" OPTION
do
case $OPTION in
h)
usage
exit 1
;;
f)
file=$OPTARG
;;
F)
numGiven=$OPTARG
;;
t)
title=$OPTARG
;;
a)
label1=$OPTARG
;;
b)
label2=$OPTARG
;;
c)
label3=$OPTARG
;;
d)
label4=$OPTARG
;;
g)
label5=$OPTARG
;;
n)
numList=$OPTARG
;;
l)
labelList=$OPTARG
;;
p)
prefix=$OPTARG
;;
C)
color_v=$OPTARG
;;
w)
uwid=$OPTARG
;;
u)
vhig=$OPTARG
;;
r)
res=$OPTARG
;;
E)
ext=$OPTARG
;;
e)
execute=$OPTARG
;;
i)
ist=$OPTARG
;;
?)
usage
exit 1
;;
esac
done

if [ -z $file ]; then
usage
exit 1
fi

if [ -z ${prefix} ]; then
mid='.vennDiagram'
else
    mid='.'${prefix}'.vennDiagram'
fi

cat <<END >${file}${mid}.r

if ($ist){
    install.packages("VennDiagram", repo="http://cran.us.r-project.org")
}
library(VennDiagram)


if ("${ext}" == "png") {
    png(filename="${file}${mid}.png", width=$uwid, height=$vhig,
        res=$res, units="cm")
} else if ("${ext}" == "eps") {
    postscript(file="${file}${mid}.eps", onefile=FALSE, horizontal=FALSE,
               paper="special", width=10, height = 12, pointsize=10)
} else if ("${ext}" == "pdf") {
    pdf(file="${file}${mid}.pdf", width=$uwid, height=$vhig)
}else if ("${ext}" == "svg") {
    svg(filename="${file}${mid}.svg", width=$uwid, height=$vhig,
        pointsize=10)
} else {
    print("This format is currently unsupported. Please check the file <Rplots.pdf> in current directory.")
}


if (! ${numGiven}) {

    data <- read.table(file="$file", sep="\t", quote="")

    num <- 0

    if("${label1}" != "CHENTONG"){
        $label1 <- data[grepl("\\\\<${label1}\\\\>",data[,2]),1]
        num <- num + 1
    }

    if("${label2}" != "CHENTONG"){
        $label2 <- data[grepl("\\\\<${label2}\\\\>",data[,2]),1]
        num <- num + 1
    }

    if("${label3}" != "CHENTONG"){
        $label3 <- data[grepl("\\\\<${label3}\\\\>",data[,2]),1]
        num <- num + 1
    }

    if("${label4}" != "CHENTONG"){
        $label4 <- data[grepl("\\\\<${label4}\\\\>",data[,2]),1]
        num <- num + 1
    }

    if("${label5}" != "CHENTONG"){
        $label5 <- data[grepl("\\\\<${label5}\\\\>",data[,2]),1]
        num <- num + 1
    }

    color_v <- c(${color_v})[1:num]
    #label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),

    if(num == 5){
        p <- venn.diagram(
            x = list($label1=$label1, $label4=$label4,
                     $label5=$label5, $label2=$label2,
                     $label3=$label3),
            filename = NULL, col = "black", lwd = 1,
            fill = color_v,
            alpha = 0.50, main="${title}",
            label.col = c("black"),
            cex = 1, fontfamily = "Helvetica",
            cat.col = color_v,cat.cex = 0.6, margin=0.3,
            cat.fontfamily = "Helvetica"
        )
    }else if(num == 4){
        p <- venn.diagram(
            x = list($label1=$label1, $label4=$label4, $label2=$label2,
                     $label3=$label3),
            filename = NULL, col = "black", lwd = 1,
            fill = color_v,
            alpha = 0.50, main="${title}",
            label.col = c("black"),
            cex = 1, fontfamily = "Helvetica",
            cat.col = color_v, cat.cex = 0.8, margin=0.2,
            cat.fontfamily = "Helvetica",
        )
    } else if (num==3) {
        p <- venn.diagram(
            x = list($label1=$label1, $label2=$label2, $label3=$label3),
            filename = NULL, col = "transparent",
            fill = color_v,
            alpha = 0.50, main="${title}",
            label.col = c("black", "black", "black", "black", "black", "black", "black"),
            cex = 1, fontfamily = "Helvetica", cat.default.pos="text",
            cat.pos=0,  magrin=0.1,
            cat.col = color_v,cat.cex = 1,cat.fontfamily = "Helvetica"
        )
    } else if (num==2) {
        p <- venn.diagram(
            x = list($label1=$label1, $label2=$label2),
            filename = NULL, col = "transparent",
            fill = color_v,
            alpha = 0.50, main="${title}",
            label.col = c("black"),
            cex = 1, fontfamily = "Helvetica",
            cat.default.pos="outer",
            cat.pos=0, margin=0.1,
            cat.col = color_v,cat.cex = 1,cat.fontfamily = "Helvetica"
        )
    }
    grid.draw(p)
} else {
    #---venn plot for given numbers---------
    numList <- c(${numList})
    labelList <- c(${labelList})

    num <- length(labelList)
    color_v <- c(${color_v})[1:num]

    if (num==2) {
        draw.pairwise.venn(area1=numList[1], area2=numList[2],
                           cross.area=numList[3], category=labelList, lwd=rep(1,1),
                           lty=rep(2,2), col="transparent", fill=color_v,
                           cat.col=color_v)
    } else if (num==3) {
        draw.triple.venn(area1=numList[1], area2=numList[2],
                         area3=numList[3], n12=numList[4], n23=numList[5],
                         n13=numList[6], n123=numList[7],
                         category=labelList, col="transparent", fill=color_v,
                         cat.col=color_v, reverse=FALSE)
    }else if (num==4) {
        draw.quad.venn(area1=numList[1], area2=numList[2],
                       area3=numList[3], area4=numList[4], n12=numList[5],
                       n13=numList[6], n14=numList[7], n23=numList[8],
                       n24=numList[9], n34=numList[10], n123=numList[11],
                       n124=numList[12], n134=numList[13], n234=numList[14],
                       n1234=numList[15],
                       category=labelList, col="transparent", fill=color_v,
                       cat.col=color_v, reverse=FALSE)
    }else if (num==5){
        #draw.quintuple.venn()
    }
}

dev.off()

END

if [ "$execute" == "TRUE" ]; then
Rscript ${file}${mid}.r
if [ "$?" == "0" ]; then
/bin/rm -f VennDiagram*.log
fi
fi

