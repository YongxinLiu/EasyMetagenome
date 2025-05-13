[TOC]

# 易宏基因组流程EasyMetagenomePipeline

    # 版本: 1.21, 2024/5/17
    # 测试环境为Windows 10+ / MacOS 10+
    
    # 设置结果目录(通常为项目中的result，此处为12个样本结果result12)
    # 可选癌症cancer或年龄age各18个样本数据示例
    # mac/linux 需要修改路径格式如下， ~ 代表家目录
    # wd=~/meta/result12
    # sd=~/EasyMicrobiome/script
    wd=/d/EasyMetagenome/cancer
    # 设置脚本所在目录(Script Directory)，系统为win/mac/linux
    sd=/d/EasyMicrobiome/script
    PATH=$PATH:$sd/../win:$sd
    # 进入结果目录
    cd $wd

## 物种Metaphlan4

### 热图

    # 显示脚本帮助 help
    Rscript ${sd}/metaphlan_hclust_heatmap.R -h
    # 按指定分类汇总、排序并取Top25种绘制热图
    # -i输入metaphlan4结果转换的spf文件；
    # -t指定分类级别，可选Kingdom/Phylum/Class/Order/Family/Genus/Species/Strain(界门纲目科属种株)，推荐门，目，属
    # -n 输出物种数量，默认为25，最大值为该类型的数量
    # -w、-e指定图片的宽和高，单位为毫米(mm)
    # -o输出图pdf、表txt前缀，默认Heatmap+(-t)+(-n)

    # 科水平Top25热图
    csvtk -t stat metaphlan4/taxonomy.spf
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i metaphlan4/taxonomy.spf \
      -t Family -n 25 \
      -w 183 -e 118 \
      -o metaphlan4/HeatmapFamily

    # 属水平的Top30
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i metaphlan4/taxonomy.spf \
      -t Genus -n 30 \
      -w 183 -e 118 \
      -o metaphlan4/HeatmapGenus

### 维恩图

    ### 筛选每个样本>0.5%的分类单元，包括界门纲目科属种
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
       else {for(i=2;i<=NF;i++) if($i>0.5) print $1, a[i];}}' \
       metaphlan4/taxonomy.tsv > metaphlan4/taxonomy_high.tsv
    wc -l metaphlan4/taxonomy_high.tsv
    # 在线绘制，支持实时查看元素交集 http://www.ehbio.com/test/venn/#/
    # 引文：Mei Yang, Tong Chen, Yong-Xin Liu, Luqi Huang. 2024. Visualizing set relationships: 
    # EVenn's comprehensive approach to Venn diagrams. iMeta 3: e184. https://doi.org/10.1002/imt2.184
    # 本地5组比较:-f输入文件,-a/b/c/d/g分组名,-w/u为宽高英寸,-p输出文件名后缀
    bash ${sd}/sp_vennDiagram.sh -f metaphlan4/taxonomy_high.tsv \
      -a C1 -b C2 -c H1 -d H2 -g E1 \
      -w 4 -u 4 \
      -p C1_C2_H1_H2_E1
    
    # 求均值再两组比较
    sed -i '/^#/d' metaphlan4/taxonomy.tsv
    Rscript ${sd}/otu_mean.R --input metaphlan4/taxonomy.tsv \
      --metadata metadata.txt \
      --group Group --thre 0 \
      --scale F --all TRUE --type mean \
      --output metaphlan4/group_mean.txt    
    ### 筛选每个组>0.5%的分类单元，包括界门纲目科属种
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
       else {for(i=2;i<=NF;i++) if($i>0.5) print $1, a[i];}}' \
       metaphlan4/group_mean.txt > metaphlan4/group_high.tsv
    bash ${sd}/sp_vennDiagram.sh -f metaphlan4/group_high.tsv \
      -a Cancer -b Healthy -c Early \
      -w 4 -u 4 \
      -p Cancer_Healthy_Early      
      
## 功能HUMAnN3

### 分组聚类热图

    bash $sd/sp_pheatmap.sh
    cut -f 1-2 metadata.txt > group.txt
    # -f输入文件，-H水平聚类，u/v图片宽/高，-P添加行注释文件，-Q添加列注释
    bash $sd/sp_pheatmap.sh \
      -f humann3/pathabundance_relab_unstratified.tsv \
      -H 'TRUE' -u 20 -v 50 \
      -Q group.txt
    # 结果为 输入文件名+pheamap.r/pdf，代码和图片

    # 水平标准化，-d row，可选column
    bash $sd/sp_pheatmap.sh \
      -f humann3/pathabundance_relab_unstratified.tsv \
      -H 'TRUE' -u 20 -v 50 \
      -Q group.txt -d row

## 物种kraken2

### Alpha多样性

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # 查看帮助
    Rscript $sd/kraken2alpha.R -h    
    # -d指定最小样本量，默认0为最小值，抽平文件tax_norm.txt，alpha多样性tax_alpha.txt
    Rscript $sd/kraken2alpha.R \
      --input kraken2/tax_count.mpa \
      --depth 0 \
      --species kraken2/tax_count.txt \
      --normalize kraken2/tax_count.norm \
      --output kraken2/tax_count.alpha
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript $sd/alpha_boxplot.R -h # 查看参数
    Rscript $sd/alpha_boxplot.R \
      -i kraken2/tax_count.alpha \
      -a shannon \
      -d metadata.txt \
      -n Group \
      -o kraken2/ \
      -w 89 -e 59
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript $sd/alpha_boxplot.R -i kraken2/tax_count.alpha -a ${i} \
      -d metadata.txt -n Group -w 89 -e 59 \
      -o kraken2/
    done
    
### 热图

    # 转换为metaphalan2 spf格式，kpcofg分隔符为下划线“__”
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="unclassified";a["p"]="unclassified";a["c"]="unclassified";a["o"]="unclassified";a["f"]="unclassified";a["g"]="unclassified";a["s"]="unclassified"; \
      split($1,x,"|");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"],$0;}' \
      kraken2/tax_count.txt > temp.txt
    cut -f 1-7,9- temp.txt > kraken2/tax_count.spf
    sed -i '1 s/unclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/' \
      kraken2/tax_count.spf
    # 绘制热图，可选界、门、纲、目、科、属、种(Kingdom	Phylum	Class	Order	Family	Genus	Species)
    tax=Genus
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i kraken2/tax_count.spf \
      -t ${tax} \
      -n 25 \
      -w 118 -e 118 \
      -o kraken2/heatmap_${tax}

### 箱线图

    # 绘制属水平Top30箱线图
    Rscript $sd/metaphlan_boxplot.R \
      -i kraken2/tax_count.spf \
      -t Genus \
      -n 30 \
      -o kraken2/boxplot_Genus
    # 绘制门水平Top10箱线图
    Rscript $sd/metaphlan_boxplot.R \
      -i kraken2/tax_count.spf \
      -t Phylum \
      -n 10 -w 6 -e 4 \
      -o kraken2/boxplot_Phylum

## 物种kraken2-braken2

### Alpha多样性

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # 查看帮助
    Rscript $sd/otutab_rare.R -h    
    # -d指定最小样本量，默认0为最小值，抽平文件tax_norm.txt，alpha多样性tax_alpha.
    tax=S
    Rscript $sd/otutab_rare.R \
      --input kraken2/bracken.${tax}.0.2 \
      --depth 0 --seed 1 \
      --normalize kraken2/bracken.${tax}.norm \
      --output kraken2/bracken.${tax}.alpha
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript $sd/alpha_boxplot.R -h # 查看参数
    mkdir -p kraken2/${tax}
    Rscript $sd/alpha_boxplot.R \
      -i kraken2/bracken.${tax}.alpha \
      -a shannon -w 89 -e 59 \
      -d metadata.txt -n Group \
      -o kraken2/${tax}
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript $sd/alpha_boxplot.R -i kraken2/bracken.${tax}.alpha -a ${i} \
      -d metadata.txt -n Group -w 89 -e 59 \
      -o kraken2/${tax}
    done
    
### Beta多样性

    # Beta多样性距离矩阵计算
    # Mac 用不了面对的 usearch，使用在线平台 https://www.bic.ac.cn/BIC
    mkdir -p kraken2/beta/
    $sd/../win/usearch -beta_div kraken2/bracken.${tax}.norm \
        -filename_prefix kraken2/beta/

    # PCoA分析输入文件，选择分组，输出文件，图片尺寸mm，统计见beta_pcoa_stat.txt
    # 可选距离有 bray_curtis, euclidean, jaccard, manhattan
    dis=bray_curtis
    Rscript $sd/beta_pcoa.R \
      --input kraken2/beta/${dis}.txt \
      --design metadata.txt \
      --group Group \
      --width 89 --height 59 \
      --output kraken2/pcoa.${dis}.pdf 
      
### 堆叠柱状图

    # 以门(P)/种(S)水平为例，结果包括output.sample/group.pdf两个文件
    tax=P
    Rscript ${sd}/tax_stackplot.R \
      --input kraken2/bracken.${tax}.0.2 --design metadata.txt \
      --group Group --output kraken2/bracken.${tax}.stackplot \
      --legend 10 --width 89 --height 59

## 基因组进化树注释table2itol

    cd ../binning/result/itol/
    ## 方案1. 分类彩带、数值热图、种标签
    # -a 找不到输入列将终止运行（默认不执行）-c 将整数列转换为factor或具有小数点的数字，-t 偏离提示标签时转换ID列，-w 颜色带，区域宽度等， -D输出目录，-i OTU列名，-l 种标签替换ID
    Rscript ${sd}/table2itol.R -a -c double -D plan1 -i ID -l Species -t %s -w 0.5 annotation.txt
    # 生成注释文件中每列为单独一个文件

    ## 方案2. 数值柱形图，树门背景色，属标签
    Rscript ${sd}/table2itol.R -a -d -c none -D plan2 -b Phylum -i ID -l Genus -t %s -w 0.5 annotation.txt

    ## 方案3.分类彩带、整数为柱、小数为热图
    Rscript ${sd}/table2itol.R -c keep -D plan3 -i ID -t %s annotation.txt

    ## 方案4. 将整数转化成因子生成注释文件
    Rscript ${sd}/table2itol.R -a -c factor -D plan4 -i ID -l Genus -t %s -w 0 annotation.txt

# 附录
 
## windows换行符处理
    
    # 查看行尾是否有^M
    cat -A $sd/metaphlan_hclust_heatmap.R | head
    # 转换Windows为Linux换行符
    dos2unix $sd/metaphlan_hclust_heatmap.R   

## 表格合并

    # 按末列注释，必须名叫KO
    Rscript $sd/mat_gene2ko.R  \
      -i temp/dbcan2/gene_fam.count \
      -o dbcan2/fam_merge.count
