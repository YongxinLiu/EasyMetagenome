[TOC]

# 2EasyMetagenome Visualization (2易宏基因组统计与可视化)

    # Authors(作者): Yong-Xin Liu(刘永鑫), Defeng Bai(白德凤), Tong Chen(陈同) et al.
    # Version(版本): 1.25, 2026/8/4
    # Operation System(操作系统): Linux Ubuntu 22.04+ / CentOS 7.7+ 
    # Homepage(主页): https://github.com/YongxinLiu/EasyMetagenome
    # Cititon(引文): Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

    # Set work directory(设置工作目录)
    wd=/d/Download/github/EasyMetagenome/result
    sd=/d/Download/github/EasyMicrobiome/script
    PATH=$PATH:$sd/../win:$sd
    cd $wd

## MetaPhlAn4 taxonomic composition (物种组成)

### Alpha diversity (α多样性)

    # index calculation (多样性指数计算)
    sed -i '/^#/d' metaphlan4/taxonomy.tsv
    Rscript ${sd}/metaphlan4_alpha.R -h
    Rscript $sd/metaphlan4_alpha.R \
      -i metaphlan4/taxonomy.tsv \
      -g metadata.txt \
      -t 7 \
      -o metaphlan4/alpha

    # Plot alpha diversity boxplot (绘制多样性指数箱线图)
    # result: metaphlan4/boxplot_{richness/shannon/shannon/invsimpson/Pielou_evenness}.pdf, sig using different letter a/b/c
    Rscript $sd/alpha_boxplot.R -h
    head -n1 metaphlan4/alpha.txt
    Rscript $sd/alpha_boxplot.R \
      -i metaphlan4/alpha.txt \
      -a shannon \
      -d metadata.txt \
      -n Group \
      -o metaphlan4/ \
      -w 89 -e 59
  
    # 5 Alpha diversity boxplot
    for i in `head -n1 metaphlan4/alpha.txt|cut -f 2-`;do
    Rscript $sd/alpha_boxplot.R -i metaphlan4/alpha.txt -a ${i} \
      -d metadata.txt -n Group -w 89 -e 59 \
      -o metaphlan4/
    done

    # Alpha diversity boxplot with pvalue
    # result: p_boxplot_${alpha_diversity}.pdf
    Rscript $sd/alpha_boxplot_new.R \
      -i metaphlan4/alpha.txt \
      -a shannon \
      -d metadata.txt \
      -n Group \
      -o metaphlan4/ \
      -w 49 -e 79
  
    # 5 Alpha diversity boxplot
    for i in `head -n1 metaphlan4/alpha.txt|cut -f 2-`;do
    Rscript $sd/alpha_boxplot_new.R -i metaphlan4/alpha.txt -a ${i} \
      -d metadata.txt -n Group -w 49 -e 79 \
      -o metaphlan4/
    done

    # Venn diagram (维恩图)
    # Filter abundance > 0.5 in each sample 筛选每个样本>0.5%的分类单元，包括界门纲目科属种
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
       else {for(i=2;i<=NF;i++) if($i>0.5) print $1, a[i];}}' \
       metaphlan4/taxonomy.tsv > metaphlan4/taxonomy_high.tsv
    wc -l metaphlan4/taxonomy_high.tsv
    # EVenn (https://www.bic.ac.cn/test/venn/#/) supports 6 groups and real-time viewing of element intersections (支持实时查看元素交集)
    # Ref：Yang, Tong Chen, Yong-Xin Liu, Luqi Huang. 2024. Visualizing set relationships: EVenn's comprehensive approach to Venn diagrams. iMeta 3: e184. https://doi.org/10.1002/imt2.184
    # script support 2-5 group: -f input, -a/b/c/d/g group, -w/u width and length, -p output filename
    # 本地2-5组比较:-f输入文件,-a/b/c/d/g分组名,-w/u为宽高英寸,-p输出文件名后缀
    bash ${sd}/sp_vennDiagram.sh -f metaphlan4/taxonomy_high.tsv \
      -a C1 -b C2 -c Y1 -d Y2 -g Y3 \
      -w 4 -u 4 \
      -p C1_C2_Y1_Y2_Y3

    # Group mean (求均值再两组比较)
    Rscript ${sd}/otu_mean.R --input metaphlan4/taxonomy.tsv \
      --metadata metadata.txt \
      --group Group --thre 0 \
      --scale F --all TRUE --type mean \
      --output metaphlan4/group_mean.txt    
    # Filter abundance > 0.5 in each group (筛选每个组>0.5%的分类单元，包括界门纲目科属种)
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
       else {for(i=2;i<=NF;i++) if($i>0.5) print $1, a[i];}}' \
       metaphlan4/group_mean.txt > metaphlan4/group_high.tsv
    bash ${sd}/sp_vennDiagram.sh -f metaphlan4/group_high.tsv \
      -a Centenarians -b Young -c All \
      -w 8 -u 6 \
      -p Centenarians_Young_All      

### Beta diversity (β多样性)
  
    # Distance calculation (距离计算)
    # -t taxonomy level: 1 kingdom; 2 phylum; 3 order; 4 class; 5 family; 6 genus; 7 species
    # -m distance type："bray", "euclidean", "jaccard", "manhattan" etc.
    # e.g. Claculate level 7 (species) in bray-curtis distance
    Rscript ${sd}/metaphlan4_beta.R -h
    Rscript $sd/metaphlan4_beta.R \
      -i metaphlan4/taxonomy.tsv \
      -g metadata.txt \
      -t 7 \
      -m bray \
      -o metaphlan4/beta

    # Beta diversity in Principal Coordinates Analysis (PCoA) 
    # input distance matrix, group name, width/height in mm, stat in beta_pcoa_stat.txt
    Rscript $sd/beta_pcoa.R \
      --input metaphlan4/beta_bray.txt \
      --design metadata.txt \
      --group Group \
      --width 89 --height 59 \
      --output metaphlan4/pcoa.bray_curtis.pdf


### Taxonomic composition (物种组成)

    # Heatmap (热图)
    # -h: Show help (显示脚本帮助)
    Rscript ${sd}/metaphlan_hclust_heatmap.R -h
    # stat  (显示脚本帮助)
    csvtk -t stat metaphlan4/taxonomy.spf
    # Taxonomy level: Order, Top n taxa to plot heatmap (指定分类级、取前n种绘制热图）
    # Order Top 20, Family Top 25, Genus Top 30
    tax=Order
    n=25
    # -i: input metaphlan4 result in spf format
    # -t: taxa in Kingdom/Phylum/Class/Order/Family/Genus/Species/Strain (界门纲目科属种株)，推荐门，目，属
    # -o: output figure in pdf, data in txt, default name Heatmap+($tax).pdf/txt
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i metaphlan4/taxonomy.spf \
      -t ${tax} -n ${n} \
      -w 183 -e 118 \
      -o metaphlan4/heatmap${tax}

    # Boxplot (箱线图)
    for tax in Phylum Family Genus Species; do
    Rscript $sd/metaphlan_boxplot.R \
          -i metaphlan4/taxonomy.spf \
          -t ${tax} \
          -n ${n} \
          -o metaphlan4/boxplot_${tax};done
      
    # Compare group boxplot *_compare.pdf (组间比较箱线图)
    for tax in Phylum Family Genus Species; do
    Rscript $sd/metaphlan4_boxplot_compare.R \
          -i metaphlan4/taxonomy.spf \
          -t ${tax} \
          -n 30 \
          -o metaphlan4/boxplot_${tax};done

    # Plot stackplot (绘制不同分类层级堆叠柱状图)
    # Prepare file (准备每个分类层级的文件)
    bash ${sd}/stack_data_prepare.sh
    for tax in Phylum Genus Species; do
    Rscript ${sd}/tax_stackplot.R \
          --input metaphlan4/${tax}.txt --design metadata.txt \
          --group Group --output metaphlan4/${tax}.stackplot \
          --legend 10 --width 120 --height 70; done

    # Sorted stackplot (排序分面堆叠柱状图)
    for tax in Phylum Genus Species; do
    Rscript ${sd}/tax_stackplot_order.R \
          --input metaphlan4/${tax}.txt --design metadata.txt \
          --group Group --output metaphlan4/${tax}.stackplot \
          --legend 10 --width 120 --height 70; done

### Different compare (差异比较)

    ## STAMP
    # If you encounter errors or no results, lower the threshold and p-value. 
    # The example sample size is small cause no significant, change it to P<0.1.
    # 如果遇见报错无结果调低threshold和pvalue值，此处样本少无显著改为P<0.1
    # Set compared groups, using '-' as seperate
    compare="Centenarians-Young"
    Rscript ${sd}/compare_stamp.R \
          --input metaphlan4/Genus.txt --metadata metadata.txt \
          --group Group --compare ${compare} --threshold 0.01 \
          --method "t.test" --pvalue 0.1 --fdr "none" \
          --width 300 --height 150 \
          --output metaphlan4/stamp_${compare}

    ### MaAsLin2差异物种分析火山图
    # 同样适用于功能通路差异分析
    Rscript ${sd}/compare_MaAsLin2.R \
          --input metaphlan4/Species.txt \
          --metadata metadata.txt \
          --output metaphlan4/
      
    # 根据差异分析结果绘制火山图
    Rscript ${sd}/compare_volcano.R \
          --input metaphlan4/MaAsLin2_overall_difference.csv \
          --trend metaphlan4/MaAsLin2_enriched_depleted.csv \
          --width 100 --height 80 \
          --output metaphlan4/
          
### Network analysis (网络分析)

    # 数据准备
    Rscript ${sd}/SparCC_data_processing.R \
          --input metaphlan4/Species.txt \
          --group metadata.txt \
          --output metaphlan4/
          
    # SparCC相关性(Correlation)和显著性(p value)计算
    chmod +x ${sd}/run_sparcc3.sh
    ## 运行如果出现Sparcity assumption violated的警告信息，可能是样本量较小导致的
    ## Centenarians
    bash ${sd}/run_sparcc3.sh \
          ./metaphlan4/Centenarians_sparcc.txt \
          ./SparCC3/R_Centenarians.txt \
          ./SparCC3/P_Centenarians.txt
          
    ## Young
    bash ${sd}/run_sparcc3.sh \
          ./metaphlan4/Young_sparcc.txt \
          ./SparCC3/R_Young.txt \
          ./SparCC3/P_Young.txt
    
    # 网络分析可视化，可自定义指定分组
    # 方式一：基于Gephi软件进行可视化 https://gephi.org/
    # 获取.graphml文件
    # Centenarians
    Rscript ${sd}/get_graphml.R \
         --Correlation ${wd}/SparCC3/R_Centenarians.txt \
         --Pvalue ${wd}/SparCC3/P_Centenarians.txt \
         --Group Centenarians \
         --output ${wd}/metaphlan4/     
    # Young
    Rscript ${sd}/get_graphml.R \
         --Correlation ${wd}/SparCC3/R_Young.txt \
         --Pvalue ${wd}/SparCC3/P_Young.txt \
         --Group Young \
         --output ${wd}/metaphlan4/     
         
    # 下载安装Gephi软件在本地带入.graphml文件进行可视化展示
    
    # 方式二：基于R代码进行可视化
    # Centenarians
    Rscript ${sd}/SparCC_visualization.R \
          --Correlation ${wd}/SparCC3/R_Centenarians.txt \
          --Pvalue ${wd}/SparCC3/P_Centenarians.txt \
          --Group Centenarians \
          --output ${wd}/metaphlan4/      
    # Young
    Rscript ${sd}/SparCC_visualization.R \
          --Correlation ${wd}/SparCC3/R_Young.txt \
          --Pvalue ${wd}/SparCC3/P_Young.txt \
          --Group Young \
          --output ${wd}/metaphlan4/      

## HUMAnN4 functional composition (功能组成)

### Clustering heatmap (分组聚类热图）

    bash $sd/sp_pheatmap.sh
    cut -f 1,3 metadata.txt > group.txt
    # -f输入文件，-H水平聚类，u/v图片宽/高，-P添加行注释文件，-Q添加列注释
    # 水平标准化，-d row，可选column
    bash $sd/sp_pheatmap.sh \
      -f humann4/path_relab_unstratified.tsv \
      -H 'TRUE' -u 20 -v 50 \
      -Q group.txt -d row
    # 结果为 输入文件名+pheamap.r/pdf，代码和图片


## Kraken2 Taxonomic composition (物种组成)

### Alpha diversity (多样性)

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
    head -6 kraken2/tax_count.txt
    sed -i 's/\r//' kraken2/tax_count.*
    
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

### Taxonomic composition (物种组成)

    # 转换为metaphalan spf格式，分隔符为下划线“__”
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="unclassified";a["p"]="unclassified";a["c"]="unclassified";a["o"]="unclassified";a["f"]="unclassified";a["g"]="unclassified";a["s"]="unclassified"; \
      split($1,x,"|");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"],$0;}' \
      kraken2/tax_count.txt > temp.txt
    cut -f 1-7,9- temp.txt > kraken2/tax_count.spf
    sed -i '1 s/unclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/' \
      kraken2/tax_count.spf
    # 绘制热图，可选域、门、纲、目、科、属、种(Domain	Phylum	Class	Order	Family	Genus	Species)
    # 单个运行
    tax=Genus
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i kraken2/tax_count.spf \
      -t ${tax} \
      -n 25 \
      -w 118 -e 118 \
      -o kraken2/heatmap_${tax}

    # 批量运行  
    for tax in Phylum Family Genus Species; do
    Rscript $sd/metaphlan_hclust_heatmap.R \
          -i kraken2/tax_count.spf \
          -t ${tax} \
          -n 10 \
          -w 118 -e 118 \
          -o kraken2/heatmap_${tax};done

### Taxonomic distribution boxplot (物种分布箱线图)

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
    # Batch other taxonomic level (批量绘制不同层级箱线图)
    for tax in Class Order Family Species; do
    Rscript $sd/metaphlan_boxplot.R \
          -i kraken2/tax_count.spf \
          -t ${tax} \
          -n 30 \
          -o kraken2/boxplot_${tax};done


## Braken2: Kraken2 taxonomic abundance reestimate (物种丰度重估)

### Alpha diversity (多样性)

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # 查看帮助
    Rscript $sd/otutab_rare.R -h    
    # -d指定最小样本量，默认0为最小值，抽平文件tax_norm.txt，alpha多样性tax_alpha.
    tax=S
    Rscript $sd/otutab_rare.R \
      --input kraken2/bracken.${tax}.txt \
      --depth 0 --seed 1 \
      --normalize kraken2/bracken.${tax}.norm \
      --output kraken2/bracken.${tax}.alpha
      
    # 批量
    for tax in D P G S; do
    Rscript $sd/otutab_rare.R \
          --input kraken2/bracken.${tax}.0.2 \
          --depth 0 --seed 1 \
          --normalize kraken2/bracken.${tax}.norm \
          --output kraken2/bracken.${tax}.alpha;done
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript $sd/alpha_boxplot.R -h # 查看参数
    Rscript $sd/alpha_boxplot.R \
      -i kraken2/bracken.${tax}.alpha \
      -a richness \
      -d metadata.txt \
      -n Group \
      -o kraken2/${tax} \
      -w 89 -e 59
    
    # 批量运行
    for tax in S; do
    for i in richness chao1 ACE shannon simpson invsimpson; do
    mkdir -p kraken2/${tax}/alpha/
    Rscript $sd/alpha_boxplot.R \
          -i kraken2/bracken.${tax}.alpha \
          -a ${i} \
          -d metadata.txt \
          -n Group \
          -o kraken2/${tax}/alpha/ \
          -w 89 -e 59; done
    done
    mv alpha_boxplot_TukeyHSD.txt  kraken2/${tax}_bracken_alpha

## Beta diversity (β多样性)

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
    
    # 批量运行  
    for tax in S; do
    mkdir -p kraken2/${tax}/beta/
    $sd/../win/usearch -beta_div kraken2/bracken.${tax}.norm     \
        -filename_prefix kraken2/${tax}/beta/
    for dis in bray_curtis euclidean jaccard manhatten; do
        Rscript $sd/beta_pcoa.R \
          --input kraken2/${tax}/beta/${dis}.txt \
          --design metadata.txt \
          --group Group \
          --width 89 --height 59 \
          --output kraken2/${tax}/beta/pcoa.${dis}.pdf; done
    mv beta_pcoa_stat.txt kraken2/${tax}/beta/;done
    

### 堆叠柱状图

    # 以门(P)/种(S)水平为例，结果包括output.sample/group.pdf两个文件
    tax=P
    Rscript ${sd}/tax_stackplot.R \
      --input kraken2/bracken.${tax}.norm --design metadata.txt \
      --group Group --output kraken2/bracken.${tax}.stackplot \
      --legend 10 --width 89 --height 59
  
    # 批量运行
    for tax in D G P S; do
    Rscript ${sd}/tax_stackplot.R \
      --input kraken2/bracken.${tax}.0.2-H --design metadata.txt \
      --group Group --output kraken2/bracken.${tax}.stackplot \
      --legend 10 --width 120 --height 70; done


### STAMP图

    compare="Centenarians-Young"
    tax=S
    Rscript ${sd}/compare_stamp.R \
      --input kraken2/bracken.${tax}.0.2-H --metadata metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method "wilcox.test" --pvalue 0.1 --fdr "none" \
      --width 89 --height 59 \
      --output kraken2/stamp_${compare}


## 基因组进化树注释table2itol

    #cd ../binning/result/itol/
    cd ../result/itol/
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


## MAGs分析

### MAGs与样本稀疏曲线

    cd ${wd}
    Rscript ${sd}/MAG_sample_rare.R \
      --input ./coverm/abundance.tsv \
      --taxonomy ./itol/tax.txt \
      --output ./coverm/


### MAGs不同功能数据库注释结果比较UpSet图
    
    # 数据处理得到data_venn1.txt
    bash $sd/extract_eggnog_anno01.sh
    # 绘图
    Rscript ${sd}/MAG_protiens_UpSet.R \
        --i eggnog/data_venn1.txt \
        --o eggnog/
        
    ## 单个功能数据库分组比较(COGs, KEGG, ECs, GOs), 这里以COGs注释结果作为示例
    # 分组比较UpSet图
    # 处理数据得到data_venn2.txt
    bash $sd/extract_eggnog_anno02.sh 
    # 绘图
    Rscript ${sd}/common_upset.R \
        --i eggnog/data_venn2.txt \
        --o eggnog/

    # 提取交集并绘制交集组成饼图，支持自定义分组
    Rscript ${sd}/upset_intersection_composition.R \
        --input eggnog/data_venn2.txt \
        --groups Centenarians,Young \
        --output eggnog/


### MAGs功能通路桑基图

    # 需要提前对感兴趣的功能通路进行分析获得data_sankey.txt文件
    # 下载KEGG Orthology （KO）分层文件
    bash $sd/KEGG_Orthology_data_extract.sh
    # 选择感兴趣的或者于研究项目密切相关的Level3通路
    # 比如这里选择Tryptophan metabolism，Butanoate metabolism，Sulfur metabolism这3个Level3通路作为案例
    # 从生成的KEGG_KO_Annotation.tsv文件中提取Level3_Name是Tryptophan metabolism，Butanoate metabolism，Sulfur metabolism的数据形成新的表格
    bash ${sd}/extract_pathways.sh \
        "Tryptophan metabolism" \
        "Butanoate metabolism" \
        "Sulfur metabolism"
    # 绘制桑基图
    Rscript ${sd}/function_sankey.R \
        --i humann4/data_sankey2.txt \
        --o humann4/


### MAGs完整性和污染率关系图

    # 生成绘图数据
    bash ${sd}/merged_tax.sh
    sed 's/ /_/g' checkm2/taxonomy_merge.txt > checkm2/taxonomy_merge2.txt
    
    # 绘图
    Rscript ${sd}/MAG_quality_taxonomy.R \
        --i checkm2/taxonomy_merge2.txt \
        --o checkm2/

### MAGs系统发育分析树状图

    # 肠道微生物样本示例，这里以百岁老人数据为例
    # 获取annotation.txt数据
    bash ${sd}/phylogenetic_anno_data.sh
    Rscript ${sd}/phylogenetic_tree.R \
        --input ../temp/gtdb_infer/tax.unrooted.tree \
        --annotation ../temp/gtdb_infer/annotation.txt \
        --output ../temp/gtdb_infer/

    # 环境样本示例
    # 由于1Pipeline.sh中只用了肠道微生物样本作为示例数据，这里仅提供参考
    # 实际运行代码和肠道微生物样本案例相似
    Rscript ${sd}/phylogenetic_tree_env.R \
        --input ../temp/gtdb_infer/tax.unrooted2.tree \
        --annotation ../temp/gtdb_infer/annotation3.txt \
        --output ../temp/gtdb_infer/
        
        
### Two Competing Guilds (TCG) Network

    # 参考Wu G, Xu T, Zhao N, et al. A core microbiome signature as an indicator of health. Cell. 2024.
    # 分析肠道核心菌群的“蹊跷板”竞争模式
    
    # HCC-TCG网络输入数据准备
    # 1.GTDBTK result：GTDBTK_tax.txt
    # 获取GTDBTK.txt文件
    bash ${sd}/GTDBTK_tax.sh \
        -i itol/tax.txt \
        -o itol/GTDBTK_tax.txt
    
    # FastSpar方法 (FastSpar方法比直接用SparCC计算更快, 部分步骤运行需要切换Bash (Windows Subsystem for Linux))
    ############################################################################
    # 2.FastSpar计算输入文件准备: MAG_for_FastSpar_Centenarians.txt; MAG_for_FastSpar_Young.txt
    # 数据准备,abundance_cutoff和prevalence_cutoff可根据实际进行调整
    Rscript ${sd}/FastSpar_input01.R \
          --abundance coverm/abundance.tsv \
          --taxonomy itol/tax.txt \
          --metadata metadata.txt \
          --group1 Centenarians \
          --group2 Young \
          --abundance_cutoff 0.0007 \
          --prevalence_cutoff 0 \
          --outdir ./coverm
          
    ## txt文件转为tsv文件; 最左上角head名称为#OTU ID, 符合输入数据运行规范
    chmod +x ${sd}/txt2tsv.sh
    # Centenarians
    bash ${sd}/txt2tsv.sh \
        ./coverm/MAG_for_FastSpar_Centenarians.txt \
        ./coverm/MAG_for_FastSpar_Centenarians_fixed.tsv
    # Young
    bash ${sd}/txt2tsv.sh \
        ./coverm/MAG_for_FastSpar_Young.txt \
        ./coverm/MAG_for_FastSpar_Young_fixed.tsv

    # 3.SparCC相关性(Correlation)和显著性(p value)计算
    # 进入Bash (Windows Subsystem for Linux)环境, Bash (Windows Subsystem for Linux)比Git Bash路径多/mnt, 需要修改路径重新运行后进入工作目录才能正常运行
    # 安装FastSpar
    conda create -n fastspar -c conda-forge -c bioconda fastspar
    conda activate fastspar
    fastspar --help
    
    # 计算相关性(Correlation)和显著性(p value)
    sed -i 's/\r$//' ${sd}/run_FastSpar.sh
    bash ${sd}/run_FastSpar.sh \
          -c $wd/coverm/MAG_for_FastSpar_Centenarians_fixed.tsv \
          -y $wd/coverm/MAG_for_FastSpar_Young_fixed.tsv \
          -o $wd/coverm/FastSpar_result \
          -t 6 \
          -b 1000 \
          -i 1000
          
    # 删除第一行开头的 #OTU ID 并输出到指定文件
    chmod +x ${sd}/remove_otu_header.sh
    sed -i 's/\r$//' ${sd}/remove_otu_header.sh
    # Centenarians correlation
    bash ${sd}/remove_otu_header.sh \
         $wd/coverm/FastSpar_result/Centenarians/Centenarians_correlation.tsv \
         $wd/coverm/R_Centenarians2.txt
    # Centenarians pvalue
    bash ${sd}/remove_otu_header.sh \
         $wd/coverm/FastSpar_result/Centenarians/Centenarians_pvalues.tsv \
         $wd/coverm/P_Centenarians2.txt
    # Young correlation
    bash ${sd}/remove_otu_header.sh \
         $wd/coverm/FastSpar_result/Young/Young_correlation.tsv \
         $wd/coverm/R_Young2.txt
    # Young pvalue
    bash ${sd}/remove_otu_header.sh \
         $wd/coverm/FastSpar_result/Young/Young_pvalues.tsv \
         $wd/coverm/P_Young2.txt
    ############################################################################
    
    # 切换到Git Bash环境中
    # 4.获取stable.edge.txt文件和stable.node.txt文件
    Rscript ${sd}/TCG_FromNetwork.R \
          --gtdb ${wd}/itol/GTDBTK_tax.txt \
          --group1_r ${wd}/coverm/R_Centenarians2.txt \
          --group1_p ${wd}/coverm/P_Centenarians2.txt \
          --group2_r ${wd}/coverm/R_Young2.txt \
          --group2_p ${wd}/coverm/P_Young2.txt \
          --output ${wd}/coverm
    
    # 5.TCG网络绘制
    # 下载安装Cytoscape软件网址https://cytoscape.org/
    # 将获取的stable.edge2.txt和stable.node2.txt导入Cytoscape
    # Import Network from File System导入stable.edge2.txt文件 -> 菜单栏Layout选择Edge-weighted Spring Embedded Layout -> 选择NorPv2
    # 左下角选中Edge对象，左上角选中style功能，Stroke Color选项Column选择NorP，Mapping Type选择Discrete Mapping，设置不同分组颜色
    
    # 6.绘制主要聚类群(Cluster)的层次聚类树
    # ccCluster.txt需要根据自己的数据手动生成, C1为ccCluster.txt中cluster列的列名
    # stable.node2.txt和stable.edge2.txt替换为上一步生成的结果
    # GTDBTK.txt替换为生成的数据
    Rscript ${sd}/cluster_from_Cytoscape.R \
          --ccCluster ${wd}/coverm/ccCluster.txt \
          --node ${wd}/coverm/stable.node2.txt \
          --edge ${wd}/coverm/stable.edge2.txt \
          --gtdb ${wd}/coverm/GTDBTK.txt \
          --cluster C1 \
          --output ${wd}/coverm


### 单菌基因组注释和绘图

# 注释：https://bakta.computational.bio/submit
# 注释和绘图：https://proksee.ca


### 宏基因组数据泛基因组分析绘图

# 绘图： https://anvi-server.org/


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


## CoverM result mean

    # 按组求均值，需要metadata中有3列且每个组有多个样本
    Rscript ${db}/EasyMicrobiome/script/otu_mean.R --input result/coverm/abundance.tsv \
      --metadata result/metadata.txt \
      --group Group --thre 0 \
      --scale TRUE --zoom 100 --all TRUE --type mean \
      --output result/coverm/group_mean.txt
    # https://www.bic.ac.cn/ImageGP/ 直接选择热图可视化


## graphlan_plot (物种组成)

    # Graphlan输入文件获取
    # Add heatmap as an example, if need to plot bar use "--type bar" instead
    Rscript ${sd}/graphlan_plot.R -h
    Rscript ${sd}/graphlan_plot.R \
      --input metaphlan4/taxonomy.spf \
      --design metadata.txt \
      --type heatmap \
      --output metaphlan4/graphlan
      
      
### MAGs耐药基因丰度热图

    Rscript ${sd}/functional_gene_heatmap.R \
        --input binning/temp/ARG2/ARG_final.txt \
        --group binning/temp/ARG2/group_ARG.txt \
        --output binning/temp/      
      
