# R语言统计绘图

    # 进入工作目录
    cd /c/meta
    # 设置脚本所在目录(Script Directory)
    sd=/c/public/script

## 2.3 HUMAnN2

### 2.3.1 物种组成热图

    # 显示脚本帮助 help
    Rscript $sd/metaphlan_hclust_heatmap.R -h
    # 按指定列合并、排序并取Top25种绘制热图
    # -i输入MetaPhlAn2文件；
    # -t 分类级别，可选Kingdom/Phylum/Class/Order/Family/Genus/Species/Strain，界门纲目科属种株，推荐门，目，属
    # -n 输出物种数量，默认为25，最大为合并后的数量
    # -o输出图表前缀，默认根据输入文件、物种级别和数量自动生成；
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i result/metaphlan2/taxonomy.spf \
      -t Family -n 25 \
      -w 183 -e 118 \
      -o result/metaphlan2/heatmap_Family

    # 属水平的Top30
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i result/metaphlan2/taxonomy.spf \
      -t Genus -n 30 \
      -w 183 -e 118 \
      -o result/metaphlan2/heatmap_Genus


## 2.7 kraken2物种注释


### 2.7.3 物种多样性分析

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # 查看帮助
    Rscript $sd/otutab_rare.R -h    
    # -d指定最小样本量，默认0为最小值，抽平文件tax_norm.txt，alpha多样性tax_alpha.
    tax=S
    Rscript $sd/otutab_rare.R \
      --input result/kraken2/bracken.${tax}.txt \
      --depth 0 --seed 1 \
      --normalize result/kraken2/bracken.${tax}.norm \
      --output result/kraken2/bracken.${tax}.alpha
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript $sd/alpha_boxplot.R -h # 查看参数
    mkdir -p result/kraken2/${tax}
    Rscript $sd/alpha_boxplot.R \
      -i result/kraken2/bracken.${tax}.alpha \
      -a shannon \
      -d result/metadata.txt \
      -n Group \
      -o result/kraken2/${tax} \
      -w 89 -e 59
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript $sd/alpha_boxplot.R -i result/kraken2/bracken.${tax}.alpha -a ${i} \
      -d result/metadata.txt -n Group -w 89 -e 59 \
      -o result/kraken2/${tax}
    done

    # Beta多样性距离矩阵计算
    mkdir -p result/kraken2/beta/
    usearch -beta_div result/kraken2/bracken.${tax}.norm \
        -filename_prefix result/kraken2/beta/

    # PCoA分析输入文件，选择分组，输出文件，图片尺寸mm，统计见beta_pcoa_stat.txt
    # 可选距离有 bray_curtis, euclidean, jaccard, manhattan
    dis=bray_curtis
    Rscript $sd/beta_pcoa.R \
      --input result/kraken2/beta/${dis}.txt \
      --design result/metadata.txt \
      --group Group \
      --width 89 --height 59 \
      --output result/kraken2/pcoa.${dis}.pdf 
      
### 2.7.4 物种组成

    # 以门(P)/种(S)水平为例，结果包括output.sample/group.pdf两个文件
    tax=S
    Rscript ${sd}//tax_stackplot.R \
      --input result/kraken2/bracken.${tax}.txt --design result/metadata.txt \
      --group Group --output result/kraken2/bracken.${tax}.stackplot \
      --legend 8 --width 89 --height 59


# 附录

## windows换行符处理
    
    # 查看行尾是否有^M
    cat -A $sd/metaphlan_hclust_heatmap.R | head
    # 转换Windows为Linux换行符
    dos2unix $sd/metaphlan_hclust_heatmap.R   

## 表格合并(按末列注释，必须名叫KO)

    Rscript $sd/mat_gene2ko.R  \
      -i temp/dbcan2/gene_fam.count \
      -o result/dbcan2/fam_merge.count

## kraken2的alpha多样性

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # 查看帮助
    Rscript $sd/kraken2alpha.R -h    
    # -d指定最小样本量，默认0为最小值，抽平文件tax_norm.txt，alpha多样性tax_alpha.txt
    Rscript $sd/kraken2alpha.R \
      --input result/kraken2/tax_count.mpa \
      --depth 0 \
      --species result/kraken2/tax_count.txt \
      --normalize result/kraken2/tax_count.norm \
      --output result/kraken2/tax_count.alpha
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript $sd/alpha_boxplot.R -h # 查看参数
    Rscript $sd/alpha_boxplot.R \
      -i result/kraken2/tax_count.alpha \
      -a shannon \
      -d result/metadata.txt \
      -n Group \
      -o result/kraken2/ \
      -w 89 -e 59
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript $sd/alpha_boxplot.R -i result/kraken2/tax_count.alpha -a ${i} \
      -d result/metadata.txt -n Group -w 89 -e 59 \
      -o result/kraken2/
    done
    
## kraken2的物种组成热图和箱线图

    # 转换为metaphalan2 spf格式，分隔符为下划线“_”
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["d"]="unclassified";a["p"]="unclassified";a["c"]="unclassified";a["o"]="unclassified";a["f"]="unclassified";a["g"]="unclassified";a["s"]="unclassified"; \
      split($1,x,"|");for(i in x){split(x[i],b,"_");a[b[1]]=b[2];} \
      print a["d"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"],$0;}' \
      result/kraken2/tax_count.txt > temp/temp.txt
    cut -f 1-7,10- temp/temp.txt > result/kraken2/tax_count.spf
    sed -i '1 s/unclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified/Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/' \
      result/kraken2/tax_count.spf
    # 绘制热图，可选域、门、纲、目、科、属、种(Domain	Phylum	Class	Order	Family	Genus	Species)
    tax=Species
    Rscript $sd/metaphlan_hclust_heatmap.R \
      -i result/kraken2/tax_count.spf \
      -t ${tax} \
      -n 25 \
      -w 118 -e 118 \
      -o result/kraken2/heatmap_${tax}

    # 绘制属水平Top30箱线图
    Rscript $sd/metaphlan_boxplot.R \
      -i result/kraken2/tax_count.spf \
      -t Genus \
      -n 30 \
      -o result/kraken2/boxplot_Genus
    # 绘制门水平Top10箱线图
    Rscript $sd/metaphlan_boxplot.R \
      -i result/kraken2/tax_count.spf \
      -t Phylum \
      -n 10 -w 6 -e 4 \
      -o result/kraken2/boxplot_Phylum
