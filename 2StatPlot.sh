[TOC]

# 易宏基因组流程EasyMetagenomePipeline

```
# 版本: 1.22, 2024/11/26
# 测试环境为Windows 10+ / MacOS 10+

# 设置结果目录(通常为项目中的result，此处为12个样本结果result12)
# 可选癌症cancer或年龄age各18个样本数据示例
# mac/linux 需要修改路径格式如下， ~ 代表家目录
# wd=~/meta/result12
# sd=~/EasyMicrobiome/script
wd=/d/EasyMetagenome/result12
# 设置脚本所在目录(Script Directory)，系统为win/mac/linux
sd=/d/EasyMicrobiome/script
PATH=$PATH:$sd/../win:$sd
# 进入结果目录
cd $wd

```

## 物种Metaphlan4

```
### Alpha多样性指数计算(Metaphlan4)
Rscript ${sd}/metaphlan4_alpha.R -h
Rscript $sd/metaphlan4_alpha.R \
  -i metaphlan4/taxonomy.tsv \
  -g metadata.txt \
  -t 7 \
  -o metaphlan4/alpha
  
### 绘制alpha多样性指数箱线图
# 绘制Alpha多样性指数，结果为输入文件+类型richness/shannon/shannon/invsimpson/Pielou_evenness
# Rscript $sd/alpha_boxplot.R -h # 查看参数
# 样式1：用字母a,b标明显著性
Rscript $sd/alpha_boxplot.R \
  -i metaphlan4/alpha.txt \
  -a shannon \
  -d metadata.txt \
  -n Group \
  -o metaphlan4/ \
  -w 89 -e 59
  
### 批量计算6种指数的箱线图+统计
for i in observed_species shannon simpson invsimpson Pielou_evenness;do
Rscript $sd/alpha_boxplot.R -i metaphlan4/alpha.txt -a ${i} \
  -d metadata.txt -n Group -w 89 -e 59 \
  -o metaphlan4/
done

# 样式2：用p值标明显著性  
Rscript $sd/alpha_boxplot_new.R \
  -i metaphlan4/alpha.txt \
  -a shannon \
  -d metadata.txt \
  -n Group \
  -o metaphlan4/ \
  -w 49 -e 79
  
# 批量计算6种指数的箱线图+统计
for i in observed_species shannon simpson invsimpson Pielou_evenness;do
Rscript $sd/alpha_boxplot_new.R -i metaphlan4/alpha.txt -a ${i} \
  -d metadata.txt -n Group -w 49 -e 79 \
  -o metaphlan4/
done
  


### Beta多样性距离矩阵计算(Metaphlan4)
# 可选计算分类级别-t：1-界；2-门；3-纲；4-目；5-科；6-属；7-种
# 可选距离-m："bray", "euclidean", "jaccard", "manhattan"等
# 这里以物种水平和bray-curtis距离举例
Rscript ${sd}/metaphlan4_beta.R -h
Rscript $sd/metaphlan4_beta.R \
  -i metaphlan4/taxonomy.tsv \
  -g metadata.txt \
  -t 7 \
  -m bray \
  -o metaphlan4/beta

### Beta多样性PCoA分析  
# PCoA分析输入文件，选择分组，输出文件，图片尺寸mm，统计见beta_pcoa_stat.txt
# 可选距离有 bray_curtis, euclidean, jaccard, manhattan
# 此处可能报错“duplicated row.names”,需要给beta.txt增加行名后运行
Rscript $sd/beta_pcoa.R \
  --input metaphlan4/beta_bray.txt \
  --design metadata.txt \
  --group Group \
  --width 89 --height 59 \
  --output metaphlan4/pcoa.bray_curtis.pdf 
```


### 热图

```
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

```

### 不同分类级别树状图

```
# 此步骤需要在linux环境中运行，请在linux环境中安装humann2，在humann2环境下运行
# 处理 taxonomy.spf 文件，生成 taxonomy_modified.spf
${sd}/taxonomy_modified.sh
#整体结构图
#安装humann2
pip install humann2
conda activate humann2
Rscript ${sd}/graphlan_plot_new.r --input metaphlan4/taxonomy_modified.spf\
		--design metadata.txt --type heatmap --output metaphlan4/Heat_Structures
```


### 维恩图

```
# 筛选每个样本>0.5%的分类单元，包括界门纲目科属种
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
   else {for(i=2;i<=NF;i++) if($i>0.5) print $1, a[i];}}' \
   metaphlan4/taxonomy.tsv > metaphlan4/taxonomy_high.tsv
wc -l metaphlan4/taxonomy_high.tsv
# 在线绘制，支持实时查看元素交集 http://www.ehbio.com/test/venn/#/
# 引文：Mei Yang, Tong Chen, Yong-Xin Liu, Luqi Huang. 2024. Visualizing set relationships: 
# EVenn's comprehensive approach to Venn diagrams. iMeta 3: e184. https://doi.org/10.1002/imt2.184
# 本地5组比较:-f输入文件,-a/b/c/d/g分组名,-w/u为宽高英寸,-p输出文件名后缀
bash ${sd}/sp_vennDiagram.sh -f metaphlan4/taxonomy_high.tsv \
  -a C1 -b C2 -c N1 -d N2 -g N3 \
  -w 4 -u 4 \
  -p C1_C2_N1_N2_N3

# 求均值再两组比较
sed -i '/^#/d' metaphlan4/taxonomy.tsv
Rscript ${sd}/otu_mean.R --input metaphlan4/taxonomy.tsv \
  --metadata metadata.txt \
  --group Group --thre 0 \
  --scale F --all TRUE --type mean \
  --output metaphlan4/group_mean.txt    
# 筛选每个组>0.5%的分类单元，包括界门纲目科属种
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
   else {for(i=2;i<=NF;i++) if($i>0.5) print $1, a[i];}}' \
   metaphlan4/group_mean.txt > metaphlan4/group_high.tsv
bash ${sd}/sp_vennDiagram.sh -f metaphlan4/group_high.tsv \
  -a Cancer -b Normal -c All \
  -w 4 -u 4 \
  -p Cancer_Normal_All      

```

### 丰度箱线图(metaphlan4)

```
# 整体丰度箱线图
for tax in Phylum Family Genus Species; do
Rscript $sd/metaphlan_boxplot.R \
      -i metaphlan4/taxonomy.spf \
      -t ${tax} \
      -n 30 \
      -o metaphlan4/boxplot_${tax};done
      
# 组间比较箱线图
for tax in Phylum Family Genus Species; do
Rscript $sd/metaphlan4_boxplot_compare.R \
      -i metaphlan4/taxonomy.spf \
      -t ${tax} \
      -n 30 \
      -o metaphlan4/boxplot_${tax};done
```

### 堆叠柱状图

```
# 准备每个分类层级的文件
chmod +x ${sd}/stack_data_prepare.sh
${sd}/stack_data_prepare.sh

# 绘制不同分类层级堆叠柱状图
for tax in Kingdom Phylum Genus Species; do
Rscript ${sd}/tax_stackplot.R \
      --input metaphlan4/${tax}.txt --design metadata.txt \
      --group Group --output metaphlan4/${tax}.stackplot \
      --legend 10 --width 120 --height 70; done
      
# 排序分面堆叠柱状图
for tax in Kingdom Phylum Genus Species; do
Rscript ${sd}/tax_stackplot_order.R \
      --input metaphlan4/${tax}.txt --design metadata.txt \
      --group Group --output metaphlan4/${tax}.stackplot \
      --legend 10 --width 120 --height 70; done

```


### STAMP组间比较图
# 如果遇见报错，可尝试调整threshold和pvalue值
#Error in UseMethod("gather") :
#  no applicable method for 'gather' applied to an object of class "factor"
#Calls: %>% -> summarise -> group_by -> <Anonymous>
```
compare="Normal-Cancer"
Rscript ${sd}/compare_stamp.R \
      --input metaphlan4/Genus.txt --metadata metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method "t.test" --pvalue 0.2 --fdr "none" \
      --width 100 --height 300 \
      --output metaphlan4/stamp_${compare}
```


### MaAsLin2差异物种分析火山图

# 同样适用于功能通路差异分析

```
# 差异分析
Rscript ${sd}/compare_MaAsLin2.R \
      --i metaphlan4/Species.txt \
      --m metadata.txt \
      --o metaphlan4/
      
# 根据差异分析结果绘制火山图
Rscript ${sd}/compare_valcano.R \
      --i metaphlan4/MaAsLin2_overall_difference.csv \
      --d metaphlan4/MaAsLin2_enriched_depleted.csv \
      --o metaphlan4/

```


### 物种稀疏相关网络(SparCC)分析

```
## 处理物种数据用于后续SparCC网络分析
# 此处以其中一组进行举例，其它分组重复运行即可
Rscript ${sd}/SparCC_data_processing.R \
      --input result12/metaphlan4/Species.txt \
      --group result12/metadata.txt \
      --output result12/metaphlan4/

## SparCC相关性(Correlation)和显著性(p value)计算
# 以下命令行需在linux环境中运行
git clone https://github.com/JCSzamosi/SparCC3.git
cd SparCC3
mkdir -p data
cp result12/metaphlan4/Cancer_sparcc.txt SparCC3/data/
cp result12/metaphlan4/Normal_sparcc.txt SparCC3/data/
# Step 1 - Compute correlations
python SparCC.py data/Cancer_sparcc.txt -i 20 \
      --cor_file=example/basis_corr/sxtr_sparcc_Cancer.tsv \
      > example/basis_corr/sxtr_sparcc_Cancer.log
      
# Step 2 - Compute bootstraps
# 此处测试节省时间设置为100，建议设置为1000
mkdir -p example/pvals_cancer
python MakeBootstraps.py data/Cancer_sparcc.txt \
      -n 100 -t bootstrap_#.txt \
      -p example/pvals_cancer/ >> sxtr_sparcc_Cancer.log

# Step 3 - Compute p-values
# 如果上面设置1000，这里99更改为999
for n in {0..99}; do python SparCC.py example/pvals_cancer/bootstrap_${n}.txt -i 20 \
      --cor_file=example/pvals_cancer/bootstrap_cor_${n}.txt >> sxtr_sparcc_Cancer.log; done

# 如果最开始设置1000，这里100更改为1000      
python PseudoPvals.py example/basis_corr/sxtr_sparcc_Cancer.tsv example/pvals_cancer/bootstrap_cor_#.txt 100 \
      -o example/pvals_cancer/pvals.two_sided_cancer.txt \
      -t two_sided >> sxtr_sparcc_cancer.log
      
# step 4 - Rename file
mv example/pvals_cancer/pvals.two_sided_cancer.txt sxtr_pvals_cancer.two_sided.tsv
mv example/basis_corr/sxtr_sparcc_Cancer.tsv sxtr_cov_mat_Cancer.tsv

## 可视化
Rscript ${sd}/SparCC_visualization.R \
      --Correlation result12/metaphlan4/sxtr_cov_mat_Cancer.tsv \
      --Pvalue result12/metaphlan4/sxtr_pvals_cancer.two_sided.tsv \
      --output result12/metaphlan4/

```


## 功能HUMAnN3

### 分组聚类热图

```
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

```

## 物种kraken2

### Alpha多样性

```
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

```

### 热图

```
# 转换为metaphalan2 spf格式，分隔符为下划线“_”
awk 'BEGIN{OFS=FS="\t"}{delete a; a["d"]="unclassified";a["p"]="unclassified";a["c"]="unclassified";a["o"]="unclassified";a["f"]="unclassified";a["g"]="unclassified";a["s"]="unclassified"; \
  split($1,x,"|");for(i in x){split(x[i],b,"_");a[b[1]]=b[2];} \
  print a["d"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"],$0;}' \
  kraken2/tax_count.txt > temp.txt
cut -f 1-7,9- temp.txt > kraken2/tax_count.spf
sed -i '1 s/unclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified/Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/' \
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
```

### 箱线图(kraken2)

```
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
# 批量绘制不同层级箱线图
for tax in Phylum Family Genus Species; do
Rscript $sd/metaphlan_boxplot.R \
      -i kraken2/tax_count.spf \
      -t ${tax} \
      -n 30 \
      -o kraken2/boxplot_${tax};done

```

## 物种kraken2-braken2

### Alpha多样性

```
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
for tax in P F G S; do
Rscript $sd/otutab_rare.R \
      --input kraken2/bracken.${tax}.0.01-H \
      --depth 0 --seed 1 \
      --normalize kraken2/bracken.${tax}.norm \
      --output kraken2/bracken.${tax}.alpha;done

# 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
# Rscript $sd/alpha_boxplot.R -h # 查看参数
mkdir -p kraken2/${tax}
Rscript $sd/alpha_boxplot.R \
  -i kraken2/bracken.${tax}.alpha \
  -a shannon \
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
mv alpha_boxplot_TukeyHSD.txt  kraken2/${tax}/alpha;done

```

### Beta多样性

```
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

```

### 堆叠柱状图

```
# 以门(P)/种(S)水平为例，结果包括output.sample/group.pdf两个文件
tax=P
Rscript ${sd}/tax_stackplot.R \
  --input kraken2/bracken.${tax}.txt --design metadata.txt \
  --group Group --output kraken2/bracken.${tax}.stackplot \
  --legend 10 --width 89 --height 59
  
# 批量运行
for tax in D G P S; do
Rscript ${sd}/tax_stackplot.R \
      --input kraken2/bracken.${tax}.0.01-H --design metadata.txt \
      --group Group --output kraken2/bracken.${tax}.stackplot \
      --legend 10 --width 120 --height 70; done

```

### STAMP图

```
compare="Normal-Cancer"
tax=S
Rscript ${sd}/compare_stamp.R \
      --input kraken2/bracken.${tax}.0.01-H --metadata metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method "t.test" --pvalue 0.1 --fdr "none" \
      --width 100 --height 200 \
      --output kraken2/stamp_${compare}
```


## 基因组进化树注释table2itol

```
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

```

## MAGs分析

### MAGs与样本稀疏曲线

```
cd ${wd}
cd ..
Rscript ${sd}/MAG_sample_rare.R \
      --input binning/result/coverm/abundance.tsv \
      --taxonomy binning/result/coverm/taxonomy_MAG.txt \
      --output binning/result/coverm/
```


### MAGs不同功能数据库注释结果比较UpSet图

```
Rscript ${sd}/MAG_protiens_UpSet.R \
        --i result/eggnog/data_venn1.txt \
        --o result/eggnog/
        
## 单个功能数据库分组比较(COGs, KEGG, ECs, GOs), 这里以COGs注释结果作为示例
# 分组比较UpSet图
Rscript ${sd}/common_upset.R \
        --i result/eggnog/data_venn2.txt \
        --o result/eggnog/

# 提取交集并绘制交集组成饼图
Rscript ${sd}/upset_intersection_composition.R \
        --input result/eggnog/data_venn2.txt \
        --composition result/eggnog/COGs_data.txt \
        --output result/eggnog/

```

### MAGs功能通路桑基图

```
Rscript ${sd}/function_sankey.R \
        --i result/humann3/data_sankey.txt \
        --o result/humann3/
```


### MAGs完整性和污染率关系图

```
bash ${sd}/merged_tax.sh
sed 's/ /_/g' result/checkm2/taxonomy_merge.txt > result/checkm2/taxonomy_merge2.txt
Rscript ${sd}/MAG_quality_taxonomy.R \
        --i result/checkm2/taxonomy_merge2.txt \
        --o result/checkm2/
```

### MAGs系统发育分析树状图

```
# 肠道微生物样本示例
Rscript ${sd}/phylogenetic_tree.R \
        --input result/gtdb_95/tax.unrooted4.tree \
        --annotation result/gtdb_95/annotation4.txt \
        --output result/gtdb_95/

# 环境样本示例
Rscript ${sd}/phylogenetic_tree_env.R \
        --input result/gtdb_95/tax.unrooted2.tree \
        --annotation result/gtdb_95/annotation3.txt \
        --output result/gtdb_95/
        

Rscript ${sd}/phylogenetic_tree_env.R \
        --input result/gtdb_95/tax.unrooted4.tree \
        --annotation result/gtdb_95/annotation4.txt \
        --output result/gtdb_95/        

```

### MAGs耐药基因丰度热图

```
Rscript ${sd}/functional_gene_heatmap.R \
        --input binning/temp/ARG2/ARG_final.txt \
        --group binning/temp/ARG2/group_ARG.txt \
        --output binning/temp/

```


### 单菌基因组注释和绘图

# 注释：https://bakta.computational.bio/submit
# 注释和绘图：https://proksee.ca



### 宏基因组数据泛基因组分析绘图

# 绘图： https://anvi-server.org/



# 附录

## windows换行符处理

```
# 查看行尾是否有^M
cat -A $sd/metaphlan_hclust_heatmap.R | head
# 转换Windows为Linux换行符
dos2unix $sd/metaphlan_hclust_heatmap.R   

```

## 表格合并

    # 按末列注释，必须名叫KO
    Rscript $sd/mat_gene2ko.R  \
      -i temp/dbcan2/gene_fam.count \
      -o dbcan2/fam_merge.count

