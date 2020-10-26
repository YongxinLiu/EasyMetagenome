[TOC]

# 宏基因组分析流程 Pipeline of metagenomic analysis

    # 版本 Version：1.09, 2020/10/16
    # 系统要求 System: Linux Ubuntu 18.04+ / CentOS 7+
    # 依赖软件 Sofware: Rstudio server 1.1+、KneadData v0.7.4、MetaPhlAn v3.0.1、HUMAnN2 v3.0.0.alpha.3 ......
    # Windows/Mac访问Rstudio服务器，推荐使用Chrome浏览器，可选Edge，不要用IE因为存在兼容性问题

# 一、数据预处理 Data preprocessing

## 1.1 准备工作 Prepare
    
    # 首次使用请参照`soft_db.sh`脚本，安装软件和数据库(大约3-5天)
    # 学员U盘或服务器家目录(~)有项目文件夹meta, 
    # 1. 包含测序数据(seq/*.fq.gz)和和实验设计(result/metadata.txt);若无可手动从U盘上传。
    # 2. 使用谷歌浏览器访问RStudio服务器，具体IP地址课上通知，服务器可使用1个月
    # 3. Terminal中新建工作目录 mkdir -p meta ，右侧File中进入meta目录并上传流程pipeline.sh流程文件，单击打开

### 1.1.1 环境变量设置(每次开始分析前必须运行)

    # 设置数据库、软件和工作目录
    # 公共数据库database(db)位置，如多人共用可由管理员布置至/db，而个人可下载至~/db
    db=/db
    # Conda软件software安装目录，如/conda2或~/miniconda2
    soft=/conda2
    # 设置工作目录work directory(wd)，如meta
    wd=~/meta
    
    # 导入指定虚拟环境，如配置时叫metagenome_env
    export PATH=${soft}/envs/metagenome_env/bin/:$PATH
    source ${soft}/bin/activate metagenome_env
    # 创建并进入工作目录
    mkdir -p $wd 
    cd $wd
    # 指定某个R语言环境
    alias Rscript="/anaconda2/bin/Rscript --vanilla"

### 1.1.2 起始文件——序列和元数据


    # 创建3个常用子目录：序列，临时文件和结果
    mkdir -p seq temp result

    # 用户使用filezilla上传测序文件至seq目录，本次从其它位置复制，或从网络下载
    
    # 12个样本数据，只取了10万条PE100数据(20MB)作为演示，通常测序数据单样本>2千万条的PE150数据(6GB)

    # 上传实验设计(元数据) metadata.txt 于结果result目录
    # 检查文件格式，^I为制表符，$为Linux换行，^M$为Windows回车，^M为Mac换行符
    cat -A result/metadata.txt
    # 转换Windows回车为Linux换行
    sed -i 's/\r//' result/metadata.txt
    cat -A result/metadata.txt

    # 从其它目录复制测序数据
    # cp -rf /db/meta/seq/*.gz seq/
    
    # 从GSA数据库下载测序数据
    # 根据元数据GSA的CRA(批次)和CRR(样品)编号下载，每个双端测序样本两个文件，并按SampleID改名
    wget -c ftp://download.big.ac.cn/gsa/CRA002355/CRR117732/CRR117732_f1.fq.gz -O seq/C1_1.fq.gz
    wget -c ftp://download.big.ac.cn/gsa/CRA002355/CRR117732/CRR117732_r2.fq.gz -O seq/C1_2.fq.gz
    
    # 按实验设计批量下载并改名
    awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$6"/"$7"/"$7"_f1.fq.gz -O seq/"$1"_1.fq.gz")}' \
        <(tail -n+2 result/metadata.txt)
    awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$6"/"$7"/"$7"_r2.fq.gz -O seq/"$1"_2.fq.gz")}' \
        <(tail -n+2 result/metadata.txt)
    ls -lsh seq

### 1.1.3 了解工作目录和文件

    # 显示文件结构
    tree 
    # .
    # ├── pipeline.sh
    # ├── result
    # │   └── metadata.txt
    # ├── seq
    # │   ├── C1_1.fq.gz
    # │   ├── .._2.fq.gz
    # │   └── N6_2.fq.gz
    # ├── soft_db.sh
    # └── temp
    # pipeline.sh是分析流程代码；
    # seq目录中有12个样本双端测序，共24个序列文件；
    # temp是临时文件夹，存储分析中间文件，结束可全部删除节约空间
    # result是重要节点文件和整理化的分析结果图表，
    # 实验设计metadata.txt也在此
    
    # 查看压缩测序数据前8行
    zcat seq/C1_1.fq.gz | head -n8


## 1.2 (可选)FastQC质量评估Quality access

    # time统计运行时间，fastqc质量评估，
    # *.gz为原始数据，-t指定多线程，1m22s
    time fastqc seq/*.gz -t 1
    # 结果见seq目录，解读见 [数据的质量控制软件——fastQC](https://mp.weixin.qq.com/s/tDMih7ISLJcL4F4sWBq3Vw)
    
    # 生成多样品报告比较
    multiqc -d seq/ -o result/qc
    # 查看右侧result/qc目录中multiqc_report.html，View in Web Browser查看可交互式报告
    # 正常N和癌症C组GC含量组间存在差异
    

## 1.3 KneadData质控和去宿主

    # kneaddata是流程，它依赖trimmomatic质控和去接头，bowtie2比对宿主并筛选非宿主序列
    # 可只选一行中部分代码点击Run，如选中下行中#号后面命令查看程序帮助
    # kneaddata -h # 显示帮助
    # (可选) 序列预处理，解决NCBI SRA数据双端ID一致问题
    gunzip seq/*.gz
    sed -i '1~4 s/$/\\1/g' seq/*_1.fq
    sed -i '1~4 s/$/\\2/g' seq/*_2.fq
    pigz seq/*.fq
    
### 1.3.1 (可选)单样品质控

    # 请务必根据自己软件和数据库安装位置，conda env list 查看软件安装位置
    # 多行注释命令运行，可全选，按Ctrl+Shift+C进行注释的取消和添加
    # 10万条序列质控，时间10s~3m
    time kneaddata -i seq/C1_1.fq.gz -i seq/C1_2.fq.gz \
      -o temp/qc -v -t 3 --remove-intermediate-output \
      --trimmomatic /conda2/envs/metagenome_env/share/trimmomatic/ \
      --trimmomatic-options 'ILLUMINACLIP:/conda2/envs/metagenome_env/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' -db ${db}/kneaddata/human_genome/Homo_sapiens
    
### 1.3.1 多样品并行质控

    # 现实中是有一大堆样品，你可以逐个修改样品名运行，但并行队列管理是最优选择
    # 打will cite承诺引用并行软件parallel
    parallel --citation 
    
    # 每步分析产生多个文件时新建子文件夹
    mkdir -p temp/qc
    # 根据样本列表:::并行处理，并行j=2个任务，每个任务t=3个线程，2~7m
    time parallel -j 2 --xapply \
      "kneaddata -i seq/{1}_1.fq.gz \
      -i seq/{1}_2.fq.gz \
      -o temp/qc -v -t 3 --remove-intermediate-output \
      --trimmomatic /conda2/envs/metagenome_env/share/trimmomatic/ \
      --trimmomatic-options 'ILLUMINACLIP:/conda2/envs/metagenome_env/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' \
      -db ${db}/kneaddata/human_genome/Homo_sapiens" \
      ::: `tail -n+2 result/metadata.txt|cut -f1`

    # 质控结果汇总
    kneaddata_read_count_table \
      --input temp/qc \
      --output result/01kneaddata_sum.txt
    cat result/01kneaddata_sum.txt


## 1.4 (可选)质控后质量再评估 

    # 1-2m
    fastqc temp/qc/*_1_kneaddata_paired_* -t 2
    multiqc -d temp/qc/ -o result/qc/
    # 整理bowtie2, trimmomatic, fastqc报告，接头和PCR污染率一般小于1%



# 二、基于读长分析 Read-based (HUMAnN2)

中文教程：https://mp.weixin.qq.com/s/XkfT5MAo96KgyyVaN_Fl7g

英文教程：https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)


## 2.1 合并质控文件为HUMAnN2输入

    # HUMAnN2要求双端序列合并的文件作为输入
    mkdir -p temp/concat
    
    # for循环根据实验设计样本名批量双端序列合并
    # 注意星号和问号，分别代表多个和单个字符，当然大家更不能溜号~~~
    for i in `tail -n+2 result/metadata.txt | cut -f 1`;do 
      cat temp/qc/${i}*_1_kneaddata_paired_?.fastq \
      > temp/concat/${i}.fq; 
    done
    
    # 查看样品数量和大小
    ls -sh temp/concat/*.fq 


## 2.2 HUMAnN2计算物种和功能组成

    # 物种组成调用MetaPhlAn2, bowtie2比对至核酸序列；
    # 功能组成为humann2调用diamond比对至蛋白库11Gb
    # 输入文件：temp/concat/*.fq 每个样品质控后双端合并后的fastq序列
    # 输出文件：temp/humann2/ 目录下
    #           C1_pathabundance.tsv
    #           C1_pathcoverage.tsv
    #           C1_genefamilies.tsv
    # 整合后的输出：
    #           result/metaphlan2/taxonomy.tsv 物种丰度表
    #           result/metaphlan2/taxonomy.spf 物种丰度表（用于stamp分析）
    #           result/humann2/pathabundance_relab_stratified.tsv 通路丰度表
    #           result/humann2/pathabundance_relab_unstratified.tsv 通路丰度表
    #           stratified(每个菌的功能组成)和unstratified(功能组成)

    ### (可选)单样本运行，~1h
    time humann2 --input temp/concat/C3.fq  \
      --output temp/ --threads 3

    # 多样本并行计算，1-6h
    mkdir -p temp/humann2
    time parallel -j 2 \
      'humann2 --input {}  \
      --output temp/humann2/ ' \
      ::: temp/concat/*.fq > temp/log

## 2.3 物种组成表

    mkdir -p result/metaphlan2

### 2.3.1 样品结果合并

    merge_metaphlan_tables.py temp/humann2/*_humann2_temp/*_metaphlan_bugs_list.tsv | \
      sed 's/_metaphlan_bugs_list//g' > result/metaphlan2/taxonomy.tsv
    head -n3 result/metaphlan2/taxonomy.tsv

### 2.3.2 转换为stamp的spf格式

    metaphlan_to_stamp.pl result/metaphlan2/taxonomy.tsv \
      > result/metaphlan2/taxonomy.spf
    head -n3 result/metaphlan2/taxonomy.spf
    # 下载metadata.txt和taxonomy.spf使用stamp分析

### 2.3.3 (可选) Python绘制热图

    # c设置颜色方案，top设置物种数量，minv最小相对丰度，s标准化方法，log为取10为底对数，xy为势图宽和高，图片可选pdf/png/svg格式
    time metaphlan_hclust_heatmap.py \
      --in result/metaphlan2/taxonomy.tsv \
      --out result/metaphlan2/heatmap.pdf \
      -c jet --top 15 --minv 0.1 \
      -s log -x 0.4 -y 0.2
    # 帮助文档见 metaphlan_hclust_heatmap.py -h
    # (可选)Excel筛选数据ImageGP绘图

### 2.3.4 (可选) R绘制热图

    # 推荐在Windows中运行
    # 若在服务器上，请保持工作目录不变，并修改R脚本目录位置

    # 本地分析目录
    cd /c/meta
    # 调置脚本所有目录，服务器位于 /db/script/
    R=/c/meta/db/script
    
    # 显示脚本帮助 help
    Rscript ${R}/metaphlan_hclust_heatmap.R -h
    # 按指定列合并、排序并取Top25种绘制热图
    # -i输入MetaPhlAn3文件；-t 分类级别，可选Kingdom/Phylum/Class/Order/Family/Genus/Species，界门纲目科属种株，推荐门，目，属
    # -n 输出物种数量，默认为25，最大为合并后的数量
    # -o输出图表前缀，默认根据输入文件、物种级别和数量自动生成；
    Rscript db/script/metaphlan_hclust_heatmap.R \
      -i result/metaphlan3/taxonomy.spf \
      -t Family -n 25 \
      -w 183 -e 118 \
      -o result/metaphlan3/heatmap_Family

    # 属水平的Top30
    Rscript db/script/metaphlan_hclust_heatmap.R \
      -i result/metaphlan3/taxonomy.spf \
      -t Genus \
      -n 30 \
      -o result/metaphlan3/heatmap_Genus
    
## 2.4 功能组成分析

    mkdir -p result/humann3

### 2.4.1 功能组成合并、标准化和分层

    # 合并通路丰度(pathabundance)，含功能和物种组成
    humann_join_tables \
      --input temp/humann3 \
      --file_name pathabundance \
      --output result/humann3/pathabundance.tsv
    # 可选基因家族(genefamilies 太多)，通路覆盖度(pathcoverage)层面
    # 删除列名多余信息
    sed -i 's/_Abundance//g' result/humann3/pathabundance.tsv
    # 预览文件头尾格式
    head -n3 result/humann3/pathabundance.tsv
    tail -n3 result/humann3/pathabundance.tsv

    # 标准化为相对丰度(1)relab或百万比cpm
    humann_renorm_table \
      --input result/humann3/pathabundance.tsv \
      --units relab \
      --output result/humann3/pathabundance_relab.tsv
    head -n5 result/humann3/pathabundance_relab.tsv

    # 分层结果：功能和对应物种表(stratified)和功能组成表(unstratified)
    humann_split_stratified_table \
      --input result/humann3/pathabundance_relab.tsv \
      --output result/humann3/ 
    # 可以使用stamp进行统计分析

### 2.4.2 添加分组和差异比较

    # 在通路丰度中添加分组
    ## 提取样品列表
    head -n1 result/humann3/pathabundance.tsv | sed 's/# Pathway/SampleID/' | tr '\t' '\n' > temp/header
    ## 对应分组
    awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' result/metadata.txt temp/header | tr '\n' '\t'|sed 's/\t$/\n/' > temp/group
    # 合成样本、分组+数据
    cat <(head -n1 result/humann3/pathabundance.tsv) temp/group <(tail -n+2 result/humann3/pathabundance.tsv) > result/humann3/pathabundance.pcl

    # 组间比较，样本量少无差异
    humann_associate --input result/humann3/pathabundance.pcl  \
        --focal-metadatum Group --focal-type categorical \
        --last-metadatum Group --fdr 0.05 \
        --output result/humann3/associate.txt

### 2.4.3 通路物种组成柱状图

    # barplot展示腺苷核苷酸合成的物种组成
    humann_barplot \
        --input result/humann3/pathabundance.pcl \
        --focal-feature PWY-7219 \
        --focal-metadatum Group --last-metadatum Group \
        --output result/humann3/barplot1.pdf
    # --sort sum 按丰度排序
    humann_barplot --sort sum \
        --input result/humann3/pathabundance.pcl \
        --focal-feature PWY-7219 \
        --focal-metadatum Group --last-metadatum Group \
        --output result/humann3/barplot2.pdf
    # --sort sum metadata 按丰度排序
    humann_barplot --sort sum metadata \
        --input result/humann3/pathabundance.pcl \
        --focal-feature PWY-7219 \
        --focal-metadatum Group --last-metadatum Group \
        --output result/humann3/barplot3.pdf


## 2.5 GraPhlAn图

    # metaphlan2 to graphlan
    export2graphlan.py --skip_rows 1,2 -i result/metaphlan2/taxonomy.tsv \
      --tree temp/merged_abundance.tree.txt \
      --annotation temp/merged_abundance.annot.txt \
      --most_abundant 1000 --abundance_threshold 20 --least_biomarkers 10 \
      --annotations 3,4 --external_annotations 7
    # 参数说明见PPT，或运行 export2graphlan.py --help
    # graphlan annotation
    graphlan_annotate.py --annot temp/merged_abundance.annot.txt \
      temp/merged_abundance.tree.txt  temp/merged_abundance.xml
    # output PDF figure, annoat and legend
    graphlan.py temp/merged_abundance.xml result/metaphlan2/graphlan.pdf \
      --external_legends 


## 2.6 LEfSe差异分析和Cladogram

    # 准备输入文件，修改样本品为组名，**非通用代码**，可手动修改
    # 只保留组名C, N, 删除数字重复编号，去除注释行
    head -n3 result/metaphlan2/taxonomy.tsv
    sed '1 s/[0-9]*//g' result/metaphlan2/taxonomy.tsv \
      | grep -v '#' > result/metaphlan2/lefse.txt
    head -n3 result/metaphlan2/lefse.txt
    
    # LEfSe本地代码供参考，可选Xshell下运行或ImageGP在线分析
    # 格式转换为lefse内部格式
    lefse-format_input.py \
      result/metaphlan2/lefse.txt \
      temp/input.in -c 1 -o 1000000
    # 运行lefse
    run_lefse.py temp/input.in \
      temp/input.res
    
    # 绘制物种树注释差异
    lefse-plot_cladogram.py temp/input.res \
      result/metaphlan2/lefse_cladogram.pdf --format pdf
    
    # 绘制所有差异features柱状图
    lefse-plot_res.py temp/input.res \
      result/metaphlan2/lefse_res.pdf --format pdf
        
    # 绘制单个features柱状图
    # 查看显著差异features，按丰度排序
    grep -v '-' temp/input.res | sort -k3,3n 
    # 手动选择指定feature绘图，如Firmicutes
    lefse-plot_features.py -f one \
      --feature_name "k__Bacteria.p__Firmicutes" \
      --format pdf \
      temp/input.in temp/input.res \
      result/metaphlan2/lefse_Firmicutes.pdf
    
    # 批量绘制所有差异features柱状图
    lefse-plot_features.py -f diff \
      --archive none --format pdf \
      temp/input.in temp/input.res \
      result/metaphlan2/lefse_


## 2.7 kraken2物种注释reads

    # 还可以进行contig、gene、bin层面的序列物种注释

### 2.7.1 物种注释

    mkdir -p temp/kraken2
    
    # 单样本注释，5m
    time kraken2 --db ${db}/kraken2 --paired temp/qc/C1_1_kneaddata_paired*.fastq \
      --threads 3 --use-names --use-mpa-style --report-zero-counts \
      --report temp/kraken2/C1_report \
      --output temp/kraken2/C1_output
    
    # 多样本并行，12个样本2个同时各3线程，
    time parallel -j 2 \
      "kraken2 --db ${db}/kraken2 --paired temp/qc/{1}_1_kneaddata_paired*.fastq \
      --threads 3 --use-names --use-mpa-style --report-zero-counts \
      --report temp/kraken2/{1}_report \
      --output temp/kraken2/{1}_output" \
      ::: `tail -n+2 result/metadata.txt | cut -f 1`
    # 屏幕会输出各样品注释比例，和运行时间 20~30m

### 2.7.2 汇总样品物种组成表

    mkdir -p result/kraken2
    # 输出结果行数相同，但不一定顺序一致，要重新排序
    parallel -j 1 \
      'sort temp/kraken2/{1}_report | cut -f 2 | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count ' \
      ::: `tail -n+2 result/metadata.txt | cut -f 1`
    # 提取第一样本品行名为表行名
    header=`tail -n 1 result/metadata.txt | cut -f 1`
    echo $header
    sort temp/kraken2/${header}_report | cut -f 1 | \
      sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
    head -n3 temp/kraken2/0header_count
    # paste合并样本为表格
    ls temp/kraken2/*count
    paste temp/kraken2/*count > result/kraken2/taxonomy_count.txt

### 2.7.3 物种多样性分析

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # Rscript ${db}/script/kraken2alpha.R -h # 查看帮助和参数详细
    # 输入count文件，按最小值抽平为norm文件，计算alpha多样性
    Rscript ${db}/script/kraken2alpha.R \
      --input result/kraken2/taxonomy_count.txt \
      --depth 0 \
      --normalize result/kraken2/tax_norm.txt \
      --output result/kraken2/tax_alpha.txt
    
    # 绘制Alpha多样性箱线图
    # 可选richness/chao1/ACE/shannon/simpson
    Rscript ${db}/script/alpha_boxplot.R \
      -i result/kraken2/tax_alpha.txt \
      -t shannon \
      -d result/metadata.txt \
      -n Group \
      -o result/kraken2/alpha_shannon \
      -w 4 -e 2.5
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript ${db}/script/alpha_boxplot.R -i result/kraken2/tax_alpha.txt -t ${i} \
      -d result/metadata.txt -n Group -w 4 -e 2.5 \
      -o result/kraken2/alpha_${i}
    done

### 2.7.4 物种组成

    # 转换为metaphalan2 spf格式，但ncbi注释不完整
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["d"]="unclassified";a["p"]="unclassified";a["c"]="unclassified";a["o"]="unclassified";a["f"]="unclassified";a["g"]="unclassified";a["s"]="unclassified";a["S"]="unclassified"; \
      split($1,x,"|");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print a["d"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"],a["S"],$0;}' \
      result/kraken2/tax_norm.txt > temp.txt
    cut -f 1-8,10- temp.txt > result/kraken2/tax_norm.spf
    sed -i '1 s/unclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain/' \
      result/kraken2/tax_norm.spf

    # 绘制热图
    Rscript db/script/metaphlan_hclust_heatmap.R \
      -i result/kraken2/tax_norm.spf \
      -t Genus \
      -n 25 \
      -w 183 -e 118 \
      -o result/kraken2/heatmap_Genus

    # 绘制属水平Top30箱线图
    Rscript db/script/metaphlan_boxplot.R \
      -i result/kraken2/tax_norm.spf \
      -t Genus \
      -n 30 \
      -o result/kraken2/boxplot_Genus
    # 绘制门水平Top10箱线图
    Rscript db/script/metaphlan_boxplot.R \
      -i result/kraken2/tax_norm.spf \
      -t Phylum \
      -n 10 -w 6 -e 4 \
      -o result/kraken2/boxplot_Phylum
  

# 三、组装分析流程 Assemble-based

## 3.1 拼接 Assembly

### 3.1.1 MEGAHIT拼接

    # 删除旧文件夹，否则megahit无法运行
    rm -rf temp/megahit
    # 组装，6~30m，TB级数据需几天至几周
    time megahit -t 6 \
        -1 `tail -n+2 result/metadata.txt|cut -f 1|sed 's/^/temp\/qc\//;s/$/_1_kneaddata_paired_1.fastq/'| tr '\n' ','|sed 's/,$//'` \
        -2 `tail -n+2 result/metadata.txt|cut -f 1|sed 's/^/temp\/qc\//;s/$/_1_kneaddata_paired_2.fastq/'| tr '\n' ','|sed 's/,$//'` \
        -o temp/megahit 
    # 检查点：查看拼接结果，8.2M，通常300M~5G
    ls -sh temp/megahit/final.contigs.fa
    # 预览重叠最前6行，前60列字符
    head -n6 temp/megahit/final.contigs.fa | cut -c1-60

### 3.1.2 (可选) metaSPAdes精细拼接

    # 精细但使用内存和时间更多，15~65m
    time metaspades.py -t 6 -m 100 \
      `tail -n+2 result/metadata.txt|cut -f 1|sed 's/^/temp\/qc\//;s/$/_1_kneaddata_paired_1.fastq/'|sed 's/^/-1 /'| tr '\n' ' '` \
      `tail -n+2 result/metadata.txt|cut -f 1|sed 's/^/temp\/qc\//;s/$/_1_kneaddata_paired_2.fastq/'|sed 's/^/-2 /'| tr '\n' ' '` \
      -o temp/metaspades
    # 23M，contigs体积更大
    ls -sh temp/metaspades/contigs.fasta

### 3.1.3 QUAST评估

    mkdir -p result/megahit/
    ln temp/megahit/final.contigs.fa result/megahit/
    quast.py result/megahit/final.contigs.fa -o result/megahit/quast -t 2
    # 生成report文本tsv/txt、网页html、PDF等格式报告
    
    # (可选) megahit和metaspades比较
    time quast.py --label "megahit,metapasdes" \
        temp/megahit/final.contigs.fa \
        temp/metaspades/contigs.fasta \
        -o temp/quast
    
    # (可选)metaquast评估，更全面，但需下载相关数据库，受网速影响可能时间很长
    # metaquast based on silva, and top 50 species genome to access
    time metaquast.py result/megahit/final.contigs.fa -o result/megahit/metaquast

## 3.2 基因预测、去冗余和定量 Gene prediction, cluster & quantitfy

### 3.2.1 metaProdigal基因预测

    mkdir -p temp/prodigal
    # prodigal的meta模式预测基因，7s，>和2>&1记录分析过程至gene.log
    time prodigal -i  result/megahit/final.contigs.fa \
        -d temp/prodigal/gene.fa \
        -o temp/prodigal/gene.gff \
        -p meta -f gff > temp/prodigal/gene.log 2>&1 
    # 查看日志是否运行完成，有无错误
    tail temp/prodigal/gene.log
    # 统计基因数量19185
    grep -c '>' temp/prodigal/gene.fa 

### 3.2.2 基因聚类/去冗余cd-hit

    mkdir -p result/NR
    # aS覆盖度，c相似度，G局部比对，g最优解，T多线程，M内存0不限制
    # 2万基因2m，2千万需要2000h，可多线程加速
    time cd-hit-est -i temp/prodigal/gene.fa \
        -o result/NR/nucleotide.fa \
        -aS 0.9 -c 0.95 -G 0 -g 0 -T 0 -M 0
    # 统计非冗余基因数量19146，单次拼接结果数量下降不大，多批拼接冗余度高
    grep -c '>' temp/prodigal/nucleotide.fa
    # 翻译核酸为对应蛋白序列，emboss
    transeq -sequence result/NR/nucleotide.fa -outseq result/NR/protein.fa -trim Y 
    # 序列名自动添加了_1，为与核酸对应要去除
    sed -i 's/_1 / /' result/NR/protein.fa

### 3.2.3 基因定量salmon

    mkdir -p temp/salmon
    # 直接运行salmon找不到lib库，可使用程序完整路径解决问题(error while loading shared libraries: liblzma.so.0)
    alias salmon="/conda2/envs/metagenome_env/share/salmon/bin/salmon"
    # 建索引, -t序列, -i 索引，10s
    time salmon index -t result/NR/nucleotide.fa -p 9 \
        -i temp/salmon/index 
    # 定量，l文库类型自动选择，p线程，--meta宏基因组模式, 2个任务并行12个样，共42s
    time parallel -j 2 \
      "/conda2/envs/metagenome_env/share/salmon/bin/salmon quant \
        -i temp/salmon/index -l A -p 3 --meta \
        -1 temp/qc/{1}_1_kneaddata_paired_1.fastq \
        -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
        -o temp/salmon/{1}.quant" \
        ::: `tail -n+2 result/metadata.txt | cut -f 1`
    # 合并
    mkdir -p result/salmon
    salmon quantmerge \
        --quants temp/salmon/*.quant \
        -o result/salmon/gene.TPM
    salmon quantmerge \
        --quants temp/salmon/*.quant \
        --column NumReads -o result/salmon/gene.count
    sed -i '1 s/.quant//g' result/salmon/gene.*
    # 预览结果表格
    head -n3 result/salmon/gene.*

## 3.3 功能基因注释

    # 输入数据：上一步预测的蛋白序列 result/NR/protein.fa
    # 中间结果：temp/eggnog/protein.emapper.seed_orthologs
    #           temp/eggnog/output.emapper.annotations
    #           temp/eggnog/output
    
    # COG定量表：result/eggnog/cogtab.count
    #            result/eggnog/cogtab.count.spf (用于STAMP)
    
    # KO定量表：result/eggnog/kotab.count
    #           result/eggnog/kotab.count.spf  (用于STAMP)
    
    # CAZy碳水化合物注释和定量：result/dbcan2/cazytab.count
    #                           result/dbcan2/cazytab.count.spf (用于STAMP)
    
    # 抗生素抗性：result/resfam/resfam.count
    #             result/resfam/resfam.count.spf (用于STAMP)
    
    # 这部分可以拓展到其它数据库

### 3.3.1 基因注释eggNOG(COG/KEGG/CAZy)

    # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2
    
    # 进入Python3环境，emapper2必在python3中运行
    source ${soft}/bin/activate humann3
    emapper.py --version
    
    # diamond比对基因至eggNOG数据库, 1~5h
    mkdir -p temp/eggnog
    time emapper.py -m diamond --no_annot --no_file_comments \
      --data_dir ${db}/eggnog5 --cpu 3 -i result/NR/protein.fa \
      -o temp/eggnog/protein --override

    # 比对结果功能注释, 1h
    time emapper.py --annotate_hits_table \
      temp/eggnog/protein.emapper.seed_orthologs --no_file_comments \
    	-o temp/eggnog/output --cpu 3 --data_dir ${db}/eggnog --override

    # 结果注释表头, 重点1序列名，9KO，16CAZy，21COG分类，22注释
    mkdir -p result/eggnog
    sed '1 i Name\tortholog\tevalue\tscore\ttaxonomic\tprotein\tGO\tEC\tKO\tPathway\tModule\tReaction\trclass\tBRITE\tTC\tCAZy\tBiGG\ttax_scope\tOG\tbestOG\tCOG\tdescription\t' \
      temp/eggnog/output.emapper.annotations > temp/eggnog/output
    head -n3 temp/eggnog/output


Python3脚本整理COG/KO/CAZy数据表

    summarizeAbundance.py \
      -i result/salmon/gene.TPM \
      -m temp/eggnog/output \
      -c '9,16,21' -s ',+,+*' -n raw \
      -o result/eggnog/eggnog
    # eggnog.CAZy.raw.txt  eggnog.COG.raw.txt  eggnog.KO.raw.txt
    # STAMP的spf格式，结合metadata.tsv进行COG/KO/Description差异比较
    # KO
    sed -i 's/^ko://' result/eggnog/eggnog.KO.raw.txt
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      ${db}/eggnog/KO.anno \
      result/eggnog/eggnog.KO.raw.txt| \
      sed 's/^\t/Description\t/' > result/eggnog/eggnog.KO.TPM.spf
    # CAZy
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
       ${db}/dbcan2/fam_description.txt result/eggnog/eggnog.CAZy.raw.txt | \
      sed 's/^\t/Description\t/' > result/eggnog/eggnog.CAZy.TPM.spf
    # COG
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
      ${db}/eggnog/COG.anno \
      result/eggnog/eggnog.COG.raw.txt | sed '1 s/^/Level1\tLevel2/'> \
      result/eggnog/eggnog.COG.TPM.spf

### 3.3.2 碳水化合物dbCAN2(可选)

    # 比对CAZy数据库, 用时10m; 加--sensitive更全但更慢1h
    mkdir -p temp/dbcan2
    time diamond blastp --db ${db}/dbCAN2/CAZyDB.07312018 --query result/NR/protein.fa \
    	--outfmt 6 --threads 2 --max-target-seqs 1 --quiet -e 1e-5 \
    	--out temp/dbcan2/gene_diamond.f6 
    # 整理比对数据为表格
    mkdir -p result/dbcan2
    # 提取基因对应基因家族，同一基因存在1对多，只取第一个
    cut -f 1,2 temp/dbcan2/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | \
    	sed 's/_[0-9]*$//' |sed '1 i Name\tKO' > temp/dbcan2/gene_fam.list
    # 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/dbcan2/gene_fam.list \
      result/salmon/gene.count | sed '/\t$/d' > temp/dbcan2/gene_fam.count
    # 按基因家族合并
    Rscript ${db}/script/mat_gene2ko.R -i temp/dbcan2/gene_fam.count -o result/dbcan2/cazytab
    # 结果中添加FAM注释，spf格式用于stamp分析
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' ${db}/dbCAN2/fam_description.txt \
    	result/dbcan2/cazytab.count > result/dbcan2/cazytab.count.spf

### 3.3.3 抗生素抗性ResFam

    mkdir -p temp/resfam result/resfam
    # 比对至抗生素数据库 1m
    time diamond blastp --db ${db}/resfam/Resfams-proteins --query result/NR/protein.fa \
    	--outfmt 6 --threads 2 --max-target-seqs 1 --quiet -e 1e-5 --sensitive \
    	--out temp/resfam/gene_diamond.f6
    # 提取基因对应基因家族
    cut -f 1,2 temp/resfam/gene_diamond.f6 | uniq | \
      sed '1 i Name\tResGeneID' > temp/resfam/gene_fam.list
    # 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      temp/resfam/gene_fam.list result/salmon/gene.count | \
    	sed '/^\t/d' > result/resfam/resfam.count
    # 统计注释基因的比例, 488/19147=2.5%
    wc -l result/salmon/gene.count result/resfam/resfam.count 
    # 结果中添加FAM注释，spf格式用于stamp分析
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4"\t"$3"\t"$2} NR>FNR{print a[$1],$0}' \
      ${db}/resfam/Resfams-proteins_class.tsv  result/resfam/resfam.count \
      > result/resfam/resfam.count.spf



# 四、挖掘单菌基因组/分箱(Binning)

## 4.1 MetaWRAP

    # 主要使用MetaWRAP，演示基于官方测试数据
    # 主页：https://github.com/bxlab/metaWRAP
    # 挖掘单菌基因组，需要研究对象复杂度越低、测序深度越大，结果质量越好。要求单样本6GB+，复杂样本如土壤推荐数据量30GB+，至少3个样本
    # 上面的演示数据12个样仅140MB，无法获得单菌基因组，这里使用官方测序数据演示讲解
    # 软件和数据库布置需2-3天，演示数据分析过程超10h，标准30G样也需3-30天，由服务器性能决定。

### 4.1.1 准备数据和环境变量

    # 准备原始数据从头分析，详见公众号或官网教程
    # 这里我们从质控后数据和拼接结果开始
    cd ${wd}
    mkdir -p binning && cd binning
    mkdir -p temp
    # 这里基于质控clean数据和拼接好的contigs，自己链接自上游分析
    # 7G质控数据，输入数据文件名格式必须为*_1.fastq和*_2.fastq
    mkdir -p seq
    cd seq
    # 方法1. 下载测序数据
    for i in `seq 7 9`;do
        wget -c http://210.75.224.110/share/meta/metawrap/ERR01134${i}_1.fastq.gz
        wget -c http://210.75.224.110/share/meta/metawrap/ERR01134${i}_2.fastq.gz
    done
    # gunzip *.gz # 解压文件
    # rename .fq .fastq *.fq # 批量修改扩展名
    # 方法2. 复制准备好的数据
    ln -s ${db}/metawrap/*.fastq ./
    cd ..
    # megahit拼接结果
    mkdir -p megahit
    cd megahit
    # wget -c http://210.75.224.110/share/meta/metawrap/final.contigs.fa.gz
    # gunzip *.gz
    ln -s ${db}/metawrap/*.fa ./
    cd ../..
    
    # 加载运行环境
    cd ${wd}/binning
    conda activate metawrap


### 4.1.2 运行三种分箱软件

    # 输入文件为contig和clean reads
    # 调用三大主流binning程序cococt, maxbin2, metabat2
    # 8p线程2h，24p耗时1h
    # nohup 和 & 保证任务在后台不被中断，且记录输出内容到 nohup.out(可选)
    nohup metawrap binning -o temp/binning -t 2 -a temp/megahit/final.contigs.fa \
      --metabat2 --maxbin2 --concoct temp/qc/ERR*.fastq &
    # 用自己的文件，替换输出文件名为 *1_kneaddata_paired*.fastq 
    # 输出文件夹 temp/binning 包括3种软件结果和中间文件

### 4.1.3 Bin提纯

    # 8线程2h， 24p 1h
    cd ${wd}/binning
    # rm -rf temp/bin_refinement
    metawrap bin_refinement \
      -o temp/bin_refinement \
      -A temp/binning/metabat2_bins/ \
      -B temp/binning/maxbin2_bins/ \
      -C temp/binning/concoct_bins/ \
      -c 50 -x 10 -t 2
    # 查看高质量Bin的数量，10个，见temp/bin_refinement/metawrap_50_10_bins.stats目录
    wc -l temp/bin_refinement/metawrap_50_10_bins.stats
    # 结果改进程度见temp/bin_refinement/figures/目录


### 4.1.4 Bin定量

    # 使用salmon计算每个bin在样本中相对丰度
    # 耗时3m，系统用时10m，此处可设置线程，但salmon仍调用全部资源
    # 需要指定输出文件夹，包括4.3中的参数的输出目录
    metawrap quant_bins -b temp/bin_refinement/metawrap_50_10_bins -t 8 \
      -o temp/bin_quant -a temp/megahit/final.contigs.fa temp/qc/ERR*.fastq
    # 文件名字改变
    # 结果包括bin丰度热图`temp/bin_quant/bin_abundance_heatmap.png`
    # 如果想自己画图，原始数据位于`temp/bin_quant/bin_abundance_table.tab`
    ls -l temp/bin_quant/bin_abundance_heatmap.png

### 4.1.5 Bin注释

    # Taxator-tk对每条contig物种注释，再估计bin整体的物种，11m (用时66 min)
    metawrap classify_bins -b temp/bin_refinement/metawrap_50_10_bins \
      -o temp/bin_classify -t 2
    # 注释结果见`temp/bin_classify/bin_taxonomy.tab`
    
    # 基于prokka基因注释，4m
    metaWRAP annotate_bins -o temp/bin_annotate \
      -b temp/bin_refinement/metawrap_50_10_bins  -t 2
    # 每个bin基因注释的gff文件bin_funct_annotations, 
    # 核酸ffn文件bin_untranslated_genes，
    # 蛋白faa文件bin_translated_genes
    

## (可选)MetaWRAP单样本分别组装和分箱

多样本受硬件、计算时间限制无法完成时，需要单样本组装、分析。或想进一步提高组装质量，减少污染和杂合度，也可以单样本组装。

### 参数设定

    # 样本名
    i=ERR011347
    # 线程数
    p=12
    # 任务数
    j=2
    # 定义完整度和污染率的阈值(50, 5; Finn NBT 2020; 50, 10, Bowers NBT 2017)
    c=50
    x=10
    
输和文件在seq目录

    mkdir -p seq
    ln -s `pwd`/temp/qc/*.fastq seq/
    
 ### 1 megahit组装

单样本并行组装，13m，314m

    rm -rf temp/megahit_*
    time parallel -j ${j} \
    "metawrap assembly \
        -1 seq/{}_1.fastq \
        -2 seq/{}_2.fastq \
        -o temp/megahit_{} \
        -m 100 -t ${p} --megahit" \
     ::: `ls seq/|cut -f1 -d '_'|uniq`  

### 2 运行三种bin软件

    # 192p, 15m (concoct会使用所有线程)
    parallel -j ${j} \
    "metawrap binning \
        -o temp/binning_{} -t ${p} \
        -a temp/megahit_{}/final_assembly.fasta \
        --metabat2 --maxbin2 --concoct \
        seq/{}_*.fastq" \
    ::: `ls seq/|cut -f1 -d '_'|uniq`
     
### 3 Bin提纯

    # 24p，10h
    parallel -j ${j} \
    "metawrap bin_refinement \
      -o temp/bin_refinement_{} -t ${p} \
      -A temp/binning_{}/metabat2_bins/ \
      -B temp/binning_{}/maxbin2_bins/ \
      -C temp/binning_{}/concoct_bins/ \
      -c ${c} -x ${x}" \
    ::: `ls seq/|cut -f1 -d '_'|uniq`

## 4.2 dRep去冗余种/株基因组集

    # 进入虚拟环境，没用用conda安装
    conda activate drep

合并所有bin至同一目录

    mkdir -p temp/drep_in
    # 混合组装分箱并重命名
    ln -s `pwd`/temp/bin_refinement/metawrap_50_10_bins/bin.* temp/drep_in/
    rename 's/bin/mix_all/' temp/drep_in/bin.*
    # 单样品组装分箱结果重命名
    for i in `ls seq/|cut -f1 -d '_'|uniq`;do
        ln -s `pwd`/temp/bin_refinement_${i}/metawrap_50_10_bins/bin.* temp/drep_in/
        rename "s/bin./s_${i}./" temp/drep_in/bin.*
    done
    # 统计混合和单样本来源数据，10个混，5个单
    ls temp/drep_in/|cut -f 1 -d '_'|uniq -c
    # 统计混合批次/单样本来源
    ls temp/drep_in/|cut -f 2 -d '_'|cut -f 1 -d '.' |uniq -c

按种水平去冗余：15个为10个，8个来自混拼，2个来自单拼

    mkdir -p temp/drep95
    # 15个，10min
    dRep dereplicate temp/drep95/ \
      -g temp/drep_in/*.fa \
      -sa 0.95 -nc 0.30 -comp 50 -con 10 -p 24
    
主要结果：

- 非冗余基因组集：dereplicated_genomes/*.fa
- 聚类信息表：data_tables/Cdb.csv
- 聚类和质量图：firgures/*clustering*

(可选)按株水平汇总

    mkdir -p temp/drep99
    dRep dereplicate temp/drep99/ \
      -g temp/drep_in/*.fa \
      -sa 0.99 -nc 0.30 -comp 50 -con 10 -p 24

## 4.3 GTDB-tk物种注释和进化树

启动软件所在虚拟环境

    conda activate gtdbtk

细菌基因组物种注释

以上面鉴定的10个种为例，注意扩展名要与输入文件一致，可使用压缩格式gz。主要结果文件描述：此9个细菌基因组，结果位于tax.bac120开头的文件，如物种注释 tax.bac120.summary.tsv。古菌结果位于tax.ar122开头的文件中。

    mkdir -p temp/gtdb_classify
    # 10个基因组，24p，26min
    gtdbtk classify_wf \
        --genome_dir temp/drep95/dereplicated_genomes \
        --out_dir temp/gtdb_classify \
        --extension fa \
        --prefix tax \
        --cpus 2


多序列对齐结果建树

    # 以9个细菌基因组的120个单拷贝基因建树，1s
    mkdir -p temp/gtdb_infer
    gtdbtk infer \
        --msa_file temp/gtdb_classify/tax.bac120.user_msa.fasta \
        --out_dir temp/gtdb_infer \
        --prefix tax \
        --cpus 2

树文件可使用iTOL在线美化，也可使用GraphLan本地美化。

## 4.4 table2itol制作树注释文件

以gtdb-tk物种注释(tax.bac120.summary.tsv)和drep基因组评估(Widb.csv)信息为注释信息

    mkdir -p result/itol
    # 制作分类学表
    tail -n+2 temp/gtdb_classify/tax.bac120.summary.tsv|cut -f 1-2|sed 's/;/\t/g'|sed '1 s/^/ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' > result/itol/tax.txt
    # 基因组评估信息
    sed 's/,/\t/g;s/.fa//' temp/drep95/data_tables/Widb.csv|cut -f 1-7,11|sed '1 s/genome/ID/' > result/itol/genome.txt
    # 整合注释文件
    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result/itol/genome.txt result/itol/tax.txt|cut -f 1-8,10- > result/itol/annotation.txt

table2itol制作注释文件

    cd result/itol/
    # 设置脚本位置
    db=~/meta/script/table2itol
    
    ## 方案1. 分类彩带、数值热图、种标签
    # -a 找不到输入列将终止运行（默认不执行）-c 将整数列转换为factor或具有小数点的数字，-t 偏离提示标签时转换ID列，-w 颜色带，区域宽度等， -D输出目录，-i OTU列名，-l 种标签替换ID
    Rscript ${db}/table2itol.R -a -c double -D plan1 -i ID -l Species -t %s -w 0.5 annotation.txt
    # 生成注释文件中每列为单独一个文件

    ## 方案2. 数值柱形图，树门背景色，属标签
    Rscript ${db}/table2itol.R -a -d -c none -D plan2 -b Phylum -i ID -l Genus -t %s -w 0.5 annotation.txt

    ## 方案3.分类彩带、整数为柱、小数为热图
    Rscript ${db}/table2itol.R -c keep -D plan3 -i ID -t %s annotation.txt

    ## 方案4. 将整数转化成因子生成注释文件
    Rscript ${db}/table2itol.R -a -c factor -D plan4 -i ID -l Genus -t %s -w 0 annotation.txt

    
# 附录1. 测试版humann3

教程：https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(humann3)

    # 进入Python3环境，humann3必在python3中运行
    source ${soft}/bin/activate humann3

## 2.1 合并质控文件为humann3输入

    # humann3也要求双端序列合并的文件作为输入
    mkdir -p temp/concat
    
    # for循环根据实验设计样本名批量双端序列合并
    # 注意星号和问号，分别代表多个和单个字符，当然大家更不能溜号~~~
    for i in `tail -n+2 result/metadata.txt | cut -f 1`;do 
      cat temp/qc/${i}*_1_kneaddata_paired_?.fastq \
      > temp/concat/${i}.fq; 
    done
    
    # 查看样品数量和大小
    ls -sh temp/concat/*.fq 


## 2.2 humann3计算物种和功能组成

    注意：humann3的命令为humann，以区别于humann2

    # 物种组成调用MetaPhlAn3, bowtie2比对至核酸序列；
    # 功能组成为humann3调用diamond比对至蛋白库11Gb
    # 输入文件：temp/concat/*.fq 每个样品质控后双端合并后的fastq序列
    # 输出文件：temp/humann3/ 目录下
    #           C1_pathabundance.tsv
    #           C1_pathcoverage.tsv
    #           C1_genefamilies.tsv
    # 整合后的输出：
    #           result/metaphlan2/taxonomy.tsv 物种丰度表
    #           result/metaphlan2/taxonomy.spf 物种丰度表（用于stamp分析）
    #           result/humann3/pathabundance_relab_stratified.tsv 通路丰度表
    #           result/humann3/pathabundance_relab_unstratified.tsv 通路丰度表
    #           stratified(每个菌的功能组成)和unstratified(功能组成)

    # (可选)单样本运行，~1h
    time humann --input temp/concat/C2.fq  \
      --output temp/ --threads 16

    # 多样本并行计算，1-6h
    mkdir -p temp/humann3
    time parallel -j 2 \
      'humann --input {}  \
      --output temp/humann3/ ' \
      ::: temp/concat/*.fq > temp/log

## 2.3 物种组成表

    mkdir -p result/metaphlan3

    # 样品结果合并
    merge_metaphlan_tables.py temp/humann3/*_humann_temp/*_metaphlan_bugs_list.tsv | \
      sed 's/_metaphlan_bugs_list//g' > result/metaphlan3/taxonomy.tsv
    head -n3 result/metaphlan3/taxonomy.tsv
    # 结果比metaphlan2多了一行注释，一列NCBI_tax_id，调整为标准表格
    sed '/^#/d' result/metaphlan3/taxonomy.tsv | cut -f 1,3- > result/metaphlan3/taxonomy.txt

    # 转换为stamp的spf格式
    metaphlan_to_stamp.pl result/metaphlan3/taxonomy.txt \
      > result/metaphlan3/taxonomy.spf
    head -n3 result/metaphlan3/taxonomy.spf
    # 下载metadata.txt和taxonomy.spf使用stamp分析

## 2.4 功能组成分析

    mkdir -p result/humann3

### 2.4.1 功能组成合并、标准化和分层

    # 合并通路丰度(pathabundance)，含功能和物种组成
    humann_join_tables \
      --input temp/humann3 \
      --file_name pathabundance \
      --output result/humann3/pathabundance.tsv
    # 可选基因家族(genefamilies 太多)，通路覆盖度(pathcoverage)层面
    # 删除列名多余信息
    sed -i 's/_Abundance//g' result/humann3/pathabundance.tsv
    # 预览文件头尾格式
    head -n3 result/humann3/pathabundance.tsv
    tail -n3 result/humann3/pathabundance.tsv

    # 标准化为相对丰度(1)relab或百万比cpm
    humann_renorm_table \
      --input result/humann3/pathabundance.tsv \
      --units relab \
      --output result/humann3/pathabundance_relab.tsv
    head -n5 result/humann3/pathabundance_relab.tsv

    # 分层结果：功能和对应物种表(stratified)和功能组成表(unstratified)
    humann_split_stratified_table \
      --input result/humann3/pathabundance_relab.tsv \
      --output result/humann3/ 
    # 可以使用stamp进行统计分析

### 2.4.2 添加分组和差异比较

    # 在通路丰度中添加分组
    ## 提取样品列表
    head -n1 result/humann3/pathabundance.tsv | sed 's/# Pathway/SampleID/' | tr '\t' '\n' > temp/header
    ## 对应分组
    awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' result/metadata.txt temp/header | tr '\n' '\t'|sed 's/\t$/\n/' > temp/group
    # 合成样本、分组+数据
    cat <(head -n1 result/humann3/pathabundance.tsv) temp/group <(tail -n+2 result/humann3/pathabundance.tsv) > result/humann3/pathabundance.pcl

    # 组间比较，样本量少无差异
    humann_associate --input result/humann3/pathabundance.pcl  \
        --focal-metadatum Group --focal-type categorical \
        --last-metadatum Group --fdr 0.05 \
        --output result/humann3/associate.txt

### 2.4.3 通路物种组成柱状图

    # barplot展示腺苷核苷酸合成的物种组成
    humann_barplot \
        --input result/humann3/pathabundance.pcl \
        --focal-feature PWY-7219 \
        --focal-metadatum Group --last-metadatum Group \
        --output result/humann3/barplot1.pdf
    # --sort sum 按丰度排序
    humann_barplot --sort sum \
        --input result/humann3/pathabundance.pcl \
        --focal-feature PWY-7219 \
        --focal-metadatum Group --last-metadatum Group \
        --output result/humann3/barplot2.pdf
    # --sort sum metadata 按丰度排序
    humann_barplot --sort sum metadata \
        --input result/humann3/pathabundance.pcl \
        --focal-feature PWY-7219 \
        --focal-metadatum Group --last-metadatum Group \
        --output result/humann3/barplot3.pdf

# 附录2. emapper+eggNOG4.5(COG/KEGG)

    # diamond比对基因至eggNOG数据库, 1~5h
    mkdir -p temp/eggnog
    time emapper.py -m diamond --no_annot --no_file_comments \
      --data_dir ${db}/eggnog --cpu 6 -i result/NR/protein.fa \
      -o temp/eggnog/protein --override

    # 比对结果功能注释, 1h
    time emapper.py --annotate_hits_table \
      temp/eggnog/protein.emapper.seed_orthologs --no_file_comments \
    	-o temp/eggnog/output --cpu 2 --data_dir ${db}/eggnog --override

    # 结果注释表头, 重点1序列名，7KO，12COG分类，13注释
    mkdir -p result/eggnog
    sed '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' \
      temp/eggnog/output.emapper.annotations > temp/eggnog/output

    # 1.整理COG表
    # 提取12列COG分类
    cut -f 1,12 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' \
      >temp/eggnog/1cog.list
    # 基因丰度矩阵末尾添加对应cog编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' \
      temp/eggnog/1cog.list result/salmon/gene.count | \
    	sed '/\t$/d' | sed '1 s/COG/KO/' > temp/eggnog/gene_cog.count
    # 按COG类型合并表格，输出count和RPM值，n设置标准化单位，默认1M，可选100/1
    Rscript ${db}/script/mat_gene2ko.R -i temp/eggnog/gene_cog.count \
      -o result/eggnog/cogtab -n 1000000
    # STAMP的spf格式，结果metadata.tsv进行KO或Description差异比较
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
      ${db}/eggnog/COG.anno result/eggnog/cogtab.count > \
      result/eggnog/cogtab.count.spf

    # 2. 整理KO表
    # 提取基因KO表，基因1对多个KO时只提取第一个KO
    cut -f 1,7 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' \
      > temp/eggnog/2ko.list
    # 基因丰度矩阵末尾添加对应KO编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' \
      temp/eggnog/2ko.list \
      result/salmon/gene.count | sed '/\t$/d' > temp/eggnog/gene_ko.count
    # 合并基因表为KO表，输出count值和tpm值
    Rscript ${db}/script/mat_gene2ko.R -i temp/eggnog/gene_ko.count \
      -o result/eggnog/kotab -n 1000000
    # STAMP的spf格式，结果metadata.txt进行KO或Description差异比较
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      ${db}/eggnog/KO.anno result/eggnog/kotab.count | \
      sed 's/^\t/Undescription\t/' > result/eggnog/kotab.count.spf
      
# 版本更新记录

## 1.08 2020.7.20

1. KneadData提供数据预处理双端标签唯一命令，兼容最新版；
2. 提供HUMAnN3测试版的安装和分析流程(附录1)；
3. eggNOG升级为emapper 2.0和eggNOG 5.0流程，结果列表从13列变为22列，新增CAZy注释。emapper 1.0版本见附录2。

## 1.08 2020.10.16

1. 新增二、三代混合组装OPERA-MS软件使用 (31Megahit)
2. 新增eggNOG-mapper结果COG/KO/CAZy整理脚本summarizeAbundance.py，删除旧版Shell+R代码 (32Annotation)
3. 新增MetaWRAP单样本分箱流程 (33Binning)
4. 新增dRep实现基因组去冗余 (34Genomes)
5. 新增GTDB-Tk基因组物种注释和进化树构建 (34Genomes)



