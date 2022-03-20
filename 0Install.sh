[TOC]

# 易宏基因组软件和数据库 EasyMetagenome software & database

    # 版本: 1.14, 2022/3/25
    # 测试环境为Linux Ubuntu 20.04 / CentOS 7.7

## 安装前准备：软件和数据库位置

    # 数据库安装位置，默认~/db目录(无需管理权限)，管理员可安装至/db，方便大家使用
    db=~/db
    mkdir -p ${db} && cd ${db}
    # 软件安装位置，也可能为miniconda3
    soft=~/miniconda3

### EasyMetagenome依赖软件和数据库(db)

EasyMetagenome依赖流程，包括很多脚本、常用小软件和数据库的合集，网址：https://github.com/YongxinLiu/EasyMicrobiome
    
    # 三种下载数据库的方法：任选其一即可
    
    # 方法1. 可使用wget或在网址 https://github.com/YongxinLiu/EasyMicrobiome 中Code Download ZIP下载压缩包，并解压
    
    # 方法2. http备用链接下载
    wget -c http://210.75.224.110/db/EasyMicrobiome.zip
    unzip EasyMicrobiome.zip

    # 方法3. git下载，需要安装git并配置好帐号
    git clone https://github.com/YongxinLiu/EasyMicrobiome
    # 可选旧版更新 cd EasyMicrobiome && git pull

    # 添加linux命令可执行权限
    chmod +x EasyMicrobiome/linux/*
    
    # 添加环境变量
    echo 'export PATH="$PATH:~/db/EasyMicrobiome/linux"' >> ~/.bashrc

### 软件管理器Conda

    # 下载最新版miniconda3，~49M
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # 安装，-b批量，-f无提示，-p目录，许可协议打yes
    bash Miniconda3-latest-Linux-x86_64.sh -b -f
    # 激活，然后关闭终端重开，提示符前出现(base)即成功
    ~/miniconda3/condabin/conda init
    # 查看版本，conda 4.11.0, python 3.9.7
    conda -V 
    python --version
    # 添加常用频道
    conda config --add channels conda-forge
    # 添加生物学软件频道，http://bioconda.github.io/ 查询软件
    conda config --add channels bioconda
    # conda默认配置文件为 ~/.condarc 查看配置文件位置
    conda config --show-sources
    
中文详细安装教程参考：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)

(可选)如果找不到conda，可手动临时添加conda至环境变量。可以添加至~/.bashrc文件中永久环境变量，需将${soft}替换为你的安装目录，如

    # export PATH="${soft}/bin:$PATH"

查看虚拟环境列表 

    conda env list

创建虚拟环境示例，防污染环境变量示例(如果软件安装卡在Solving environment步骤，可尝试新建环境)

    conda create -n meta
    # 加载环境
    conda activate meta
    # 退出环境
    conda deactivate

## conda环境打包pack和安装(推荐)

conda安装经常会卡在Collecting package metadata或Solving environment。我们推荐一种直接下载解压的安装方式，可以加速环境的部署。

下载安装包并解压安装：

    # 指定环境名称，如 meta, humann2, kraken2, eggnog, rgi, metawrap1.3, drep, gtdbtk1.5, humann3, qiime2-2021.2
    soft=~/miniconda3
    n=metawrap1.3
    # 下载
    wget -c http://210.75.224.110/db/conda/${n}.tar.gz
    # 指定安装目录
    mkdir -p ${soft}/envs/${n}
    tar -xvzf ${n}.tar.gz -C ${soft}/envs/${n}
    # 方法1. 启动环境
    conda activate $n
    # 初始化环境
    conda unpack
    # 退出环境
    conda deactivate
    # 方法2. 绝对目录激活环境
    # source ${soft}/envs/${n}/bin/activate

(可选)安装好的环境下打包导出，以宏基因组主环境meta为例

    # 安装conda-pack实现打包
    conda install conda-pack -c conda-forge
    # conda 安装不成功可用pip安装
    # pip conda-pack
    
    # conda环境包统一存放
    cd ~/db/conda/
    # 设置环境名，如meta, humann2, kraken2, eggnog, rgi, metawrap1.3, drep, gtdbtk1.5
    n=meta
    conda pack -n ${n} -o ${n}.tar.gz
    # 导出软件安装列表
    conda activate ${n}
    conda env export > ${n}.yml
    # 添加权限，方便下载和别人使用
    chmod 755 *
    

## 常用数据库快速下载(示例)

### 基因组-kneaddata去宿主

 人类基因组，为Kneaddata作者自定义构建好的索引，可直接下载使用。下面是我的备份链接，官方安装方法见下方“标准化安装”段落
 
    mkdir -p ${db}/kneaddata/human_genome && cd  ${db}/kneaddata/human_genome
    wget -c http://210.75.224.110/db/kneaddata/human_genome/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz
    tar xvzf Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz

自定义基因组构建索引，大多数基因组可在ensembl genome下载。此处以拟南芥为例，访问 http://plants.ensembl.org/index.html ，选择Arabidopsis thaliana —— Download DNA sequence (FASTA)，选择toplevel右键复制链接，填入下面链接处

    # 创建子目录
    mkdir -p ${db}/kneaddata/ath && cd  ${db}/kneaddata/ath
    # 下载
    wget -c http://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    # 解压
    gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    # 简化文件名
    mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa tair10.fa
    # bowtiew建索引，输入文件，输出文件前缀，9线程2分
    time bowtie2-build -f tair10.fa tair10 --threads 9 --seed 1

### humann2数据库

    mkdir -p ${db}/humann2 && cd ${db}/humann2 
    # 泛基因组下载、指定目录解压
    wget -c http://210.75.224.110/db/humann2/full_chocophlan_plus_viral.v0.1.1.tar.gz
    mkdir -p chocophlan
    tar xvzf full_chocophlan_plus_viral.v0.1.1.tar.gz -C chocophlan/
    # 功能注释
    wget -c http://210.75.224.110/db/humann2/full_mapping_1_1.tar.gz
    mkdir -p utility_mapping
    tar xvzf full_mapping_1_1.tar.gz -C utility_mapping/
    # 功能基因蛋白序列
    wget -c http://210.75.224.110/db/humann2/uniref90_annotated_1_1.tar.gz
    mkdir -p uniref
    tar xvzf uniref90_annotated_1_1.tar.gz -C uniref/
    # 物种标记基因数据库
    wget -c http://210.75.224.110/db/humann2/metaphlan2.tar.gz
    mkdir -p metaphlan2
    tar xvzf metaphlan2.tar.gz -C metaphlan2/
    
    humann2_config --update database_folders utility_mapping ${db}/humann2/utility_mapping
    humann2_config --update database_folders nucleotide ${db}/humann2/chocophlan
    humann2_config --update database_folders protein ${db}/humann2/uniref
    # metaphlan2数据库默认位于~/miniconda3/envs/biobakery/bin/ 程序所在目录的db_v20 (v0.11)和databases(2.8.1)下
    humann2_config --print
    # 例如，链接下载的metaphlan2：
    ln -s `pwd`/metaphlan2 ~/miniconda3/envs/biobakery/bin/databases

# 标准化安装

以上移植方法无法使用时，或没有的软件，可参考以下代码。软件更新频率，有问题可先查看软件的官方最新安装教程是否更新，以最新版为准

## 1质控软件

### 质量评估fastqc

    # =为指定版本，-c指定安装源，均可加速安装
    # -y为同意安装
    conda install fastqc=0.11.9 -c bioconda -y
    fastqc -v

### 评估报告汇总multiqc

    # 注1.7为Python2环境，1.8/9新版本需要Python3的环境
    conda install multiqc=1.9 -c bioconda -y 
    multiqc --version

### 质量控制流程kneaddata

    conda install kneaddata=0.7.4 -c bioconda -y 
    kneaddata --version
    trimmomatic -version # 0.39
    bowtie2 --version # 2.4.2

    # 查看可用数据库
    kneaddata_database
    # 包括人基因组bowtie2/bmtagger、人类转录组、核糖体RNA和小鼠基因组
    # 下载人基因组bowtie2索引 3.44 GB
    mkdir -p ${db}/kneaddata/human_genome
    kneaddata_database --download human_genome bowtie2 ${db}/kneaddata/human_genome

    
## 2有参分析流程MetaPhlAn2、HUMAnN2

### HUMAnN2安装

方法1. 新建虚拟环境安装：安装MetaPhlAn2、HUMAnN2和所有依赖关系，并记录主要软件版本

    # conda install humann2=2.8.1 -c bioconda -y
    conda create -n humann2 humann2=2.8.1 -c bioconda -y
    conda activate humann2
    # 记录核心软件版本
    humann2 --version # humann2 v2.8.1
    metaphlan2.py -v # MetaPhlAn version 2.7.5 (6 February 2018)
    diamond help #  v0.8.22.84

方法2. Conda导入和导出环境

    # 安装conda-pack
    conda install -c conda-forge conda-pack
    #conda 安装不成功可用pip安装
    # pip conda-pack
    # 安装好的环境下打包导出
    conda pack -n humann2 -o humann2.tar.gz
    # 下载
    wget -c http://210.75.224.110/db/humann2/humann2.tar.gz

    # 新建文件夹存放humann2环境
    mkdir -p ~/miniconda3/envs/humann2
    tar -xzf humann2.tar.gz -C ~/miniconda3/envs/humann2
    # 激活环境
    conda activate humann2
    # 初始化
    conda unpack

测试流程是否可用

    humann2_test

### HUMAnN2物种和功能数据库

显示可用分类、泛基因组和功能数据库

    humann2_databases

安装数据库(注：数据库下载慢或失败，附录有国内备份链接)

    cd ${db}
    mkdir -p ${db}/humann2 # 建立下载目录
    # 输助比对数据库 593MB
    humann2_databases --download utility_mapping full ${db}/humann2
    # 微生物泛基因组 5.37 GB
    humann2_databases --download chocophlan full ${db}/humann2
    # 功能基因diamond索引 10.3 GB
    humann2_databases --download uniref uniref90_diamond ${db}/humann2
    
    # 设置数据库位置
    # 显示参数
    humann2_config --print
    # 如修改线程数，推荐3-8，根据实际情况调整
    humann2_config --update run_modes threads 3
    humann2_config --update database_folders utility_mapping ${db}/humann2/utility_mapping
    humann2_config --update database_folders nucleotide ${db}/humann2/chocophlan
    humann2_config --update database_folders protein ${db}/humann2/uniref
    # metaphlan2数据库默认位于~/miniconda3/envs/biobakery/bin/ 程序所在目录的db_v20 (v0.11)和databases(2.8.1)下
    humann2_config --print
    # 例如，链接下载的metaphlan2：
    ln -s `pwd`/metaphlan2 ~/miniconda3/envs/biobakery/bin/databases

### Microbiome helper

主页：https://github.com/mlangill/microbiome_helper

下载并安装

    # 下载、解压 、添加环境变量
    wget -c https://github.com/LangilleLab/microbiome_helper/archive/master.zip
    unzip master.zip
    export PATH=`pwd`/microbiome_helper-master:$PATH
    # 写入bashrc永久添加环境
    echo "export PATH=`pwd`/microbiome_helper-master:\$PATH" >> ~/.bashrc
    
    # metaphlan_to_stamp.pl 这个脚本有修改，修改后的脚本在db/script/下，也在QQ群文件中
    #----修改的内容如下----------------------------------------
    # my @taxa_ranks=("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Strain");
    > 
    # #start with assumption that this is not metaphlan2
    # my $metaphlan2_version=1;
    #------------------------------------------------------
    
### 物种组成美化GraPhlAn

    # GraPhlAn核心程序包
    conda install graphlan=1.1.3 -c bioconda -y
    graphlan.py --version # 1.1.3 (5 June 2018)
    # GraPhlAn输入文件制作程序，如转换LEfSe、Metaphlan2结果格式为GraPhlAn用于绘图
    conda install export2graphlan # 38 KB
    export2graphlan.py -h # 0.22 of 05 May 2020

### 生物标记鉴定和可视化LEfSe

    conda install lefse -c bioconda -y # 76.3 MB, 1.0.8.post1
    
    # Rstudio中运行命令调用R版本问题的解决
    # 在Rstudio中默认调用Rstudio的R，具体写在/etc/rstudio/rserver.conf
    # 或在R中用Sys.getenv()["R_HOME"]
    # 在rpy2中print(robjects.r)可以查看其调用的r版本
    
    # 指定lefse调用的R版本，需根据conda实际目录修改
    sed -i "2 i os.environ['R_HOME'] = '/conda/envs/meta/lib/R/'" \
      /conda/envs/meta/share/lefse-1.0.8.post1-1/lefse.py

### 物种注释Kraken2

物种注释：基于LCA算法的物种注释kraken2  https://ccb.jhu.edu/software/kraken/

- 软件安装

安装方法1. 直接安装并查看版本

    conda install kraken2 -c bioconda -y
    kraken2 --version # 2.1.1
    conda install bracken=2.6.0 -c bioconda

安装方法2. 新建环境安装并启动软件环境

    conda create -n kraken2 -y -c bioconda kraken2
    conda activate kraken2
    conda install bracken=2.6.0 -c bioconda
    # krakentools 0.1 补充脚本
    conda install krakentools -c bioconda
    # krona绘图
    conda install krona -c bioconda

安装方法3. conda导入和导出

    # 从安装好的环境打包，此处为Ubuntu 20.04LTS
    # conda pack -n kraken2 -o ~/db/kraken2/kraken2.tar.gz
    # 添加权限，否则别人无法下载
    # chmod 777 ~/db/kraken2/kraken2.tar.gz
    # 下载压缩包
    wget -c http://210.75.224.110/db/kraken2/kraken2.tar.gz
    # 新建文件夹存放kraken2环境
    mkdir -p ~/miniconda3/envs/kraken2
    # 解压环境到指定目录
    tar -xzf kraken2.tar.gz -C ~/miniconda3/envs/kraken2
    # 激活环境和初始化
    conda activate kraken2
    conda unpack
    # 启动指定目录中的环境
    # source ~/miniconda3/envs/kraken2/bin/activate

- 数据库

下载数据库(NCBI每2周更新一次)，记录下载日期和大小。需根据服务器内存、使用目的选择合适方案。--standard标准模式下只下载5种数据库：古菌archaea、细菌bacteria、人类human、载体UniVec_Core、病毒viral。也可选直接下载作者构建的索引，还包括bracken的索引。

- 方法1. 数据库下载

下载标准+原生动物+真菌+植物8GB(PlusPFP-8)数据库，包括kraken2和bracken2的索引。更多版本数据库详见：https://benlangmead.github.io/aws-indexes/k2 。

    mkdir -p ~/db/kraken2 && cd ~/db/kraken2

方案1. 迷你库(8G，低配推荐)

    # 压缩包5.2G，
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_8gb_20210517.tar.gz
    # 备用地址：
    # wget -c http://210.75.224.110/db/kraken2/k2_pluspf_8gb_20210517.tar.gz
    tar xvzf k2_pluspf_8gb_20210517.tar.gz

方案2. 完整库(70G，高配推荐)

压缩包70G，解压后100G。指定解压目录，包括时间和类型。201202为更新时间，pfp指标准库+原生动物+真菌+植物。注：但我使用中发现仍然没有真菌。

    d=201202pfp
    mkdir ${d}
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20201202.tar.gz
    tar xvzf k2_pluspfp_20201202.tar.gz -C ${d}

方案3. 第3方个性数据库(人和病毒含新冠)，仅用于病毒检测

https://genexa.ch/sars2-bioinformatics-resources/

    wget -c https://storage.googleapis.com/sars-cov-2/kraken2_h%2Bv_20200319.tar.gz
    mkdir -p 200319hv/
    tar xvzf kraken2_h+v_20200319.tar.gz -C 200319hv/


- (可选)方法2. 数据库安装

方案1. 标准库安装，下载数据~100GB，时间由网速决定，索引5h，多线程可加速至1h完成
    
    cd ${db}
    d=210808
    mkdir -p kraken2/$d && cd kraken2/$d
    kraken2-build --standard --threads 24 --db ./
    
方案2. 自定义微生物数据库，如标准+真菌+原生动物+质粒+植物

    cd ${db}
    d=210808fpf
    mkdir -p kraken2/$d && cd kraken2/$d
    # 显示帮助
    kraken2-build -h
    # 下载物种注释
    kraken2-build --download-taxonomy --threads 24 --db ./
    # 下载数据库，需要12-24小时
    for i in archaea bacteria UniVec_Core viral human fungi plasmid protozoa plant; do
        kraken2-build --download-library $i --threads 24 --db ./
    done
    # 确定的库建索引，4p,4h
    time kraken2-build --build --threads 48 --db ./
    # bracken索引，长度推荐100/150, 24p,1h;
    time bracken-build -d ./ -t 24 -k 35 -l 100
    time bracken-build -d ./ -t 24 -k 35 -l 150


## 3基因组拼接、注释和定量

### 拼接megahit/metaSPAdes

    # 快速组装megahit
    conda install megahit
    megahit -v # v1.2.9
    
    # metaSPAdes拼接，只是spades系列中的一个定制脚本
    conda install spades # 12.8 MB
    metaspades.py -v # v3.14.0 [metaSPAdes mode]
    
### 组装评估QUEST

    conda install quast -y # 87.2 MB
    metaquast.py -v # QUAST v5.0.2 (MetaQUAST mode)
    
### 基因预测prokka

    # 细菌基因组注释prokka
    conda install prokka -y # 352.8 MB
    prokka -v # 1.14.6
    
### 非冗余基因集cd-hit

    conda install cd-hit # 790 KB
    cd-hit -v # 4.8.1 (built on May 14 2019)
    
    # emboss transeq工具，93.9 MB
    conda install emboss -y
    embossversion # 6.6.0.0
    
### 定量工具salmon

    conda install salmon -y # 17.8 MB
    salmon -v # 1.3.0


## 4基因功能注释

### KEGG层级注释

https://www.kegg.jp/kegg-bin/show_brite?ko00001.keg 下载htext

    # 转换ABCD为列表
    kegg_ko00001_htext2tsv.pl -i ko00001.keg -o ko00001.tsv
    # 统计行数，2021.1月版55761行，整理后为55103个条目
    wc -l ko00001.*
    # 统计各级数量, /54/527/23917
    for i in `seq 1 2 8`;do
        cut -f ${i} ko00001.tsv|sort|uniq|wc -l ; done
    # 生成KO编号和注释列表
    cut -f 7,8 ko00001.tsv|sort|uniq|sed '1 i KO\tDescription' \
      > KO_description.txt
    # KO与通路(Pathway)对应表，用于合并D级为C级
    awk 'BEGIN{FS=OFS="\t"} {print $7,$6}' ko00001.tsv | sed '1 i KO\tpathway' \
      > KO_path.list
      
      
### 蛋白同源综合注释eggNOG

eggNOG http://eggnogdb.embl.de

    # 新版要求python3.7，需要新环境
    conda create -n eggnog
    conda activate eggnog
    # 安装eggnog比对工具emapper
    conda install eggnog-mapper -y -c bioconda
    emapper.py --version # 2.1.6
    
    # 下载常用数据库，注意设置下载位置
    mkdir -p ${db}/eggnog && cd ${db}/eggnog
    # -y默认同意，-f强制下载，eggnog.db.gz 7.9G+4.9G
    download_eggnog_data.py -y -f --data_dir ./
    
    # 下载方式2(可选)：链接直接下载
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz # 6.3G
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz # 69M
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz # 4.9G
    gunzip *.gz
    
    # 备用链接
    wget -c http://210.75.224.110/db/eggnog/5.0.2/eggnog.db.gz
    wget -c http://210.75.224.110/db/eggnog/5.0.2/eggnog.taxa.tar.gz
    wget -c http://210.75.224.110/db/eggnog/5.0.2/eggnog_proteins.dmnd.gz
    gunzip *.gz
    tar xvf eggnog.taxa.tar

    # 如果内存够大，复制eggNOG至内存加速比对
    # cp eggnog.* /dev/shm

### 碳水化合物CAZy

    # dbCAN2 http://bcb.unl.edu/dbCAN2
    # 创建数据库存放目录并进入
    mkdir -p ${db}/dbcan2 && cd ${db}/dbcan2
    # 下载序列和描述
    wget -c https://bcb.unl.edu/dbCAN2/download/CAZyDB.09242021.fa
    wget -c https://bcb.unl.edu/dbCAN2/download/Databases/V10/CAZyDB.07292021.fam-activities.txt
    # 提取基因家簇对应注释
    # wget http://210.75.224.110/db/dbcan2/CAZyDB.07292021.fam-activities.txt
    grep -v '#' CAZyDB.07292021.fam-activities.txt \
      |sed 's/  //'| \
      sed '1 i CAZy\tDescription' \
      > CAZy_description.txt
    # 备用数据库下载地址并解压 
    # wget -c http://210.75.224.110/db/dbcan2/CAZyDB.09242021.fa.gz
    # gunzip CAZyDB.09242021.fa.gz
    # diamond建索引，1G，1-18m
    conda activate eggnog
    diamond --version # 2.0.13
    time diamond makedb --in CAZyDB.09242021.fa \
      --db CAZyDB.09242021
    # 压缩原始数据节约空间
    gzip CAZyDB.09242021.fa


### 抗生素抗性基因CARD

官网：https://card.mcmaster.ca

Github: https://github.com/arpcard/rgi

Bioconda: http://bioconda.github.io/recipes/rgi/README.html

软件安装

    # 方法1. Conda新建环境安装rgi(推荐)
    conda create -n rgi -c bioconda rgi=5.2.1

    # 方法2. rgi环境导出和导入
    # 安装好的环境下打包导出，262M
    conda pack -n rgi -o rgi.tar.gz
    
    # 下载软件包
    wget -c http://210.75.224.110/db/conda/rgi.tar.gz
    # 新建文件夹存放rgi环境
    mkdir -p ~/miniconda3/envs/rgi
    tar -xzf rgi.tar.gz -C ~/miniconda3/envs/rgi
    # 激活环境
    conda activate rgi
    conda unpack
    # source ~/miniconda3/envs/rgi/bin/activate
    
    # 方法4. docker安装(需要权限)，版本号见conda页面
    docker pull quay.io/biocontainers/rgi:5.2.1--pyhdfd78af_1


数据库部署

    # 下载最新版数据库，2.9M (2021-7-7, 3.1.2)
    wget -c https://card.mcmaster.ca/latest/data
    # 解压后20M
    tar -xvf data ./card.json
    # 加载数据库
    rgi load --card_json card.json
    
    # (可选)宏基因组分析扩展数据库和加载
    rgi card_annotation -i card.json
    rgi load -i card.json --card_annotation card_database_v3.1.2.fasta

    
### 抗生素抗性基因Resfam

    # http://dantaslab.wustl.edu/resfams
    mkdir -p ${db}/resfam && cd ${db}/resfam
    # 官网的数据格式非常混乱, 推荐下载我手动整理的索引和注释
    wget http://210.75.224.110/share/Resfams-proteins.dmnd # 1.5 MB
    wget http://210.75.224.110/share/Resfams-proteins_class.tsv # 304 KB

### 毒力因子数据库VFDB

    # http://www.mgc.ac.cn/VFs/ 
    mkdir -p ${db}/vfdb && cd ${db}/vfdb
    # 毒力因子描述文件
    wget -c http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
    # 核心数据库(966K)
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
    # 完整数据库(5.3M)
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
    # 解压
    gunzip *.gz


## 5分箱工具

### metawrap分箱和提纯

    # 方法1. conda在线安装物种注释和分箱流程 https://github.com/bxlab/metaWRAP
    conda create -y -n metawrap python=2.7 # 22.2MB
    conda activate metawrap
    conda config --add channels ursky
    conda install -y -c ursky metawrap-mg # 1.14 GB, v1.2

    # 方法2. 安装好的环境下打包导出，最新版1.3
    # 设置conda位置
    soft=~/miniconda3/
    # conda pack -n metawrap1.3 -o metawrap1.3.tar.gz
    # 下载最新环境
    wget -c http://210.75.224.110/db/metawrap/metawrap1.3.tar.gz
    # 新建文件夹存放metawrap1.3环境
    mkdir -p ${soft}/envs/metawrap1.3
    tar -xzf metawrap1.3.tar.gz -C ${soft}/envs/metawrap1.3
    # 激活环境
    source ${soft}/envs/metawrap1.3/bin/activate    
    
    
    # 相关数据库，大小近300GB
    # 这里我们安装数据库到`~/db`目录，保证你有权限，
    # 但要保证至少有500GB的空间。请根据你的情况修改为自己有权限且空间足够的位置。
    # 多人使用，建议管理员统一安装节省空间
    cd ${db}

    
    ## CheckM用于Bin完整和污染估计和物种注释
    mkdir -p checkm && cd checkm
    # 下载文件275 MB，解压后1.4 GB
    wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    # 国内备用链接
    # wget -c http://210.75.224.110/share/meta/checkm/checkm_data_2015_01_16.tar.gz
    tar -xvf *.tar.gz
    # rm *.gz
    # 设置数据库位置
    checkm data setRoot
    # 按提示输出你数据下载的路径或直接回车默认为当前位置
    
    ## NCBI_nt核酸序列用于bin物种注释
    # 41GB，我下载大约12h；解压后99GB
    cd ${db}
    mkdir -p NCBI_nt && cd NCBI_nt
    wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
    # 备用下载链接，或百度云下载
    # wget -c http://210.75.224.110/share/meta/NCBI_nt/filelist.txt
    # for a in `cat filelist.txt`; do wget -c http://210.75.224.110/share/meta/NCBI_nt/$a; done
    for a in nt.*.tar.gz; do tar xzf $a; done &

    ## NCBI物种信息
    # 压缩文件45M，解压后351M
    cd ${db}
    mkdir NCBI_tax
    cd NCBI_tax
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xvf taxdump.tar.gz
    
    ## 人类基因组去宿主
    cd ${db}
    mkdir -p metawrap/BMTAGGER && cd metawrap/BMTAGGER
    wget -c ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
    gunzip *fa.gz
    cat *fa > hg38.fa
    rm chr*.fa
    # 上方下载太慢，使用国内备份链接手动下载
    wget -c http://210.75.224.110/share/meta/metawrap/BMTAGGER/hg38.fa
    bmtool -d hg38.fa -o hg38.bitmask
    srprism mkindex -i hg38.fa -o hg38.srprism -M 100000
    
    ## KRAKEN物种注释数据库
    # 下载建索引需要 > 300GB以上空间，完成后占用192GB空间
    cd ${db}
    mkdir -p kraken
    kraken-build --standard --threads 24 --db kraken > log &
    kraken-build --db kraken --clean
    # 手动下载
    cd kraken
    wget -c http://210.75.224.110/share/meta/kraken/database.kdb
    wget -c http://210.75.224.110/share/meta/kraken/database.idx
    mkdir -p taxonomy && cd taxonomy
    wget -c http://210.75.224.110/share/meta/kraken/taxonomy/nodes.dmp
    wget -c http://210.75.224.110/share/meta/kraken/taxonomy/names.dmp
    # 从其它位置复制
    # cp -r /db/kraken/* ./

    ## 数据库位置设置
    which config-metawrap
    # 配置文件通常为~/miniconda3/envs/metawrap/bin/config-metawrap
    # 使用Rstudio/vim等文本编辑器来修改数据库的位置

### drep基因组去冗余

挑单菌测序的基因组存在大量冗余。

metawrap混合分箱的结果中冗余度非常低，甚至无冗余。而单样本、分批次分箱的结果中存在大量冗余，需要采用derep分箱获得非冗余的基因组。

    conda create -n drep
    conda activate drep
    conda install drep -c bioconda


### GTDB细菌基因组注释和进化分析

Github: https://github.com/Ecogenomics/GTDBTk

GTDB-Tk是一个软件工具包，用于根据基因组数据库分类法GTDB为细菌和古细菌基因组分配客观的分类法。它旨在与最近的进展一起使用，从而可以直接对环境样本中获得数百或数千个由基因组组装的基因组（MAG）进行物种分类注释。它也可以用于分离和单细胞的基因组物种注释。

本次测试版本为 GTDB-Tk v1.3.0，于2020年7月17发布，参考数据为95版。

硬件要求：

内存100Gb

硬盘27Gb

64核1小时可分析1000个细菌基因组

Conda安装：

    conda create -n gtdbtk
    conda activate gtdbtk
    # gtdbtk-1.3.0, 2020-9-27
    conda install -c bioconda gtdbtk

download-db.sh自动下载数据库，将下载至conda中的envs/gtdbtk/share/gtdbtk-1.3.0/db/：

    download-db.sh
    
(可选)手动下载和配置GTDB参考基因组最新版(测试时为95版，34Gb)

    mkdir -p ~/db/gtdb & cd ~/db/gtdb
    # 下载解压
    wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
    tar zxvf gtdbtk_data.tar.gz
    # 设置数据库位置，注意修改软件安装位置
    locate gtdbtk.sh # 查找配置文件位置
    # 修改PATH=后面的路径为数据库解压目录，如/home/meta/db/gtdb/release95/
    vim /conda/envs/gtdbtk/etc/conda/activate.d/gtdbtk.sh

此外：GTDB数据库国内下载备份链接见 http://nmdc.cn/datadownload


## 6 其它

    # Bin可视化VizBin (可选)
    sudo apt-get install libatlas3-base libopenblas-base default-jre
    curl -L https://github.com/claczny/VizBin/blob/master/VizBin-dist.jar?raw=true > VizBin-dist.jar
    mv VizBin-dist.jar /usr/local/bin # 或~/bin
    
    # 比对结果整理samtools
    conda install samtools
    
    ### CARD(选学) https://card.mcmaster.ca/download 
    # 方法1. 直接conda安装
    conda install --channel bioconda rgi
    # 方法2. 虚拟环境安装
    conda activate metawrap
    conda install --channel bioconda rgi
    rgi main -v # 4.0.3
    # rgi教程 https://github.com/arpcard/rgi

  

## 常见问题

### 数据库国内下载链接

- 方法1：国家微生物科学数据中心 —— 数据下载 (网站有时会更新和维护中，不可用)

http://nmdc.cn/datadownload

本资源由宏基因组平台发起，微生物所提供服务器，宏基因组团队负责维护的常用软件、扩增子和宏基因组数据库的国内下载链接，解决常用数据库下载慢、或无法下载的问题。同时提供定制的软件、数据库索引，节约大家下载时间，节省数据库编制索引的计算资源消耗。

- 方法2：刘永鑫Github中的备用百度云/下载链接索引

详细下载链接和说明见：https://github.com/YongxinLiu/MicrobiomeStatPlot/blob/master/Data/BigDataDownlaodList.md
    
下载的tar.gz压缩包，可放置于指定目录，使用`tar -xvzf *.tar.gz`解压 

个人网站备用下载链接，带宽有限，推荐官方、数据中心和百度云失败再尝试此法

    # humann3
    site=http://210.75.224.110/db
    wget -c ${site}/humann3/full_chocophlan.v296_201901.tar.gz
    wget -c ${site}/humann3/uniref90_annotated_v201901.tar.gz
    wget -c ${site}/humann3/full_mapping_v201901.tar.gz
    mkdir -p chocophlan uniref utility_mapping
    tar xvzf full_chocophlan.v296_201901.tar.gz -C chocophlan/
    tar xvzf uniref90_annotated_v201901.tar.gz -C uniref/
    tar xvzf full_mapping_v201901.tar.gz -C utility_mapping/

大文件的分卷压缩和解压 以kraken2 hash.k2d为例

    cd ~/db/kraken2_200706
    # https://www.cnblogs.com/wang--lei/p/9046643.html
    # 文件夹kraken2idx/打包压缩，1h
    tar -zcvf kraken2.tar.gz kraken2idx/
    # b分割为指定大小文件G/M/K，-d数字，a序列长度，输入和输出前缀
    split -b 13G -d -a 1 kraken2.tar.gz kraken2.tar.gz.
    # 一行命令
    tar -zcvf kraken2.tar.gz kraken2idx/ | split -b 19G -d -a 1 - kraken2.tar.gz.
    # 分割后解压缩
    cat kraken2.tar.gz.* | tar -zxv

    # 附tar打包用法，c创建、v输出过程，z压缩，f文件 ，x解压
    单个文件压缩打包 tar -cvzf my.tar.gz file1
    多个文件压缩打包 tar -cvzf my.tar.gz file1 file2 file*）（也可以给file*文件mv 目录在压缩）
    单个目录压缩打包 tar -cvzf my.tar.gz dir1
    多个目录压缩打包 tar -cvzf my.tar.gz dir1 dir2
    解包至当前目录：tar -xvzf my.tar.gz


### conda批量安装软件

    conda install -y fastqc multiqc kneaddata=0.6.1 humann2 graphlan export2graphlan lefse kraken2 megahit spades quast prokka cd-hit emboss salmon eggnog-mapper samtools
    
### kneaddata常见问题

#### kneaddata运行提示java版本不支持
    
	# 解决思路，新建虚拟环境，安装kneaddata，再安装对应java版本
	# 务必指定2.7，软件依赖2.7的python，但conda会自动安装3.6，运行报错
    conda create -n kneaddata python=2.7
    conda activate kneaddata
    conda install openjdk=8.0.152
    conda install kneaddata=0.6.1

#### kneaddata质控双端结果不成对
    
    # 0.7.4存在对旧格式fastq去宿主后数据极少或双端数据不对称，可指定版本安装0.6.1
    conda remove kneaddata
    conda install kneaddata=0.6.1 # 175 MB

### humann2数据库无法下载备用链接

    # 可选：链接直接下载、百度云链接(宏基因组公众号回复：数据库)或国内备份链接
    mkdir -p ${db}/humann2/chocophlan && cd ${db}/humann2/chocophlan
    wget -c http://210.75.224.110/share/meta/full_chocophlan_plus_viral.v0.1.1.tar.gz
    tar xvzf full_chocophlan_plus_viral.v0.1.1.tar.gz
    # uniref90和50任选其1，推荐uniref90更全5.9 GB
    cd ${db}/humann2
    wget -c http://210.75.224.110/share/meta/uniref90_annotated_1_1.tar.gz
    tar xvzf uniref90_annotated_1_1.tar.gz
    # 内存<32G内存选uniref5 (2.5 GB)
    # wget -c http://210.75.224.110/share/meta/uniref50_annotated_1_1.tar.gz
    # tar xvzf uniref50_annotated_1_1.tar.gz
    # 不要同一文件中有两个文件，会先比90，再比50出现混乱
    
### Metaphlan2数据库找不到

    # 下载
    mkdir -p ${db}/metaphlan2 && cd ${db}/metaphlan2
    wget -c http://210.75.224.110/db/humann2/metaphlan2.tar.gz
    tar xvzf metaphlan2.tar.gz
    # 链接到软件安装目录
    mkdir -p ${soft}/envs/metagenome_env/bin/databases
    ln -s ${db}/metaphlan2/* ${soft}/envs/metagenome_env/bin/databases/

### Kraken2

#### 定制数据库

官方教程详见 https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

设定数据库日期为版本，建立数据库目录

    d=200918
    cd ~/db
    mkdir -p kraken2/$d && cd kraken2/$d
    
下载物种名信息表，gb 1.8G; wgs 3.3G; 解压后大小为9.2/21G; taxdump 50M

    kraken2-build --download-taxonomy --threads 24 --db ./

序列数据库包括 archaea bacteria plasmid viral human fungi plant protozoa nr nt env_nr env_nt UniVec

下载单个真菌库为例

    i=fungi
    kraken2-build --download-library $i --threads 24 --db ./
    # 批量下载数据库，除默认5种外新加植物、真菌、原生生物和质粒，下载需几小时-数天 archaea bacteria UniVec_Core viral human fungi plant 
    for i in protozoa plasmid; do
    kraken2-build --download-library $i --threads 24 --db ./
    done
    
    # 建索引，4h, 40h；2h，40h
    time kraken2-build --build --threads 24 --db ./
    
    # 数据库大小，2020/4/12更新
    du -sh library/*
    918M    archaea
    74G     bacteria
    1.2G    fungi
    3.1G    human
    51G     plant
    2.0G    plasmid
    877M    protozoa
    2.0M    UniVec_Core
    310M    viral

#### Perl版本不对

常见问题：Perl版本不对，人工指定perl版本如下

    PERL5LIB=~/miniconda3/envs/kraken2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:~/miniconda3/envs/kraken2/lib/site_perl/5.26.2:~/miniconda3/envs/kraken2/lib/5.26.2/x86_64-linux-thread-multi:~/miniconda3/envs/kraken2/lib/5.26.2


### salmon手动安装和使用方法

    # 如不可用，尝试下载二进制和添加环境变量
    wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.0/salmon-0.14.0_linux_x86_64.tar.gz
    tar xvzf salmon-0.14.0_linux_x86_64.tar.gz 
    cp -rf salmon-latest_linux_x86_64/ ${soft}/envs/metagenome_env/share/salmon
    ${soft}/envs/metagenome_env/share/salmon/bin/salmon -v # 0.14.0

### MetaWRAP分箱

#### Bin可视化：BLAST Database error: Not a valid version 4 database

下载最新版的NCBI_nt，但与相应的blast版本不兼容，应该是本地软件版本太低。

方法1. 使用which config-metawrap配置为以前下载的旧版NCBI_nt库位置可用(2018年及以前下载)。

    # 备用nt下载链接
	mkdir -p NCBI_nt && cd NCBI_nt
    wget -c http://210.75.224.110/db/metawrap/NCBI_nt_181116/filelist.txt
    for a in `cat filelist.txt`; do \
	  wget -c http://210.75.224.110/db/metawrap/NCBI_nt_181116/$a; done
    for a in nt.*.tar.gz; do tar xzf $a; done
    # 配套tax下载链接
	mkdir -p NCBI_tax && cd NCBI_tax
    wget -c http://210.75.224.110/db/metawrap/NCBI_tax_181116/taxdump.tar.gz
	tar xzf taxdump.tar.gz

方法2，或升级blast为最新版。直接到https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/，下载最新版本blastn，再到conda的metawrap环境bin目录下，替换掉旧版本的blastn

## 附录

### Conda安装小工具

以下小工具已经整合至EasyMicrobiome项目中的linux文件夹，以下代码提供学习多种自主安装的参考方法，用于积累conda使用

并行计算管理rush/paprllel
    
    # conda安装rush，无依赖关系更好用的并行工具
    conda install rush -c bioconda
    # Ubuntu下安装方法 apt install parallel
    # conda安装parallel，版本有点老
    conda install parallel -c bioconda
    parallel --version # GNU parallel 20170422

表格统计工具csvtk和序列处理seqkit(可选中)

    # 方法1. conda安装，可能有点旧
    conda install csvtk -c bioconda
    conda install seqkit -c bioconda
    
    # 方法2. 直接下载最新版 https://github.com/shenwei356，如以csvtk为例手动安装
    wget -c https://github.com/shenwei356/csvtk/releases/download/v0.22.0/csvtk_linux_amd64.tar.gz
    tar xvzf csvtk_linux_amd64.tar.gz
    cp csvtk ~/miniconda3/bin/
    
### 宿主参考基因组下载

- EnsembleGenomes http://ensemblgenomes.org/ 
- 包括动物、植物、原生生物、真菌、细菌等，此外植物还 Phytozome https://phytozome-next.jgi.doe.gov/ ，以及单个物种和专用数据库

以Ensemble中拟南芥为例：Arabidopsis thaliana -- Genome assembly -- Download DNA sequence (无反应)，点TAIR链接跳转ENA，下载All Seq FASTA
    
    wget https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_000001735.1?download=true&gzip=true
    mv GCA_000001735.1\?download\=true TAIR10.fa

以Ensemble中水稻为例：Oryza sativa Japonica —— IRGSP-1.0
    
    wget https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001433935.1?download=true&gzip=true
    mv GCA_001433935.1\?download\=true IRGSP1.0.fa

### 软件和数据库清单(db)

扩增子、宏基因组分析常用软件、数据库、脚本文件


#### 软件

- linux：Linux系统下分析软件
    - microbiome_helper：metaphlan2结果转换STAMP格式(metaphlan_to_stamp)，picurst结果功能组成绘图(plot_metagenome_contributions.R)
    - miniconda3-latest-Linux-x86_64.sh：Conda安装程序
    - qiime2-2021.2.tar.gz：QIIME2安装包，解压至conda的envs目录可用
    - qiime2-2021.2-py36-linux-conda.yml：QIIME2软件安装清单，使用conda在线安装
    - sparcc.zip：sparcc网络分析python脚本
    - usearch：扩增子分析流程
    - vsearch：扩增子分析流程(免费64位版usearch)
- mac：Mac系统下分析软件
    - csvtk：表格分析工具
    - iqtree：进化树构建
    - qiime2-2021.2-py36-osx-conda.yml：QIIME2软件安装清单，使用conda在线安装
    - R-4.0.4.pkg：R语言安装包
    - RStudio-1.4.1106.dmg：RStudio安装包
    - rush：并行管理工具
    - seqkit：序列处理工具
    - taxonkit：NCBI分类处理工具
    - usearch：扩增子分析流程
    - vsearch：扩增子分析流程(免费64位版usearch)
- win：Windows系统下分析软件
    - STAMP2.1.3：微生物组图形界面差异分析工具
    - 4.0.zip：R语言常用400+包合集，解压至R包安装位置
    - Adobe_Illustrator_CC_2018_v22.1.0.314_x64_zh_CN_Portable.7z：图片拼图、模式图绘制工具
    - csvtk.exe：表格分析工具
    - Cytoscape_3_8_2_windows_64bit.exe：网络分析安装包
    - epp510_1828_64bit.exe：文本编辑器
    - FileZilla_3.49.1_win64_sponsored-setup.exe：文件上传下载
    - gephi-0.9.2-windows.exe：网络图绘制工具
    - Git-2.30.2-64-bit.exe：提供Git bash环境
    - iqtree.exe：进化树构建
    - jdk-11.0.7_windows-x64_bin.exe：Java运行环境
    - libiomp5md.dll：动态库，运行中提示缺少时可添加至系统目录
    - muscle.exe：多序列比对工具
    - npp.7.8.9.Installer.x64.exe：文本编辑器NotePad++安装包
    - R-4.0.4-win.exe：R语言安装包
    - RStudio-1.4.1106.exe：RStudio安装包
    - rtools40-x86_64.exe：R源码安装时的编绎工具
    - rush.exe：并行管理工具
    - seqkit.exe：序列处理工具
    - taxonkit.exe：NCBI分类处理工具
    - usearch.exe：扩增子分析流程
    - vsearch.exe：扩增子分析流程(免费64位版usearch)
    - wget.exe：命令行下载工具

#### 数据库

- gg：GreenGenes细菌16S数据库
    - gg_13_8_otus.tar.gz：13年8月更新OTU数据库，用于usearch有参定量和PICRUSt/BugBase功能预测、QIIME 2制作分类器
- kegg：KEGG数据库描述信息整理
    - KO_description.txt：KO编号对应的功能描述
    - KO_path.list：KO对应通路(Pathway)
    - ko00001.tsv：KO对应描述、通路、超级通路和分类信息
- usearch：usearch/vsearch物种分类sintax命令使用数据库
    - rdp_16s_v16_sp.fa.gz：16S的RDP16数据库，usearch作者整理，更16S、ITS和18S数据库见 http://www.drive5.com/usearch/manual/sintax_downloads.html
    - rdp_16s_v18.fa.gz：16S的RDP18数据库，易生信团队2021年基于RDP数据库整理
    - utax_reference_dataset_all_04.02.2020.fasta.gz：ITS注释数据库，可从UNITE下载

#### 脚本 

- 使用说明：分析常用脚本类型
    - .R文件为R脚本，使用Rscript命令执行；
    - .sh为Shell脚本，使用/bin/bash命令执行；
    - .pl为Perl脚本，使用perl命令执行；
    - .py为Python脚本，使用python执行，注意还分为python2和python3两种
- script：微生物组数据分析
    - BugBase：16S扩增子表型预测R脚本和数据库
    - FAPROTAX_1.2.4：16S扩增子元素循环预测Python脚本和数据库
    - table2itol：iTOL进化树注释文件制作R脚本
    - alpha_barplot.R：Alpha多样性指数柱状图+标准差图绘制
    - alpha_boxplot.R：Alpha多样性指数箱线图+统计绘制
    - alpha_rare_curve.R：usearch计算稀释曲线可视化
    - beta_cpcoa.R：基于距离矩阵开展限制性PCoA分析及可视化散点图+分组着色+置信椭圆，要求至少3个分组
    - beta_pcoa.R：基于距离矩阵的主坐标PCoA分析及可视化散点图+分组着色+置信椭圆+组间两两统计
    - BetaDiv.R：更多Beta多样性分析，如PCA、PCoA、NMDS、LDA、CCA、RDA等
    - compare.R：两组比较，支持t.test、wilcox、edgeR三种方法
    - compare_heatmap.R/sh：基于两组比较结果绘制热图
    - compare_manhattan.sh：基于两组比较结果绘制曼哈顿图
    - compare_volcano.R：基于两组比较结果绘制火山图
    - faprotax_report_sum.pl：FARPROTAX分析结果报告整理
    - filter_feature_table.R：按频率过滤OTU表
    - format_dbcan2list.pl：dbcan数据库注释结果整理
    - format2lefse.R：OTU表和物种注释生成LEfSe输入文件
    - format2stamp.R：OTU表和物种注释生成STAMP输入文件
    - kegg_ko00001_htext2tsv.pl：KEGG注释结果整理
    - kraken2alpha.R：Kraken2结果整理、抽平和alpha多样性指数计算
    - mat_gene2ko.R：按类型折叠表格
    - metaphlan_boxplot.R：metaphalan2结果可视化为箱线图
    - metaphlan_hclust_heatmap.R：metaphalan2结果可视化为聚类热图
    - metaphlan_to_stamp.pl：metaphalan2结果转换为STAMP格式
    - otu_mean.R：OTU表统计各组均值
    - otutab_filter_nonBac.R：16S的OTU表按sintax注释结果选择细菌、古菌且过滤叶绿体和线粒体
    - otutab_filter_nonFungi.R：ITS的OTU表选择真菌
    - otutab_freq2count.R：转换频率为伪整数，用于要求整型输入的分析，如多样性、edgeR差异分析等
    - otutab_rare.R：OTU表抽平
    - plot_metagenome_contributions.R：PICRUSt结果物种的功能组成绘制
    - sp_pheatmap.sh：绘制热图
    - sp_vennDiagram.sh：绘制维恩图
    - summarizeAbundance.py：按类型折叠大表，如基因按KEGG的KO合并
    - tax_circlize.R：物种组成圈图
    - tax_maptree.R：物种组成气泡图
    - tax_stackplot.R：物种组成堆叠柱状图

### 版本更新日志

#### 1.10 2021.1.22

1. humann2添加utility_mapping数据库，支持生成KEGG表；
2. kraken2添加最小8G索引；
3. 添加KEGG注释、层级信息及整理代码
4. 添加CARD数据库

#### 1.11 2021.5.7

1. 新增conda环境移植教程，并提供常用conda环境下载；
2. 新增软件、数据库、脚本说明文档
3. 常见问题整合到文档正文后面

#### 1.12 2021.8.20

1. 新增并行管理软件rush，比parallel更易安装，绿色版无依赖关系，整合在db/linux/目录中
2. 新增seqkit，可以统计序列数据量，支持序列长度过滤，格式转换等；
3. 新增质控软件fastp，软件fastqc更快，适合单独质控不去宿主；
4. kraken2的迷你库升级为2021年5月17日版本
5. eggNOG软件至eggNOG-mapper 2.1和配套数据库5.0.2

#### 1.13 2021.11.19

1. 新增软件conda环境下载安装方式，且作为首选
2. 新增kneaddata自定义物种基因组数据库示例

#### 1.14 2022.3.25

1. EasyMicrobiome升级为1.14
2. 升级miniconda2为miniconda3
3. dbcan2从2020/7/31的808M更新为2021/9/24版1016M，格式变化，配套format_dbcan2list.pl更新
4. 新增eggnog环境，包含emapper 2.1.6，summarizeAbundance.py含pandas (conda install sklearn-pandas)，配套更新数据库
5. rgi更新到最新版及配套代码
