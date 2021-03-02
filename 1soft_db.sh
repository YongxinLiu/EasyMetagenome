[TOC]

# 宏基因组软件和数据库安装 Metagenomic software & database

    # 测试环境为Linux Ubuntu 18.04 / CentOS 7
    # 版本: 1.10, 2021/1/22

## 安装前准备：Conda和数据库位置

    # 数据库安装位置，默认~/db目录(无需管理权限)，管理员可安装至/db，方便大家使用
    db=~/db
    mkdir -p ${db} && cd ${db}
    # 软件安装位置
    soft=~/miniconda3

### 软件管理器miniconda2

    # 下载最新版miniconda3
    wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # Conda默认软件安装目录为~/miniconda2，管理员可修改为/conda
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p ${soft}
    # 安装时许可协议打yes，-p指定安装目录为预定义的soft变量，注意安装完成后按提示激活
    ~/miniconda3/bin/conda init
    # 退出终端重新打开，提示符前出现(base)，方可使用conda
    conda -V # 查看版本 4.9.2
    python --version # 2.7.17

安装说明详见：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)

常用配置如下：

    # 添加生物频道，才能找到生物学软件
    conda config --add channels bioconda
    # http://bioconda.github.io/ 查询软件及版本
    # 添加常用频道
    conda config --add channels conda-forge

(可选)添加清华镜像加速下载(通常会加速，但有时会导致无法访问)

    site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
    conda config --add channels ${site}/pkgs/free/ 
    conda config --add channels ${site}/pkgs/main/
    conda config --add channels ${site}/cloud/conda-forge/
    conda config --add channels ${site}/pkgs/r/
    conda config --add channels ${site}/cloud/bioconda/

conda默认配置文件为 ~/.condarc 查看配置文件位置

    conda config --show-sources
    
可选：如果找不到conda，可手动临时添加conda至环境变量。可以添加至~/.bashrc文件中永久环境变量，注${soft}替换为你的安装目录，如~/miniconda2

    export PATH="${soft}/bin:$PATH"

查看虚拟环境列表 

    conda env list

创建虚拟环境，防污染环境变量，如果有的软件在Solving environment步骤数小时无法安装，可以新建环境

    conda create -n meta

加载环境

    conda activate meta

### 常用系统工具

(可选)Java运行环境jdk

    # example support trimmomatic
    conda install -c cyclus java-jdk # 156M，cyclus极慢
    java -version # 11.0.1

并行计算parallel
    
    # Ubuntu下安装方法
    # sudo apt install parallel
    # conda安装parallel，版本有点老
    conda install parallel -c bioconda
    parallel --version # GNU parallel 20170422

表格统计工具csvtk

    conda install csvtk -c bioconda

    
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
    # 数据库下载慢或失败，附录有百度云和国内备份链接
    
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

方法1. Conda导入和导出环境

    # 安装conda-pack
    conda install -c conda-forge conda-pack
    #conda 安装不成功可用pip安装
    # pip conda-pack
    # 安装好的环境下打包导出
    conda pack -n humann2 -o humann2.tar.gz

    # 新建文件夹存放humann2环境
    mkdir -p ~/miniconda3/envs/humann2
    tar -xzf humann2.tar.gz -C ~/miniconda3/envs/humann2
    # 激活环境
    source ~/miniconda3/envs/humann2/bin/activate

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
    humann2_config --update database_folders protein ${db}/humann2
    # metaphlan2数据库默认位于程序所在目录的db_v20和databases下各一份
    humann2_config --print

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

    conda install lefse # 76.3 MB, 1.0.8.post1
    
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

    # 从安装好的环境打包，此处为Ubuntu 16.04LTS
    # conda pack -n kraken2 -o ~/db/kraken2/kraken2.tar.gz
    # 添加权限，否则别人无法下载
    # chmod 777 ~/db/kraken2/kraken2.tar.gz
    # 下载压缩包
    wget -c http://210.75.224.110/db/kraken2/kraken2.tar.gz
    # 新建文件夹存放kraken2环境
    mkdir -p /conda2/envs/kraken2
    # 解压环境到指定目录
    tar -xzf kraken2.tar.gz -C /conda2/envs/kraken2
    # 激活环境
    source /conda2/envs/kraken2/bin/activate

- 数据库

下载数据库(NCBI每2周更新一次)，记录下载日期和大小。需根据服务器内存、使用目的选择合适方案。--standard标准模式下只下载5种数据库：古菌archaea、细菌bacteria、人类human、载体UniVec_Core、病毒viral。也可选直接下载作者构建的索引，还包括bracken的索引。

- 方法1. 数据库安装

方案1. 标准库安装，下载数据~100GB，时间由网速决定，索引5h，多线程可加速至1h完成
    
    cd ${db}
    d=210122
    mkdir -p kraken2/$d && cd kraken2/$d
    kraken2-build --standard --threads 24 --db ./
    
方案2. 自定义微生物数据库，如标准+真菌+原生动物+质粒

    cd ${db}
    d=210122
    mkdir -p kraken2/$d && cd kraken2/$d
    # 下载物种注释
    kraken2-build --download-taxonomy --threads 24 --db ./
    # 下载数据库
    for i in archaea bacteria UniVec_Core viral human fungi plasmid protozoa; do
        kraken2-build --download-library $i --threads 24 --db ./
    done
    # 确定的库建索引
    kraken2-build --build --threads 24 --db ./

- 方法1. 数据库下载

下载标准+原生动物+真菌+植物8GB(PlusPFP-8)数据库，包括kraken2和bracken2的索引。更多版本数据库详见：https://benlangmead.github.io/aws-indexes/k2 。

方案1. 迷你库(8G)

    # 压缩包5.2G，
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_8gb_20201202.tar.gz
    # 备用地址：
    # wget -c http://210.75.224.110/db/kraken2/k2_pluspfp_8gb_20201202.tar.gz
    tar xvzf k2_pluspfp_8gb_20201202.tar.gz

方案2. 完整库(8G)

    # 压缩包70G，
    d=20201202
    mkdir ${d} && cd ${d}
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20201202.tar.gz
    tar xvzf k2_pluspfp_20201202.tar.gz


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

    # 安装eggnog比对工具emapper
    conda install eggnog-mapper=2.0.1 -y -c bioconda
    emapper.py --version # 2.0.1
    
    # 下载常用数据库，注意设置下载位置
    mkdir -p ${db}/eggnog5 && cd ${db}/eggnog5
    # -y默认同意，-f强制下载，eggnog.db.gz 7.9G+4.9G
    download_eggnog_data.py -y -f --data_dir ./
    
    # 下载方式2(可选)：链接直接下载
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.0/eggnog.db.gz # 7.9G
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.0/eggnog_proteins.dmnd.gz # 4.9G
    gunzip *.gz
    
    # 如果内存够大，复制eggNOG至内存加速比对
    # cp eggnog.db /dev/shm
    # 手工整理COG分类注释
    wget -c http://210.75.224.110/share/COG.anno
    # 手工整理KO注释
    wget -c http://210.75.224.110/share/KO.anno



### CAZy碳水化合物

    # dbCAN2 http://bcb.unl.edu/dbCAN2
    # 创建数据库存放目录并进入
    mkdir -p ${db}/dbCAN2 && cd ${db}/dbCAN2
    # 下载序列和描述
    wget -c http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa
    wget -c http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07302020.fam-activities.txt
    # 备用数据库下载地址并解压 
    #wget -c http://210.75.224.110/db/dbcan2/CAZyDB.07312020.fa.gz
    #gunzip CAZyDB.07312020.fa.gz
    # diamond建索引，800M，1m
    diamond --version # 0.8.22/2.0.5
    time diamond makedb \
      --in CAZyDB.07312020.fa \
      --db CAZyDB.07312020
    # 压缩原始数据节约空间
    gzip CAZyDB.07312020.fa
    # 提取fam对应注释
    grep -v '#' CAZyDB.07302020.fam-activities.txt \
      |sed 's/  //'| \
      sed '1 i CAZy\tDescription' \
      > CAZy_description.txt

### 抗生素抗性基因CARD

    # 官网：https://card.mcmaster.ca
    # Bioconda: http://bioconda.github.io/recipes/rgi/README.html
    # Github: https://github.com/arpcard/rgi

软件安装

    # 方法1. Cona安装rgi
    conda install --channel bioconda rgi
    
    # 方法2. 指定环境安装rgi
    conda create -n rgi -c bioconda rgi
    
    # 方法3. 从其他安装的环境导出和导入
    # 安装好的环境下打包导出，262M
    conda pack -n rgi -o rgi.tar.gz
    # 下载软件包
    wget -c http://210.75.224.110/db/card/rgi.tar.gz
    # 新建文件夹存放rgi环境
    mkdir -p ~/miniconda2/envs/rgi
    tar -xzf rgi.tar.gz -C ~/miniconda3/envs/rgi
    # 激活环境
    source ~/miniconda2/envs/rgi/bin/activate

数据库部署

    # 下载最新版数据库，2.8M
    wget -c https://card.mcmaster.ca/latest/data
    # 解压后20M
    tar -xvf data card.json
    # wget -c http://210.75.224.110/db/card/card.json
    # 加载数据库
    rgi load --card_json card.json
    
    # 宏基因组分析扩展数据库和加载
    rgi card_annotation -i card.json
    rgi load -i card.json --card_annotation card_database_v3.1.0.fasta

    
### 抗生素抗性基因Resfam

    # http://dantaslab.wustl.edu/resfams
    mkdir -p ${db}/resfam && cd ${db}/resfam
    # 官网的数据格式非常混乱, 推荐下载我手动整理的索引和注释
    wget http://210.75.224.110/share/Resfams-proteins.dmnd # 1.5 MB
    wget http://210.75.224.110/share/Resfams-proteins_class.tsv # 304 KB


## 5分箱工具

### metawrap分箱和提纯

    # 物种注释和分箱流程 https://github.com/bxlab/metaWRAP
    conda create -y -n metawrap python=2.7 # 22.2MB
    conda activate metawrap
    conda config --add channels ursky
    conda install -y -c ursky metawrap-mg # 1.14 GB, v1.2
            
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
    # 配置文件通常为~/miniconda2/envs/metawrap/bin/config-metawrap
    # 使用Rstudio/vim等文本编辑器来修改数据库的位置

### derep基因组去冗余

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

    # 按注释类型合并脚本，链接 https://github.com/YongxinLiu/Metagenome
    cd ${db}
    mkdir script && cd script
    # 矩阵按最后列注释合并的脚本，如Gene-KO-Module-Pathway的合并
    wget http://210.75.224.110/share/meta/script/mat_gene2ko.R
    # 基于Kraken2的结果计算alpha多样性
    wget http://210.75.224.110/share/meta/script/kraken2alpha.R
    # 基于alpha多样性和分组信息绘制箱线图
    wget http://210.75.224.110/share/meta/script/alpha_boxplot.R
    chmod +x *.R


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

  

## 附录

### 数据库国内下载链接

- 方法1：国家微生物科学数据中心 —— 数据下载

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

    # 附tar打包用法
    单个文件压缩打包 tar -czvf my.tar.gz file1
    多个文件压缩打包 tar -czvf my.tar.gz file1 file2 file*）（也可以给file*文件mv 目录在压缩）
    单个目录压缩打包 tar -czvf my.tar.gz dir1
    多个目录压缩打包 tar -czvf my.tar.gz dir1 dir2
    解包至当前目录：tar -xzvf my.tar.gz


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

    # 下载并构建索引
    mkdir -p ${db}/metaphlan2 && cd ${db}/metaphlan2
    wget -c http://210.75.224.110/share/meta/metaphlan2/mpa_v20_m200.tar
    tar xvf mpa_v20_m200.tar
    bzip2 -d mpa_v20_m200.fna.bz2
    bowtie2-build mpa_v20_m200.fna mpa_v20_m200
    # 链接到软件安装目录
    mkdir -p ${soft}/envs/metagenome_env/bin/db_v20
    ln -s ${db}/metaphlan2/* ${soft}/envs/metagenome_env/bin/db_v20/
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

    PERL5LIB=~/miniconda2/envs/kraken2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:~/miniconda2/envs/kraken2/lib/site_perl/5.26.2:~/miniconda2/envs/kraken2/lib/5.26.2/x86_64-linux-thread-multi:~/miniconda2/envs/kraken2/lib/5.26.2


### salmon手动安装和使用方法

    # 如不可用，尝试下载二进制和添加环境变量
    wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.0/salmon-0.14.0_linux_x86_64.tar.gz
    tar xvzf salmon-0.14.0_linux_x86_64.tar.gz 
    cp -rf salmon-latest_linux_x86_64/ ${soft}/envs/metagenome_env/share/salmon
    ${soft}/envs/metagenome_env/share/salmon/bin/salmon -v # 0.14.0


## 1.10 2021.1.22

1. humann2添加utility_mapping数据库，支持生成KEGG表；
2. kraken2添加最小8G索引；
3. 添加KEGG注释、层级信息及整理代码
4. 添加CARD数据库
