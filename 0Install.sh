[TOC]

# EasyMetagenome software & database (易宏基因组软件和数据库)

    # Authors(作者): Yong-Xin Liu(刘永鑫), Defeng Bai(白德凤), Tong Chen(陈同) et al.
    # Version(版本): 1.24, 2025/11/12
    # Homepage(主页): https://github.com/YongxinLiu/EasyMetagenome

    # All software and databases downloaded from the website (所有软件和数据库可从官网下载)
    # Backup source add download speed and success rate (备用站点提高下载速度和成功率)
    # Backup1. Institute of Microbiology, Chinese Academy of Sciences(中科院微生物所)：ftp://download.nmdc.cn/tools/ (FileZilla访问) 
    # Backup2. Baidu Netdisk(百度网盘)：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315

# 0. Initialization and pipeline installation (零、初始化和流程安装 )

## Initialization: The following code must be run every time when installation starts (初始化：每次开始安装必须运行下面代码)

    # Software and Database Locations(软件和数据库位置)
    # The follwoing paragraph must run before(分析前必须运行)
    # Database Locations, default ~/db directory(No administrative privileges required), administrator can select /db
    # 数据库安装位置，默认~/db目录(无需管理权限)，管理员可选/db
    db=~/db
    mkdir -p ${db} && cd ${db}
    # Software installation location, default ~/miniconda3 , the test server is /anaconda3
    # 软件安装位置，默认为 ~/miniconda3 ，测试服务器为 /anaconda3
    soft=~/miniconda3
    # In the frequently used server environment, you can replace the variable ${db} and ${soft} with absolute paths, and you will no longer need to run the above environment variables every time
    # 经常使用的服务器环境，可把全文${db}和${soft}替换为绝对路径，将不再需要每次读取以上环境变量
    # Optional: Initialize environment variables, which may improve the success rate of software installation
    # 可选：初始化环境变量，可能提高软件安装成功率
    PATH=${soft}/bin:${soft}/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${db}/EasyMicrobiome/linux:${db}/EasyMicrobiome/script
    echo $PATH

### EasyMetagenome Pipeline (流程)

    # The EasyMetagenome workflow, including workflow installation, usage, and visualization scripts, 
    # as well as workflow test data and result comparisons, can be found at:
    # EasyMetagenome流程，包括流程安装、使用和可视化脚本，以及流程测试数据和结果正对照，
    # 网址：https://github.com/YongxinLiu/EasyMetagenome
    # Each site provides 2-4 download or installation methods: try them one by one and any one method work is OK.
    # 每处提供2-4种下载或安装方法：依次尝试，任意一种成功即可。
    
    cd ${db}

    # Method 1. Download ZIP archive at https://github.com/YongxinLiu/EasyMetagenome, `Code` - `Download` and upload it to the server.
    # 方法1. 网页 https://github.com/YongxinLiu/EasyMetagenome 中Code - Download ZIP下载压缩包，上传至服务器
    # Command 'unzip' not found, using `sudo apt install unzip` (解压，缺少使用sudo apt安装)
    unzip EasyMetagenome-master.zip 
    # rename 改名
    mv EasyMetagenome-master EasyMetagenome

    # Method 2. Baidu Netdisk (方法2. 百度网盘) https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 /db/soft/EasyMetagenome.tar.gz

    # Method 3. Download by Git (needs Git installed)（方法3. git下载，需安装git）
    git clone https://github.com/YongxinLiu/EasyMetagenome
    # Old version update (旧版更新)
    cd EasyMetagenome && git pull && cd ../

    # Method 4. Backup link for the Institute of Microbiology (方法4. 微生物所备用链接)
    wget -c ftp://download.nmdc.cn/tools/soft/EasyMetagenome.tar.gz
    tar xvzf EasyMetagenome.tar.gz

### EasyMicrobiome: dependencies soft and scripts for pipeline (软件和数据库合集)

    # EasyMetagenome depends on EasyMicrobiome, which includes a collection of scripts, 
    # software, and databases. (https://github.com/YongxinLiu/EasyMicrobiome)
    # EasyMetagenome依赖EasyMicrobiome，其包括众多脚本、软件和数据库的集合，
    # 网址：https://github.com/YongxinLiu/EasyMicrobiome
    
    # Method 1. Download ZIP archive at https://github.com/YongxinLiu/EasyMicrobiome, `Code` - `Download` and upload it to the server.
    # 方法1. 网页 https://github.com/YongxinLiu/EasyMicrobiome 中Code - Download ZIP下载压缩包，上传至服务器
    unzip EasyMicrobiome-master.zip
    mv EasyMicrobiome-master EasyMicrobiome
    
    # Method 2. Baidu Netdisk (方法2. 百度网盘) https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 /db/soft/EasyMicrobiome.tar.gz
   
    # Method 3. Download by Git（方法3. git下载）
    git clone https://github.com/YongxinLiu/EasyMicrobiome
    # Old version update (旧版更新)
    cd EasyMicrobiome && git pull && cd ../

    # Method 4. Backup link for the Institute of Microbiology (方法4. 微生物所备用链接)
    wget -c ftp://download.nmdc.cn/tools/soft/EasyMicrobiome.tar.gz
    tar -xvzf EasyMicrobiome.tar.gz
     
    # Add command execute permissions (添加命令可执行权限)
    chmod +x `pwd`/EasyMicrobiome/linux/* `pwd`/EasyMicrobiome/script/*
    # Add environment variables (添加环境变量)
    echo "export PATH=\"\$PATH:`pwd`/EasyMicrobiome/linux:`pwd`/EasyMicrobiome/script\"" >> ~/.bashrc
    source ~/.bashrc
    echo $PATH

### Conda: Software Manager (软件管理器)

    # Downloaded the latest version of miniconda3 v25.9.1, installed on 2025/10/22, 154.6 Mb
    # 下载最新版miniconda3 v25.9.1 , 安装日期2025/10/22, 154.6 Mb   
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # Install, -b batch, -f no prompt, -p directory, select yes for license agreement
    # 安装，-b批量，-f无提示，-p目录，许可协议打yes
    bash Miniconda3-latest-Linux-x86_64.sh -b -f 
    # Initialize and reset the environment; success show (base) at the prompt.
    # 初始化，并重置环境，提示符前出现(base)即成功
    ~/miniconda3/condabin/conda init
    source ~/.bashrc
    # Show version(查看版本), conda 25.9.1, python 3.13.9
    conda -V
    python --version
    # Add frequently used channels (添加常用频道)
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
    # conda configure file (配置文件) ~/.condarc
    conda config --show-sources

    # Install mamba parallel acceleration installation, pandas computing dependencies, and conda-pack packaging and unpack enrivonment.
    # 安装mamba并行加速安装、pandas计算依赖库和conda-pack打包移植环境
    conda install mamba -c conda-forge -c bioconda -y
    mamba install pandas -c conda-forge -c bioconda -y
    mamba install conda-pack -c conda-forge -c bioconda -y

    # View environment list （查看虚拟环境列表)
    conda env list


# 1. Data preprocessing (一、数据预处理)

## kneaddata install (安装): 

    # Note: You can choose one of the following installation methods: direct installation, download and unpack, etc. If one method fails, try another.
    # 注：直接安装、下载解压等安装方法多选一。一种方法不成功，再尝试另一种。
    # BioConda search software: https://bioconda.github.io/recipes/kneaddata/README.html

### Opt 1. Download, extract, and install kneaddata (方法1.kneaddata下载解压安装)

    # Specify conda filename (指定conda文件名)
    s=kneaddata
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    # Set installation directory and extract (指定安装目录并解压)
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    # Startup Environment (启动环境)
    conda activate ${s}
    # Initialize environment (初始化环境)
    conda unpack

### Opt 2. Conda kneaddata (方法2. conda安装kneaddata)

    # Create and activate the kneaddata environment (新建并激活环境)
    conda create -y -n kneaddata
    conda activate kneaddata
    # install KneadData host removal, FastQC quality assessment, and MultiQC summary
    # 安装kneaddata去宿主流程，fastqc质量评估，multiqc评估报告汇总
    mamba install kneaddata fastqc multiqc r-reshape2 -y 

### Record software version (记录软件版本)

    fastp --version # 1.0.1
    kneaddata --version # 0.12.3
    bowtie2 --version # 2.5.4
    fastqc -v # v0.12.1
    trimmomatic -version # 0.40
    multiqc --version  # 1.321
    
    # Optional: Software packaging -- can copy to others (可选：软件打包--复制解压使用)
    conda pack -f --ignore-missing-files -n ${s} -o ${s}.tar.gz

### kneaddata database download (数据库下载)

    # View available databases (查看可用数据库)
    kneaddata_database
    # Including human genome/transcriptome, ribosomal RNA, mouse/dog/cat genome
    # 包括人基因组/转录组、核糖体RNA、小鼠/狗/猫基因组

human genome download (人类基因组下载)

    mkdir -p ${db}/kneaddata/human

    # Opt 1. Automatically download and extract the human T2T genome bowtie2 index (3.6 GB).
    # 方法1. 自动下载解压人类T2T基因组bowtie2索引 3.6 GB
    kneaddata_database --download human_genome bowtie2 ${db}/kneaddata/human
    
    # Opt 2. Manually download the Human Genome from the official, NMDC, or Baidu Netdisk, then extract 1 minute
    # 方法2. 手动下载官网、备用链接或百度云人类基因组并解压1分钟
    cd ${db}/kneaddata/human
    wget -c https://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz
    # wget -c ftp://download.nmdc.cn/tools/meta/kneaddata/human/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz
    time tar xvzf Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz
    cd $db

mouse genome download (小鼠基因组下载)
    
    mkdir -p ${db}/kneaddata/mouse

    # Opt 1. Automatically download and extract mouse genome bowtie2 index 2.83 GB
    # 方法1. 下载小鼠基因组bowtie2索引 2.83 GB
    kneaddata_database --download mouse_C57BL bowtie2 ${db}/kneaddata/mouse
    
    # Opt 2. Manually download the Human Genome from the official, NMDC, or Baidu Netdisk, then extract 1 minute
    # 方法2. 手动下载官网、备用链接或百度云人类基因组并解压1分钟
    cd ${db}/kneaddata/mouse
    # wget -c http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz
    wget -c ftp://download.nmdc.cn/tools/meta/kneaddata/mouse/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz
    time tar xvzf mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz
    cd $db
    
### kneaddata Custom reference genome index (自定义参考基因组索引)

    # **Reference genome index by Bowtie2, if there are multiple hosts, the reference genome sequences can be merged then index.**
    # **参考基因组采用bowtie2索引即可；如果有多个宿主可以基因组序列合并后构建索引**
    # Part of the genome can be downloaded from the ensembl genome. Using Arabidopsis thaliana as an example, 
    # visit http://plants.ensembl.org/index.html, select *Arabidopsis thaliana* — Download DNA sequence (FASTA), 
    # select toplevel, right-click and copy the link, then paste it into the link field below.
    # 部分基因组可在ensembl genome下载。此处以拟南芥为例，访问 http://plants.ensembl.org/index.html ，
    # 选择Arabidopsis thaliana —— Download DNA sequence (FASTA)，选择toplevel右键复制链接，填入下面链接处

    conda activate kneaddata
    # Create directories(创建目录)
    mkdir -p ${db}/kneaddata/ath
    cd ${db}/kneaddata/ath
    # Download host genome (下载宿主基因组）
    wget -c http://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz tair10.fa.gz
    # unzip and index, 3 minutes (解压，建索引，120M 3分钟)
    gunzip tair10.fa.gz
    time bowtie2-build -f tair10.fa tair10 --threads 4
    cd ${db}

# 2. Read-based HUMAnN4/Kraken2 (二、基于读长分析)

    # HUMAnN v4.0.0.alpha.1 + MetaPhlAn4 v4.1.1, HUMAnN3and HUMAnN2 see appendix

## Read-based HUMAnN4/MetaPhlAn4 (基于读长的分析)

### Opt 1. Download install HUMAnN4 (方法1. HUMAnN4下载解压安装)

    # Specify conda filename (指定conda文件名)
    s=humann4
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选微生物所nmdc 或 百度云/db/conda/humann4.tar.gz)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    time tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s} # 1m
    conda activate ${s}
    conda unpack
    
### Method 2. Conada installation of HUMAnN4 (方法2. HUMAnN4直接安装)

    # install (安装)
    conda create -n humann4
    conda activate humann4
    conda install -c biobakery humann=4.0.0a1
    conda install -c bioconda metaphlan=4.1.1

    # Opt. Package (打包)
    n=humann4
    conda pack -f --ignore-missing-files -n ${n} -o ${n}.tar.gz

### Method 3. Installing HUMAnN4 from source code (方法3. HUMAnN4源码安装)

    # Prepare the environment (准备环境)
    conda create -n humann4 python=3.9 -y
    conda activate humann4
    conda install -c bioconda bowtie2 diamond blastp
    # clone & install
    git clone https://github.com/biobakery/humann.git
    cd humann
    git checkout humann-4.0.0.alpha.1   # 或最新分支
    python setup.py install
    cd ..

### HUMAnN4 Installation test (安装测试)

    # Record core software version (记录核心软件版本)
    humann --version # v4.0.0.alpha.1
    metaphlan -v # 4.1.1 (11 Mar 2024)
    diamond help | head -n 1 # v2.1.14.168
    # test (测试)
    humann_test

### HUMAnN4 database (数据库)

    # database directory (建立数据库安装目录 )
    mkdir -p ${db}/humann4 && cd ${db}/humann4
    mkdir -p chocophlan chocophlan_ec uniref utility_mapping
    # 显示可用分类、泛基因组和功能数据库
    humann_databases

    # Method 1. Download and extract the file using wget or Baidu (方法1. wget或百度下载并解压)
    # chocophlan full 42G, 5-7min
    wget -c http://huttenhower.sph.harvard.edu/humann_data/chocophlan/chocophlan.v4_alpha.tar.gz
    time tar xvzf chocophlan.v4_alpha.tar.gz -C ${db}/humann4/chocophlan
    # chocophlan ec_filtered, 6.5G, 1min
    wget -c http://huttenhower.sph.harvard.edu/humann_data/chocophlan/chocophlan_EC_FILTERED.v4_alpha.tar.gz
    time tar xvzf chocophlan_EC_FILTERED.v4_alpha.tar.gz -C ${db}/humann4/chocophlan_ec
    # uniref, 893M, 13s 1.6G
    wget -c http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_ec_filtered/uniref90_annotated_v4_alpha_ec_filtered.tar.gz
    time tar xvzf uniref90_annotated_v4_alpha_ec_filtered.tar.gz -C ${db}/humann4/uniref
    ls -lh ${db}/humann4/uniref/humann4_protein_database_filtered_v2019_06.dmnd
    # uniref, 2.7G, 30s
    wget -c http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v4_alpha.tar.gz
    time tar xvzf full_mapping_v4_alpha.tar.gz -C ${db}/humann4/utility_mapping

    # Method 2. The program install database (方法2. 程序自动安装数据库)
    # microbial pangenome(微生物泛基因组) 42 GB
    humann_databases --download chocophlan full ${db}/humann4
    # icrobial pangenome ec_filtered(微生物泛基因组ec过滤版) 6.5 GB
    humann_databases --download chocophlan ec_filtered ${db}/humann4
    # uniref90 diamond index (功能基因索引) 1.6G
    humann_databases --download uniref uniref90_ec_filtered_diamond ${db}/humann4 
    # functional description (功能描述) 2.7 GB
    humann_databases --download utility_mapping full ${db}/humann4

    # Set database location (设置数据库位置)
    # Display parameters (显示参数)
    humann_config --print
    # Modify the number of threads (修改线程数)
    humann_config --update run_modes threads 8
    # Set the locations for nucleic acids, proteins, and annotation database (设置核酸、蛋白和注释库位置)
    humann_config --update database_folders nucleotide ${db}/humann4/chocophlan
    humann_config --update database_folders protein ${db}/humann4/uniref
    humann_config --update database_folders utility_mapping ${db}/humann4/utility_mapping
    # Check settings results(核对设置结果)
    humann_config --print

### MetaPhlAn4 taonomic database (物种数据库)

    # humann4 need MetaPhlAn4 mpa_vOct22_CHOCOPhlAnSGB_202403 database (2.98G+19.87G)
    mkdir -p ${db}/metaphlan4 && cd ${db}/metaphlan4

    # Opt 1. Download from official website (官网下载) 2.7G, 30s; 19G, 5m;
    wget -c http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202403.tar
    time tar xvf mpa_vOct22_CHOCOPhlAnSGB_202403.tar
    wget -c http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar
    time tar xvf mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar

    # Opt 2/3. Download from Baidu NetDisk or NMDC (下载可选百度云或微生物所)
    # Baidu NetDisk: https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 db/meta/metaphlan4
    wget -c ftp://download.nmdc.cn/tools/meta/metaphlan4/mpa_vOct22_CHOCOPhlAnSGB_202403.tar
    wget -c ftp://download.nmdc.cn/tools/meta/metaphlan4/mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar.gz
    tar xvf mpa_vOct22_CHOCOPhlAnSGB_202403.tar
    time tar xvzf mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar.gz

### GraPhlAn taxonomic and phylogenetic trees (高颜值物种或进化树)

    # Method 1. Unzip installation conda package
    # 方法1. Conda包解压安装
    cd $db
    n=graphlan
    # Download from NMDC ftp://download.nmdc.cn/tools/conda or BaiduNetdisk db/conda from https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    mkdir -p ${soft}/envs/${n}
    tar -xvzf ${n}.tar.gz -C ${soft}/envs/${n}
    # activate and initial conda enviroment (启动和初始化环境)
    conda activate ${n}
    conda unpack

    # Method 2. Conda installation
    # 方法2. Conda安装
    conda create -n graphlan graphlan export2graphlan -c bioconda -y
    graphlan.py --version # GraPhlAn version 1.1.3 (5 June 2018)

## LEfSe Biomarker identification and visualization (生物标记鉴定和可视化)

    # Method 1. ImageGP 2 https://www.bic.ac.cn/BIC/#/analysis?page=b%27MzY%3D%27&tool_type=tool

    # Method 2. Download install LEfSe (方法2. LEfSe下载解压安装)
    n=lefse
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    mkdir -p ${soft}/envs/${n}
    tar -xvzf ${n}.tar.gz -C ${soft}/envs/${n}
    conda activate ${n}
    chmod 755 ${soft}/envs/lefse/lib/python3.9/site-packages/lefse/*
    conda unpack
    lefse_run.py -h # LEfSe 1.1.01

    # Method 3. Conda installation (方法3. Conda安装)
    mamba create -n lefse lefse -c bioconda -y


## Kraken2 taxonomic classification (物种注释)

    # kraken2: https://ccb.jhu.edu/software/kraken/

    # Method 1. Kraken2 Conda package download install (包本地解压安装)
    n=kraken2
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选微生物所nmdc 或 百度云/db/conda/kraken2.tar.gz)
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.gz
    mkdir -p ${soft}/envs/${n}
    time tar -xvzf ${n}.tar.gz -C ${soft}/envs/${n}
    conda activate ${n}
    conda unpack

    # Method 2. Conada installation of kraken2 (方法2. Kraken2在线安装指定版本)
    n=kraken2
    mamba create -n ${n} -y -c bioconda kraken2=2.1.6 python=3.9
    conda activate ${n}
    mamba install bracken krakentools krona r-optparse -y
    kraken2 --version # 2.1.6
    less `type bracken | cut -f2 -d '('|cut -f 1 -d ')'`|grep 'VERSION' # 2.9
    conda pack -f --ignore-missing-files -n ${n} -o ${n}.tar.gz

### Kraken2 database (数据库安装)

    # Database：https://benlangmead.github.io/aws-indexes/k2  
    # set version, current 20250714, set each type directory
    db=~/db
    v=20250714
    mkdir -p ${db}/kraken2 && cd ${db}/kraken2
    mkdir -p pluspf16g pluspf pluspfp

    # Option 1. Download the standard + protozoa + fungi, compressed 11.2GB, uncompressed 14.9GB (PlusPF-16)
    # 方案1. 下载标准+原生动物+真菌，压缩包11.2G，解压14.9GB (PlusPF-16) 
    # Download from offical, NMDC or Baidu NetDisk (下载自官网、微生物所nmdc或百度云)
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_${v}.tar.gz
    # wget -c ftp://download.nmdc.cn/tools/meta/kraken2/k2_pluspf_16_GB_${v}.tar.gz
    time tar xvzf k2_pluspf_16_GB_${v}.tar.gz -C pluspf16g # 1min

    # Option 2. Download the standard + protozoa + fungi, compressed 77.5G, uncompressed 100.6G (PlusPF)
    # 方案2. 下载标准+原生动物+真菌，压缩包77.5G，解压100.6G (PlusPF) 
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_${v}.tar.gz
    # wget -c ftp://download.nmdc.cn/tools/meta/kraken2/k2_pluspf_${v}.tar.gz
    time tar xvzf ${db}/kraken2/k2_pluspf_${v}.tar.gz -C ~/db/kraken2/pluspf # 1min
    
    # Option 3. Download full for standards, protozoa, fungi, and plants; 
    # the compressed file is 158.8GB, and the uncompressed file is 214.5GB (PlusPFP).
    # 方案3. 下载标准+原生动物+真菌+植物完整库，压缩包158.8G，解压214.5G (PlusPFP) 
    wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_${v}.tar.gz
    # wget -c ftp://download.nmdc.cn/tools/meta/kraken2/k2_pluspf_${v}.tar.gz
    time tar xvzf ${db}/kraken2/k2_pluspfp_${v}.tar.gz -C pluspfp # 6min


# 3. Assemble-based (三、组装分析流程)

## Assembly and quantification: megahit/spades/prodigal/cd-hit/salmon (组装和定量)

    cd $db
    s=megahit
    
    ### Opt 1. Download, extract, and install megahit (方法1.megahit下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/megahit
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack

    ### Opt 2. Conda install (方法2. conda安装)
    mamba create -y -n megahit megahit spades quast cd-hit emboss salmon prodigal
    conda activate megahit 

    ### test(测试)
    megahit -v # MEGAHIT v1.2.9
    metaspades.py -v # metaSPAdes v4.2.0
    metaquast.py -v # MetaQUAST v5.3.0
    cd-hit -v | grep version # CD-HIT v4.8.1 (built on Apr 24 2025) 
    embossversion # EMBOSS v6.6
    salmon -v # salmon v1.10.3
    # Package(打包)
    conda pack -f --ignore-missing-files -n ${s} -o ${s}.tar.gz

## eggNOG protein annotation (蛋白注释)

    # eggNOG: http://eggnogdb.embl.de
    s=eggnog

    ### Opt 1. Download and install eggNOG (方法1.eggNOG下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack

    ### Opt 2. Conda install (方法2. conda安装)
    mamba create -n eggnog -y
    conda activate eggnog
    mamba install eggnog-mapper -y -c bioconda -c conda-forge

### eggNOG database (数据库)

    # Set download directory (下载数据库目录建立并进入)
    mkdir -p ${db}/eggnog && cd ${db}/eggnog

    ### Opt 1. Download and install database (方法1.数据库下载解压)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/meta/eggnog/eggnog.tar.gz
    tar xvzf eggnog.tar.gz

    ### Opt 2. Office install install (方法2. 软件脚本下载)
    # : -y yes, -f force, eggnog.db.gz 6.3G+4.9G，解压后48G
    download_eggnog_data.py -y -f --data_dir ${db}/eggnog
    
    # Option. package database for others (数据库打包)
    tar czvf 5.0.2.tar.gz eggnog*
    mv 5.0.2.tar.gz eggnog.tar.gz

    # link to default directory （链接至默认目录，注意按实际情况修改)
    mkdir -p ${soft}/envs/eggnog/lib/python3.11/site-packages/data/
    ln -sf `pwd`/* ${soft}/envs/eggnog/lib/python3.11/site-packages/data/

    # Copying database to memory speeds up (复制数据至内存中加速比对)
    cp eggnog.* /dev/shm
    cd $db
    
    ### test(测试), 2.1.13, expected eggNOG DB version: 5.0.2, diamond version 2.0.15
    emapper.py --version 


## CARD/rgi antibiotic resistance gene (抗生素抗性基因)

    # CARD：https://card.mcmaster.ca
    # RGI Github: https://github.com/arpcard/rgi
    s=rgi

    ### Opt 1. Download and install RGI (方法1. RGI下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack
    
    ### Opt 2. rgi conda install
    mamba create -y -n rgi rgi=6.0.5
    conda activate rgi
    conda pack -f --ignore-missing-files -n ${s} -o ${s}.tar.gz

    # view version (查看版本) 6.0.5
    rgi main -v
    
    ### rgi database (数据库部署)
    mkdir -p ${db}/card && cd ${db}/card
    # Download latest database, 4.41M (2025-11-1, 4.0.1)
    wget -c https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    rgi load --card_json card.json
    # Metagenomic analysis expands database and loads (宏基因组分析扩展数据库和加载)
    rgi card_annotation -i card.json
    mv card_database_v4.0.1_all.fasta card.fasta
    rgi load -i card.json --card_annotation card.fasta
    

# 4. Binning (四、分箱挖掘单菌基因组)

## 4.1 MetwWRAP binning (分箱)

    # GitHub: https://github.com/bxlab/metaWRAP
    s=metawrap
    cd ${db} 

    ### Opt 1. Download and install MetwWRAP (方法1.MetwWRAP下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack

    ### Opt 2. Conda install (方法2. conda安装)
    mamba create -y --name metawrap --channel ursky -c conda-forge -c bioconda metawrap-mg=1.3.2
    conda activate metawrap
    conda pack -f --ignore-missing-files -n ${s} -o ${s}.tar.gz

    metawrap -h # 1.3.2
    
    # CheckM for Bin integrity and contamination estimation and species annotation
    # CheckM用于Bin完整和污染估计和物种注释
    mkdir -p $db/checkm && cd $db/checkm
    # download 275 MB, unzip 1.4 GB
    wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xvf *.tar.gz
    # Set database location, Enterx2 (设置数据库位置，直接2次回车默认当前)
    checkm data setRoot

    # NCBI taxonomy (物种信息) zip 68M, unzip 351M
    mkdir -p ${db}/NCBI/tax && cd ${db}/NCBI/tax
    wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xvzf taxdump.tar.gz
    
    # (Optional) NCBI Nucleic Acid and Species Information (可选)NCBI核酸和物种信息
    mkdir -p ${db}/NCBI/nt && cd ${db}/NCBI/nt
    wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz
    cd ${db}/NCBI/nt; for i in *.tar.gz; do tar xzf $i; done
    # Some libraries may not download completely. Delete them and download them again; do not resume downloading.
    # 可能会出现个别库下载不完整的情况，删了重下，不要续传
    
    ## Database location settings (数据库位置设置)
    which config-metawrap
    # usually in ~/miniconda3/envs/metawrap/bin/config-metawrap
    # Using Rstudio/vim to set database location

## 4.2 dRep genome redundance (基因组去冗余)

    # GitHub: https://github.com/MrOlm/drep
    # Conda: https://bioconda.github.io/recipes/drep/README.html
    # The genomes obtained from single-strain sequencing contain a significant amount of redundancy. 
    # Metawrap binning results show very low, or even no, redundancy. However, single-sample, batch binning 
    # results contain substantial redundancy, necessitating the use of drep to obtain non-redundant genomes.
    # 挑单菌测序的基因组存在大量冗余。metawrap混合分箱的结果中冗余度非常低，甚至无冗余。
    # 而单样本、分批次分箱的结果中存在大量冗余，需要采用drep获得非冗余的基因组。
    s=drep
    cd ${db} 
    
    ### Opt 1. Download install (方法1.下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack

    ### Opt 2. Conda install (方法2. conda安装)
    # 2025/11/4 install dRep v3.6.2
    conda create -y -n drep
    conda activate drep
    pip install drep -i https://pypi.tuna.tsinghua.edu.cn/simple
    conda install -c bioconda numpy matplotlib pysam -y
    conda install -c bioconda hmmer prodigal pplacer -y
    pip3 install checkm-genome -i https://pypi.tuna.tsinghua.edu.cn/simple
    pip install mash 
    pip install fastANI 
    pip install networkx
    conda pack -f --ignore-missing-files -n ${s} -o ${s}.tar.gz
    
    # checkm_data database set 
    echo "export CHECKM_DATA_PATH=$db/checkm/" >> ~/.bashrc
    # drep & checkm-genome version 
    dRep -h  # v3.5.0
    checkm -h # v1.2.3
    dRep check_dependencies # mash, checkm, prodigal,fastANI all good

## 4.3 CoverM quantify in genome (基因组定量)

    s=coverm
    cd ${db} 
    
    ### Opt 1. Download install (方法1.下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack
    
    ### Opt 2. Conda install & package (方法2. conda安装)
    mamba create -n coverm -c bioconda coverm -y
    conda activate coverm
    conda pack -f --ignore-missing-files -n coverm -o coverm.tar.gz

    # recored version
    coverm -V # 0.7.0

## 4.4 GTDB taxonomic classifications of prokaryote (原核基因组注释)

    # Home: https://gtdb.ecogenomic.org/
    # Github: https://github.com/Ecogenomics/GTDBTk
    s=gtdbtk

### GTDB-tk software install (软件安装)

    ### Opt 1. Download install (方法1.下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack

    ### Opt 2. Conda install & package (方法2. conda安装)
    # gtdbtk-2.5.2, 2025-11-17
    n=gtdbtk
    mamba create -y -n ${s} -c conda-forge -c bioconda gtdbtk=2.5.2
    conda activate ${s}
    # conda pack, --exclude database
    conda pack -n ${s} -o ${s}.tar.gz --exclude gtdbtk-2.5.2 --ignore-editable-packages --ignore-missing-files
    chmod 755 *

    # recored version
    gtdbtk -v # 2.5.2

### GTDB-tk database install (软件安装)

    ### Opt 1. Download install (方法1.下载解压安装)
    mkdir -p ${db}/gtdb && cd ${db}/gtdb
    wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    # Backup Download options include NMDC, Baidu NetDisk conda, etc (备用下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/meta/gtdbtk/gtdbtk_data.tar.gz
    # 手工解压，指定安装完整路径
    tar xvzf gtdbtk_data.tar.gz -C ./  --strip 1
    conda env config vars set GTDBTK_DATA_PATH="${db}/gtdb"
    
    ### Opt 2. Script install & package (方法2. conda安装)
    # download-db.sh, modify database location, wget -c for network continue
    sed -i 's#miniconda3/envs/gtdbtk2.5/share/gtdbtk-2.5.2/db#db/gtdb#;s/wget /wget -c /' ${soft}/envs/gtdbtk/bin/download-db.sh
    # Download database 132 G
    download-db.sh
    
## 4.5 CheckM2 genome accession

    # GitHub：https://github.com/chklovski/CheckM2
    # Conda：https://bioconda.github.io/recipes/checkm2/README.html
    s=checkm2

    ### Opt 1. Conda install & package (方法1. conda安装)
    mamba create --name ${s} ${s} -y
    conda activate ${s}
    ${s} -h # CheckM2 v1.1.0
    conda pack -f --ignore-missing-files -n ${s} -o ${s}.tar.gz

    ### Opt 2. Download install (方法2.下载解压安装)
    # Download options include NMDC, Baidu NetDisk conda, etc (下载，可选NMDC、百度云等)
    wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    mkdir -p ${soft}/envs/${s}
    tar -xvzf ${s}.tar.gz -C ${soft}/envs/${s}
    conda activate ${s}
    conda unpack

    # Database install (数据库安装)
    mkdir $db/$s
    checkm2 database --download --path $db/$s
    # zenodo backup link: https://zenodo.org/records/14897628 diamond 2.1.11 (2025/2/20)  Download need VPN
    
    # Backup database
    wget -c ftp://download.nmdc.cn/tools/meta/checkm2/uniref100.KO.1.dmnd.tar.gz
    tar xvzf uniref100.KO.1.dmnd.tar.gz
    
    #指定数据库
    export CHECKM2DB="/data/meta/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd"  
    # 测试
    checkm2 testrun
    
    # 运行，输入目录或文件列表
    checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder> 
    checkm2 predict --threads 30 --input ../bin1.fa ../../bin2.fna /some/other/directory/bin3.fasta --output-directory <output_folder> 
    
    #打包
    cd ~/project/EasyMetagenome/package
    n=checkm2
    conda pack -f --ignore-missing-files -n ${n} -o ${n}.tar.gz
    cd ..



# 5. Genome and virome (5. 单菌基因组、病毒组等其他软件)


# 泛基因组鉴定软件anvio-8安装

安装 https://anvio.org/install/linux/stable/

    # 软件安装
    #anvio安装 https://anvio.org/install/linux/stable/
    conda create -y --name anvio-8 python=3.10
    conda activate anvio-8

    # 安装依赖
    mamba install -y -c conda-forge -c bioconda python=3.10 \
            sqlite prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
            blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
            trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
            r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript \
            nodejs fastANI

    # install anvi'o
    curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz \
            --output temp/anvio-8.tar.gz
    mv anvio-8.tar.gz anvio8.tar.gz
    # 安装 
    wget -c 
    pip install temp/anvio8.tar.gz -i https://pypi.tuna.tsinghua.edu.cn/simple  #使用清华源加速安装
    #查看是否安装成功，弹出帮助页面即成功
    anvi-merge -h  

    #打包，方便安装
    cd ~/project/EasyMetagenome/package
    n=anvio-8
    conda pack -f --ignore-missing-files -n ${n} -o ${n}.tar.gz
    cd ..

    #安装
    cd ~/project/EasyMetagenome/package
    n=anvio-8
    tar -xvzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    conda create --name ${n} --clone ~/miniconda3/envs/${n}
    conda activate ${n}



# 常见问题

## 软件和数据库国内备份

### 国家微生物科学数据中心 —— 数据下载

http://nmdc.cn/datadownload，可以使用Filezilla直接连接 ftp://download.nmdc.cn/tools

本资源由宏基因组平台发起，微生物所提供服务器，宏基因组团队负责维护的常用软件、扩增子和宏基因组数据库的国内下载链接，解决常用数据库下载慢、或无法下载的问题。同时提供定制的软件、数据库索引，节约大家下载时间，节省数据库编制索引的计算资源消耗。

    # humann3为例
    mkdir -p ~/db/humann3 && cd ~/db/humann3
    site=ftp://download.nmdc.cn/tools
    wget -c ${site}/humann3/full_chocophlan.v296_201901.tar.gz
    wget -c ${site}/humann3/uniref90_annotated_v201901.tar.gz
    wget -c ${site}/humann3/full_mapping_v201901.tar.gz
    mkdir -p chocophlan uniref utility_mapping
    tar xvzf full_chocophlan.v296_201901.tar.gz -C chocophlan/
    tar xvzf uniref90_annotated_v201901.tar.gz -C uniref/
    tar xvzf full_mapping_v201901.tar.gz -C utility_mapping/

### 百度云备份链接

https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315

下载的tar.gz压缩包，可放置于指定目录，使用`tar -xvzf *.tar.gz`解压 

    # 大文件的分卷压缩和解压 以kraken2为例
    cd ~/db/kraken2
    # https://www.cnblogs.com/wang--lei/p/9046643.html
    # 文件夹kraken2/打包压缩，1h
    tar -zcvf kraken2.tar.gz kraken2/
    # b分割为指定大小文件G/M/K，-d数字，a序列长度，输入和输出前缀
    split -b 13G -d -a 1 kraken2.tar.gz kraken2.tar.gz.
    # 一行命令打包并分割
    tar -zcvf kraken2.tar.gz kraken2idx/ | split -b 19G -d -a 1 - kraken2.tar.gz.
    # 分割后合并及解压缩
    cat kraken2.tar.gz.* | tar -zxv
    
    # 附tar打包用法，c创建、v输出过程，z压缩，f文件 ，x解压
    单个文件压缩打包 tar -cvzf my.tar.gz file1
    多个文件压缩打包 tar -cvzf my.tar.gz file1 file2 file*
    单个目录压缩打包 tar -cvzf my.tar.gz dir1
    多个目录压缩打包 tar -cvzf my.tar.gz dir1 dir2
    解压缩至当前目录：tar -xvzf my.tar.gz

    
### Lefse在Rstudio中运行命令调用R版本问题的解决

    # 在Rstudio中默认调用Rstudio的R，具体写在/etc/rstudio/rserver.conf
    # 或在R中用Sys.getenv()["R_HOME"]，在rpy2中print(robjects.r)可以查看其调用的r版本
    # 指定lefse调用的R版本，需根据conda实际目录修改
    sed -i "2 i os.environ['R_HOME'] = '~/miniconda3/envs/meta/lib/R/'" \
      ~/miniconda3/envs/meta/share/lefse-1.0.8.post1-1/lefse.py
      
### Perl版本不对

常见问题：Perl版本不对，人工指定perl版本如下

    PERL5LIB=~/miniconda3/envs/kraken2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:~/miniconda3/envs/kraken2/lib/site_perl/5.26.2:~/miniconda3/envs/kraken2/lib/5.26.2/x86_64-linux-thread-multi:~/miniconda3/envs/kraken2/lib/5.26.2


## salmon手动安装和使用

    # 如不可用，尝试下载二进制和添加环境变量
    wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.0/salmon-0.14.0_linux_x86_64.tar.gz
    tar xvzf salmon-0.14.0_linux_x86_64.tar.gz 
    cp -rf salmon-latest_linux_x86_64/ ${soft}/envs/metagenome_env/share/salmon
    # 或者直接使用软件全路径
    ${soft}/envs/metagenome_env/share/salmon/bin/salmon -v # 0.14.0


# 附录

## Conda安装小工具

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
    
## 宿主参考基因组下载

- EnsembleGenomes http://ensemblgenomes.org/ 
- 包括动物、植物、原生生物、真菌、细菌等，此外植物还 Phytozome https://phytozome-next.jgi.doe.gov/ ，以及单个物种和专用数据库

以Ensemble中拟南芥为例：Arabidopsis thaliana -- Genome assembly -- Download DNA sequence (无反应)，点TAIR链接跳转ENA，下载All Seq FASTA
    
    wget https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_000001735.1?download=true&gzip=true
    mv GCA_000001735.1\?download\=true TAIR10.fa

以Ensemble中水稻为例：Oryza sativa Japonica —— IRGSP-1.0
    
    wget https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001433935.1?download=true&gzip=true
    mv GCA_001433935.1\?download\=true IRGSP1.0.fa

## KEGG层级注释整理

己整合至EasyMicrobiome中，自己更新请访问 https://www.kegg.jp/kegg-bin/show_brite?ko00001.keg 下载htext

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

## 毒力因子数据库VFDB

官网：http://www.mgc.ac.cn/VFs/ 数据每周更新
 
    mkdir -p ${db}/vfdb && cd ${db}/vfdb
    # 毒力因子描述文件
    wget -c http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
    # 核心数据库(1117K)
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
    # 完整数据库(4.99M)
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
    # 解压
    gunzip *.gz

### Conda usuage

    # Remove env (删除环境)
    conda env remove --name drep
    conda env remove --name qiime2-2023.7


# Change log (版本更新记录)

    # **1.24 2025.11.6**
    # 1. update fastp 0.23.4 to 1.0.1
    # 2. taxonkit v0.14.1 to v0.20.0
    # 3. dbcan3 2025 data format update format_dbcan2list.pl to format_dbcan3list.pl
    # 4. rgi 6.0.3 to 6.0.5
    # **正在开发中功能**
    # 1.  rgi应用于菌群分析及结果展示
    # 2.  antisamsh应用于菌群分析及结果展示
    # 3.  cazy应用于菌群分析及结果展示