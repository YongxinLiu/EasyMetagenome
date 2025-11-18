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

## MetwWRAP binning (分箱)

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

## drep基因组去冗余

挑单菌测序的基因组存在大量冗余。metawrap混合分箱的结果中冗余度非常低，甚至无冗余。而单样本、分批次分箱的结果中存在大量冗余，需要采用drep获得非冗余的基因组。
GitHub: https://github.com/MrOlm/drep
Conda: https://bioconda.github.io/recipes/drep/README.html

### drep 基因组去冗余解包安装

    # 下载dRep v3.5.0
    wget -c ftp://download.nmdc.cn/tools/conda/drep.tar.gz
    # 指定安装目录
    mkdir -p ${soft}/envs/drep
    tar -xvzf drep.tar.gz -C ${soft}/envs/drep
    # 启动环境
    conda activate drep
    # 初始化环境
    conda unpack
    dRep -h

### drep 基因组去冗余直接安装

    # 2024/11/13安装drep最新版 v3.5.0
    conda create -y -n drep
    conda activate drep
    
    #使用清华源加速安装
    pip install drep==3.5.0  -i https://pypi.tuna.tsinghua.edu.cn/simple
    conda install -c bioconda numpy matplotlib pysam -y
    conda install -c bioconda hmmer prodigal pplacer -y
    pip3 install checkm-genome -i https://pypi.tuna.tsinghua.edu.cn/simple
    mamba install -c bioconda mash fastANI -y
    
    
    #checkm_data数据库下载
    #checkm_data数据从https://zenodo.org/record/7401545/files/checkm_data_2015_01_16.tar.gz（需vpn）或https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz网址下载，然后放到~/miniconda3/envs/drep/checkm_data/中解压
    
    cd ~/miniconda3/envs/drep/checkm_data/ 
    wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xvf checkm_data_2015_01_16.tar.gz
    
    #指定checkm数据库位置
    export CHECKM_DATA_PATH=~/miniconda3/envs/drep/checkm_data/
    
    # 查看checkm-genome 和 drep版本
    checkm -h # v1.2.3
    dRep -h  # v3.5.0
    dRep check_dependencies #查看依赖项安装情况，确保mash、checkm、prodigal、fastANI显示all good
    
    #测试dRep
    conda activate drep
    mkdir -p ~/project/EasyMetagenome/temp/drep_in
    cd ~/project/EasyMetagenome/temp/drep_in
    cp /data1/liuyongxin/age/meta/temp/drep_in/Gp_EAB* ./
    cd ~/project/EasyMetagenome
    mkdir -p temp/drep95
    time dRep dereplicate temp/drep95/ \
        -g temp/drep_in/*.fa  \
        -sa 0.95 -nc 0.30 -comp 50 -con 10 -p 5   #log文件内显示 Finished the dereplicate operation!
    
    #打包
    cd ~/project/EasyMetagenome/package
    n=drep
    conda pack -f --ignore-missing-files -n ${n} -o ${n}.tar.gz
    cd ..

### drep 数据库构建

CheckM用于Bin完整和污染估计和物种注释，安装过metawrap已经下载完成

    mkdir -p ${db}/checkm && cd ${db}/checkm
    # 下载文件275 MB，解压后1.4 GB
    wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xvf *.tar.gz
    # 设置数据库位置，直接2次回车默认为当前位置
    checkm data setRoot `pwd`

## coverm基因组定量

conda安装

    mamba create -n coverm -c bioconda coverm -y
    conda activate coverm
    coverm -V #0.7.0
    # conda安装后打包(可选)
    conda pack -f --ignore-missing-files -n coverm -o coverm.tar.gz

压缩包安装

    # 指定conda文件名
    s=coverm
    # 下载，可选NMDC、百度云等
    # wget -c ftp://download.nmdc.cn/tools/conda/${s}.tar.gz
    # 指定安装目录
    mkdir -p ~/miniconda3/envs/${s}
    tar -xvzf ${s}.tar.gz -C ~/miniconda3/envs/${s}
    # 启动环境
    conda activate ${s}
    # 初始化环境
    conda unpack

## GTDB细菌基因组注释和进化分析

Github: https://github.com/Ecogenomics/GTDBTk
GTDB-Tk是一个软件工具包，用于根据基因组数据库分类法GTDB为细菌和古细菌基因组分配客观的分类法。它旨在与最近的进展一起使用，从而可以直接对环境样本中获得数百或数千个由宏基因组组装的基因组（MAG）进行物种分类注释。它也可以用于分离和单细胞的基因组物种注释。
本次测试版本为 gtdbtk-2.4.0，Release 09-RS220 (24th April 2024)。
硬件要求：内存>256Gb，硬盘>256Gb，64核1小时可分析1000个细菌基因组

### GTDB-Tk直接安装

    # gtdbtk-2.4.0, 2024-5-5
    n=gtdbtk2.4
    mamba create -y -n ${n} -c conda-forge -c bioconda gtdbtk=2.4.0
    # 检查版本
    conda activate ${n}
    gtdbtk -v # 2.4.0
    
    # conda pack软件打包一次, --exclude排除数据库
    conda pack -n ${n} -o ${n}.tar.gz --exclude gtdbtk-2.4.0 --ignore-editable-packages --ignore-missing-files
    chmod 755 *
    
### GTDB-Tk解包安装

    soft=~/miniconda3
    n=gtdbtk2.4
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # 指定安装目录
    mkdir -p ${soft}/envs/${n}
    tar -xvzf ${n}.tar.gz -C ${soft}/envs/${n}
    # 启动环境
    conda activate ${n}
    # 初始化环境
    conda unpack
    
### GTDB-Tks数据库安装

download-db.sh自动下载数据库，将下载至conda中的envs/gtdbtk/share/gtdbtk-2.4.0/db/，我们修改为~/db/gtdb中

    conda activate gtdbtk2.4
    # download-db.sh中，修改数据库下载位置，的 wget 建议改成wget -c 防止覆盖
    sed -i 's#miniconda3/envs/gtdbtk2.4/share/gtdbtk-2.4.0/db#db/gtdb2.4#;s/wget /wget -c /' ${soft}/envs/gtdbtk2.4/bin/download-db.sh
    # 下载数据,101G
    download-db.sh
    
(备选)上面无法下载时手动下载和配置GTDB数据库

    mkdir -p ${db}/gtdb2.4 && cd ${db}/gtdb2.4
    # 下载
    wget -c https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    
    # 备用链接
    wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    # 手工解压，指定安装完整路径
    tar xvzf gtdbtk_data.tar.gz -C ./  --strip 1
    conda env config vars set GTDBTK_DATA_PATH="${db}/gtdb/2404"

# 5 单菌基因组、病毒组等其他软件

## CheckM2

[metawrap下载安装](http://127.0.0.1:3334/milkdown/#4QPx-1732175063664)

Conda主页：https://bioconda.github.io/recipes/checkm2/README.html

软件主页：https://github.com/chklovski/CheckM2

    # 软件安装
    mamba create --name checkm2 checkm2 -y
    conda activate checkm2
    checkm2 -h # CheckM2 v1.0.2
    # 数据库安装
    mkdir ~/db/checkm2
    checkm2 database --download --path ~/db/checkm2
    # 报错：checkm2.zenodo_backpack.ZenodoConnectionException: Connection error: HTTPSConnectionPool(host='zenodo.org', port=443): Max retries exceeded with url: /record/5571251 (Caused by NewConnectionError('<urllib3.connection.HTTPSConnection object at 0x7f335c64abe0>: Failed to establish a new connection: [Errno 111] Connection refused'))
    # 在github中搜索和查找issues无解答，提新issue https://github.com/chklovski/CheckM2/issues
    # 数据库 https://zenodo.org/records/5571251  下载，需要VPN
    
    # 备用数据库下载并解压(待上传)
    wget -c ftp://download.nmdc.cn/tools/meta/CheckM2_database/uniref100.KO.1.dmnd.tar.gz
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

    
## kneaddata常见问题

### kneaddata运行提示java版本不支持
    
	# 解决思路，新建虚拟环境，安装kneaddata，再安装对应java版本
	# 务必指定2.7，软件依赖2.7的python，但conda会自动安装3.6，运行报错
    conda create -n kneaddata python=2.7
    conda activate kneaddata
    conda install openjdk=8.0.152
    conda install kneaddata=0.6.1

### 解压失败-重新下载再安装

    tar -xvzf kneaddata.tar.gz -C ~/miniconda3/envs/kneaddata

解压文件提示如下错误

    gzip: stdin: invalid compressed data--format violated
    tar: Unexpected EOF in archive
    tar: Unexpected EOF in archive
    tar: Error is not recoverable: exiting now

检查md5值确认文件是否不同

    md5sum kneaddata.tar.gz

当前为d26125bee1def1faa99d03a9715bf392
原文件为9fa47a364096b2c33be52a91850b2cde

删除当前文件并重新下载即可

    rm kneaddata.tar.gz
    wget ftp://download.nmdc.cn/tools//conda/kneaddata.tar.gz

### Lefse在Rstudio中运行命令调用R版本问题的解决

    # 在Rstudio中默认调用Rstudio的R，具体写在/etc/rstudio/rserver.conf
    # 或在R中用Sys.getenv()["R_HOME"]，在rpy2中print(robjects.r)可以查看其调用的r版本
    # 指定lefse调用的R版本，需根据conda实际目录修改
    sed -i "2 i os.environ['R_HOME'] = '~/miniconda3/envs/meta/lib/R/'" \
      ~/miniconda3/envs/meta/share/lefse-1.0.8.post1-1/lefse.py
      
## Kraken2

### 定制数据库

官方教程详见 https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

本地构建最完整索引，自定义微生物数据库，如标准+真菌+原生动物+质粒+植物

    mkdir -p ${db}/kraken2/kraken2_self
    conda activate kraken2
    # 显示帮助
    kraken2-build -h
    # 下载物种注释
    kraken2-build --download-taxonomy --threads 24 --db ${db}/kraken2/kraken2_self
    # 下载数据库，需要12-24小时
    for i in archaea bacteria UniVec_Core viral human fungi plasmid protozoa plant; do
        kraken2-build --download-library $i --threads 24 --db ${db}/kraken2/kraken2_self
    done
    # 确定的库建索引，4p,4h
    time kraken2-build --build --threads 48 --db ${db}/kraken2/kraken2_self
    # bracken索引，长度推荐100/150, 24p, 1h;
    time bracken-build -d ./ -t 24 -k 35 -l 100
    time bracken-build -d ./ -t 24 -k 35 -l 150

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

## MetaWRAP分箱

### shorten_contig_names.py报错

更新 ${soft}/envs/metawrap/bin/metawrap-scripts/shorten_contig_names.py 脚本

    #!/usr/bin/env python2.7
    import sys
    shorten=False
    for line in open(sys.argv[1]):
    	if line[0]!=">":
    		print line.rstrip()
    	else:
    		if shorten==True:
    			#print "_".join(line.rstrip().split("_")[:4])
    			lineL = line.rstrip().split("_")
    			new_line = '_'.join([lineL[0], lineL[1], lineL[3]])
    			print new_line[:20]
    		elif len(line)>20 and len(line.split("_"))>5:
    			lineL = line.rstrip().split("_")
    			new_line = '_'.join([lineL[0], lineL[1], lineL[3]])
    			#print "_".join(line.rstrip().split("_")[:4])
    			print new_line[:20]
    			shorten=True
    		else:
    			print line.rstrip()

### 绘图plot_binning_results.py报错

更新 ${soft}/envs/metawrap/bin/metawrap-scripts/plot_binning_results.py 脚本

    # 原脚本存在嵌套错误，会输出报错，修改部分如下
    # Traceback (most recent call last):
    #  File "/anaconda3/envs/metawrap-env/bin/metawrap-scripts/plot_binning_results.py", line 119, in <module>
    #   plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)
    # NameError: name 'x_pos' is not defined

        # add bin set label to plot
        for x_pos,y in enumerate(data[bin_set]):
                if y>y_pos:
                        break
                plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)
                y_pos+=y_increment
                
    # add plot and axis titles and adjust the edges
    plt.title("Bin contamination ranking", fontsize=26)
    plt.xlabel("Acending contamination rank", fontsize=16)
    plt.ylabel("Estimated bin contamination (log scale)", fontsize=16)
    plt.gcf().subplots_adjust(right=0.9)
    
    # save figure
    print "Saving figures binning_results.eps and binning_results.png ..."
    plt.tight_layout(w_pad=10)
    plt.subplots_adjust(top=0.92, right=0.90, left=0.08)
    plt.savefig("binning_results.png",format='png', dpi=300)
    plt.savefig("binning_results.eps",format='eps')
    #plt.show()
    EOF

### blast版本不兼容

更新 metawrap 中的 blast 版本，直接到https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/，下载最新版本blastn，再到conda的metawrap环境bin目录下，替换掉旧版本的blastn

    # 如果出现这个错误，BLAST Database error: Error: Not a valid version 4 database.
    # 是metawrap 中 blast 版本太老了，需要更新下
    wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz
    tar xvzf ncbi-blast-2.13.0+-x64-linux.tar.gz
    mv ncbi-blast-2.13.0+/bin/* ${soft}/envs/metawrap-env/bin/


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

## 宏基因组基于读长的分析 HUMAnN2/Metaphlan2/graphlan

### HUMAnN2解包安装

    # 下载
    wget -c ftp://download.nmdc.cn/tools/conda/humann2.tar.gz
    # 指定安装目录
    mkdir -p ${soft}/envs/humann2
    tar -xvzf humann2.tar.gz -C ${soft}/envs/humann2
    # 启动环境
    conda activate humann2
    # 初始化环境
    conda unpack

### HUMAnN2直接安装

    # mamba 是快速版本的 conda
    mamba create -n humann2 humann2 graphlan export2graphlan -c bioconda -y

### HUMAnN2安装测试

    conda activate humann2
    # 记录核心软件版本
    humann2 --version # v2.8.1
    metaphlan2.py -v # 2.7.5 (6 February 2018)
    diamond help | head -n 1 #  v0.8.36.98
    graphlan.py --version # 1.1.3 (5 June 2018)
    export2graphlan.py -h # 0.22 of 05 May

    # 测试流程是否可用
    humann2_test

### HUMAnN2物种和功能数据库

    # 显示可用分类、泛基因组和功能数据库
    humann2_databases

    # 安装数据库(注：数据库下载慢或失败，附录有国内备份链接)
    cd ${db}
    mkdir -p ${db}/humann2 # 建立下载目录
    # 输助比对数据库 593MB
    humann2_databases --download utility_mapping full ${db}/humann2
    # 微生物泛基因组 5.37 GB
    humann2_databases --download chocophlan full ${db}/humann2
    # 功能基因diamond索引 10.3 GB
    humann2_databases --download uniref uniref90_diamond ${db}/humann2

    # humann2数据库无法下载：附录备用链接下载后手动配置
    mkdir -p ${db}/humann2/chocophlan && cd ${db}/humann2/chocophlan
    tar xvzf full_chocophlan_plus_viral.v0.1.1.tar.gz
    mkdir -p ${db}/humann2/uniref && cd ${db}/humann2/uniref
    tar xvzf uniref90_annotated_1_1.tar.gz
    mkdir -p ${db}/humann2/utility_mapping && cd ${db}/humann2/utility_mapping
    tar xvzf full_mapping_1_1.tar.gz


    # 设置数据库位置
    # 显示参数
    humann2_config --print
    # 如修改线程数，推荐3-8，根据实际情况调整
    humann2_config --update run_modes threads 4
    humann2_config --update database_folders utility_mapping ${db}/humann2/utility_mapping
    humann2_config --update database_folders nucleotide ${db}/humann2/chocophlan
    humann2_config --update database_folders protein ${db}/humann2/uniref
    humann2_config --print
    
    ## metaphlan2数据库下载和配置
    mkdir -p ${db}/humann2 && cd ${db}/humann2
    wget -c ftp://download.nmdc.cn/tools/humann2/metaphlan2.tar.gz
    tar xvzf metaphlan2.tar.gz
    # 链接到软件安装目录
    mkdir -p ${soft}/envs/humann2/bin/databases
    ln -s ${db}/humann2/metaphlan2/* ${soft}/envs/humann2/bin/databases/

    n=kneaddata
    conda pack -f --ignore-missing-files -n ${n} -o ${n}.tar.gz


HUMAnN3+MetaPhlAn4为目前最新版，目前最广泛使用的HUMAnN2安装见附录

### HUMAnN3直接安装

    #1.更新--------------------
    conda create -n humann3
    conda activate humann3
    conda install metaphlan=4.1.1 humann=3.9  -c bioconda -c conda-forge -y
    humann --version # v3.9
    metaphlan -v # version 4.1.1 (11 Mar 2024)
    
    #2.测试--------------------
    # 测试
    humann_test #无报错最后显示 OK
    
    #3.打包--------------------
    #安装软件打包，f覆盖输出文件，ignore跳过修改检测
    cd package
    n=humann3
    conda pack -f --ignore-missing-files -n ${n} -o ${n}.tar.gz
    cd ..
    
### HUMAnN3解包安装

    # 下载
    wget -c ftp://download.nmdc.cn/tools/conda/humann3.tar.gz
    # 指定安装目录
    mkdir -p ${soft}/envs/humann3
    tar -xvzf humann3.tar.gz -C ${soft}/envs/humann3
    # 启动环境
    conda activate humann3
    # 初始化环境
    conda unpack

### HUMAnN3安装测试

    # 记录核心软件版本
    humann --version # v3.7
    metaphlan -v # 4.0.6 (1 Mar 2023)
    diamond help | head -n 1 #  v2.1.8.162
    # 测试
    humann_test

### HUMAnN3物种和功能数据库

    # 显示可用分类、泛基因组和功能数据库
    humann_databases
    
    # 安装数据库
    cd ${db}
    mkdir -p ${db}/humann3 # 建立下载目录
    # 微生物泛基因组 16 GB
    humann_databases --download chocophlan full ${db}/humann3
    # 功能基因diamond索引 20 GB
    humann_databases --download uniref uniref90_diamond ${db}/humann3
    # 输助比对数据库 2.6 GB
    humann_databases --download utility_mapping full ${db}/humann3
    
    # humann3数据库无法自动下载，备用链接下载安装
    wget -c ftp://download.nmdc.cn/tools/meta/humann3/full_chocophlan.v201901_v31.tar.gz
    wget -c ftp://download.nmdc.cn/tools/meta/humann3/uniref90_annotated_v201901b_full.tar.gz
    wget -c ftp://download.nmdc.cn/tools/meta/humann3/full_mapping_v201901b.tar.gz
    # 安装、解压
    mkdir -p ${db}/humann3/chocophlan
    tar xvzf full_chocophlan.v201901_v31.tar.gz -C ${db}/humann3/chocophlan
    mkdir -p ${db}/humann3/uniref
    tar xvzf uniref90_annotated_v201901b_full.tar.gz -C ${db}/humann3/uniref
    mkdir -p ${db}/humann3/utility_mapping			66	
    tar xvzf full_mapping_v201901b.tar.gz -C ${db}/humann3/utility_mapping
    
    # 设置数据库位置
    # 显示参数
    humann_config --print
    # 如修改线程数，推荐3-8，根据实际情况调整
    humann_config --update run_modes threads 8
    # 设置核酸、蛋白和注释库位置
    humann_config --update database_folders nucleotide ${db}/humann3/chocophlan
    humann_config --update database_folders protein ${db}/humann3/uniref
    humann_config --update database_folders utility_mapping ${db}/humann3/utility_mapping
    # 核对设置结果
    humann_config --print

### MetaPhlAn4物种数据库
    
    # MetaPhlAn4数据库下载2024数据和索引2.98G+19.87G
    cd ~/project/EasyMetagenome
    mkdir -p ~/db/metaphlan4
    cd ~/db/metaphlan4
    
    # 官网下载
    wget -c http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202403.tar
    wget -c http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar
    tar xvf mpa_vOct22_CHOCOPhlAnSGB_202403.tar
    tar xvf mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar
    
    # 官方没有压缩体积大下载慢，备用国内百度链接：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315 或 微生物所FTP ftp://download.nmdc.cn/tools/meta 下载压缩包
    wget -c ftp://download.nmdc.cn/tools/meta/metaphlan4/mpa_vOct22_CHOCOPhlAnSGB_202403.tar.gz
    wget -c ftp://download.nmdc.cn/tools/meta/metaphlan4/mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar.gz
    tar xvzf mpa_vOct22_CHOCOPhlAnSGB_202403.tar.gz
    tar xvzf mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar.gz
    # 可选(制作下载文件和md5值)
    gunzip mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar.gz
    md5sum mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.tar > mpa_vOct22_CHOCOPhlAnSGB_202403_bt2.md5


### Conda用法补充

    # 删除env环境
    conda env remove --name qiime2-2023.7
    conda env remove --name picrust2

    # 仓库优先级
    conda config --set channel_priority strict # 设置严格的仓库优先级（最好不要使用）
    conda config --set channel_priority flexible # 禁用仓库优先级


# 版本更新记录

**1.08 2020.7.20**

1.  KneadData提供数据预处理双端标签唯一命令，兼容最新版；
2.  提供HUMAnN3测试版的安装和分析流程(附录1)；
3.  eggNOG升级为emapper 2.0和eggNOG 5.0流程，结果列表从13列变为22列，新增CAZy注释。emapper 1.0版本见附录2。

**1.09 2020.10.16**

1.  新增二、三代混合组装OPERA-MS软件使用 (31Megahit)
2.  新增eggNOG-mapper结果COG/KO/CAZy整理脚本summarizeAbundance.py，删除旧版Shell+R代码 (32Annotation)
3.  新增MetaWRAP单样本分箱流程 (33Binning)
4.  新增dRep实现基因组去冗余 (34Genomes)
5.  新增GTDB-Tk基因组物种注释和进化树构建 (34Genomes)

**1.10 2021.1.22**

1.  增加删除中间文件部分，节约空间，防止硬盘写满；
2.  正文的补充分析方法、常见问题移至附录，按软件名、问题/方法分级索引；
3.  软件使用前，增加检查软件版本命令，方便文章方法中撰写准确版本；
4.  删除不稳定的humann3、过时的eggnog版本教程；
5.  增加kraken2新环境, 增加bracken, krakentools新工具；
6.  kraken2结果新增beta多样性PCoA，物种组成堆叠柱状图；
7.  增metaspades二、三代组装代码示例；
8.  新增KEGG层级注释整理代码；
9.  更新dbcan3中2018版为2020版；
10. 新增CARD本地分析流程；

**1.11 2021.5.7**

1.  增加prodigal基因预测并行版方法，使用seqkit split拆分后并行，数10倍加速单线程基因预测步骤；
2.  增加megahit拼装结果片段大小选择步骤，使用seqkit -m按长度筛选，并统计筛选前后变化；
3.  不常用或可选代码调整到附录
4.  两批数据快速合并去冗余cd-hit-est-2d
5.  二三代混合组装OPERA-MS的混装和3代优化代码

**1.12 2021.8.20**

1.  新增并行管理软件rush，比parallel更易安装，绿色版无依赖关系，整合在db/linux/目录中
2.  新增seqkit，可以统计序列数据量，支持序列长度过滤，格式转换等；
3.  新增质控软件fastp，软件fastqc更快，适合单独质控不去宿主；
4.  kraken2新数据库，同样大小下注释率提高明显；
5.  eggNOG软件和数据库配套升级
6.  GTDB-tk软件和数据库需要配套重新才可使用新版25万基因组数据库

**1.13 2021.11.19**

1.  陈同参与EasyMicrobiome的更新，并提交了mac版本代码
2.  新增humann2运行bowtie2出错的解决方案
3.  新增软件conda环境下载安装方式，且作为首选
4.  新增kneaddata自定义物种基因组数据库示例

**1.14 2022.3.25**

1.  EasyMicrobiome升级为1.14
2.  升级miniconda2为miniconda3
3.  dbcan3从2020/7/31的808M更新为2021/9/24版1016M，格式变化，配套format_dbcan2list.pl更新
4.  新增eggnog环境，包含emapper 2.1.6，summarizeAbundance.py含pandas (conda install sklearn-pandas)，配套更新数据库
5.  rgi更新到最新版及配套代码

**1.15 2022.5.27**

1.  陈同老师全面更新课程，并在新服务器上重新布置所有软件和数据库
2.  课题尝试改为长期：自学理论课程视频，每周线上答疑，持续2个月完成实操

**1.18 2023.4.7**

1.  课程恢复为3天连续学习模式
2.  更新所有软件和数据库为可成功安装的最新版
3.  更新软件和数据备份至微生物所和百度网盘

**1.19 2023.7.13**

1.  HUMAnN2+MetaPhlAn2为3和4
2.  Kraken2数据库不同版本统一官方名称，仍使用3月版本数据库，最新版6月数据库官方文件有不完整
3.  GTDB-tk数据库更新为214版

**1.20 2023.11.24**

1. MetaPhlAn4新增物种注释转换为GTDB、多样性计算脚本
2. 整合陈同老师用Pipeline的修正
3. CoverM定量MAGs相对丰度、结果合并和求均值，并添加到进入树注释结果中
4. drep的conda包为3.4.2缺少checkm(267M)，替换为旧版2.6.2(526M)
5. dbCAN3数据库更新为2023版，diamond新版建索引更快
6. Kneaddata质控跳过，fastp质控为必选步骤
7. mutliqc升级1.14为1.15
8. 增加第五章：单菌基因组分析流程
9. 更新Kraken2数据库为20231009版本，新增alpha, beta多样性、Krona网页、Pavian桑基图
10. 新增可选的checkm2评估

**1.21 2024.5.4**
1. format_dbcan2list.pl更新，解决结果丢失第一列结果的bug，新增了Evalue，及按Evalue筛选的参数
2. 新增viwarp软件数据库：iPHoP.latest_rw.tar.gz(116G)、METABOLIC_test_files.tgz(2G)
3. kraken2数据库更新为2024版，kraken2从2.1.2升级为2.1.3，环境名为kraken2.1.3，且bracken2.5升级为2.8，解决结果校正后大量为0的Bug；
4. 测试数据新增不同疾病程度、不同年龄组；
5. CARD更新为2024版，v3.2.9

**1.22 2024.11.9**
1. PPT更新为2024.11版
2. kraken2/tax_count.spf 格式更新，解决缺失一个样本的问题；
3. CAZy数据库更新为2024版


**1.23 2024.11.18**
1.增加了4.5功能注释-耐药基因的处理步骤
2.增加了6(可选)泛基因组分析部分

**1.24 2025.11.6**
1. update fastp 0.23.4 to 1.0.1
2. taxonkit v0.14.1 to v0.20.0
3. dbcan3 2025 data format update format_dbcan2list.pl to format_dbcan3list.pl
4. rgi 6.0.3 to 6.0.5

**正在开发中功能**
1.  rgi应用于菌群分析及结果展示
2.  antisamsh应用于菌群分析及结果展示
3.  cazy应用于菌群分析及结果展示