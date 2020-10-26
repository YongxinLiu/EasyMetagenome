[TOC]

# 宏基因组分析常见问题及补充分析

    # 测试环境为Linux Ubuntu 18.04 / CentOS 7
    # 版本: 1.09, 2020/10/16

## Conda软件和数据库位置

    # 数据库安装位置，默认~/db目录(无需管理权限)，管理员可安装至/db，方便大家使用，安装和运行必备环境变量
    db=~/db
    mkdir -p ${db} && cd ${db}
    soft=~/miniconda2

### miniconda2软件管理器

    wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    # Conda软件安装目录，如管理员可能为/conda，普通用户为~/miniconda2
    bash Miniconda2-latest-Linux-x86_64.sh -b -f -p ${soft}
    # 安装时许可协议可打yes，默认目录为~/miniconda2，默认不运行conda直接回车，这里用-p批定soft变量为目录
    conda -V # 查看版本号 4.6.8
    python --version # 2.7.15
    # 安装说明详见：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)
    # 添加生物频道，才能找到生物学软件
    conda config --add channels bioconda
    # http://bioconda.github.io/ 查询软件及版本
    conda config --add channels conda-forge
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
    # conda默认配置文件为 ~/.condarc ，如不存在使用 conda config --show-sources 查看配置文件位置
    # 添加conda至环境变量
    export PATH="${soft}/bin:$PATH"

    # 创建虚拟环境，防污染环境变量
    conda create -n meta # python=2.7 r-base=3.6
    # 加载环境
    conda activate meta

    # 如果有的软件在Solving environment步骤数小时无法安装，可以新建环境
    conda create -n metaRef
    conda activate metaRef

## 质控和去宿主

### multiqc评估报告汇总multiqc

    # 注：新版本1.9仅支持python3，需要有Python3环境下安装
    conda install multiqc # 111 MB
    multiqc --version # multiqc, version 1.8

## 物种注释

### lefse物种差异比较和绘制

    conda install lefse # 57.5 MB, 1.0.8.post1
    
    # 在Rstudio中默认调用Rstudio的R，具体写在/etc/rstudio/rserver.conf
    # 或在R中用Sys.getenv()["R_HOME"]
    # 在rpy2中 print(robjects.r)可以查看其调用的r版本
    
    # 认为指定R版本
    sed -i "2 i/os.environ['R_HOME'] = '/conda2/envs/metagenome_env/lib/R/'" \
      ${soft}/envs/metagenome_env/share/lefse-1.0.7-1/lefse.py
  
## 组装流程
    
### salmon可选安装

    # 定量工具salmon
    # conda install salmon # 15.1 MB
    wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.0/salmon-0.14.0_linux_x86_64.tar.gz
    tar xvzf salmon-0.14.0_linux_x86_64.tar.gz 
    cp -rf salmon-latest_linux_x86_64/ ${soft}/envs/metagenome_env/share/salmon
    ${soft}/envs/metagenome_env/share/salmon/bin/salmon -v # 0.14.0


### eggNOG2

    # eggNOG http://eggnogdb.embl.de
    # 安装eggnog比对工具emapper
    conda install eggnog-mapper
    emapper.py --version # eggnog-mapper-2.0.1
    
    # 下载常用数据库，注意设置下载位置
    mkdir -p ${db}/eggnog && cd ${db}/eggnog
    # -y默认同意，-f强制下载
    download_eggnog_data.py -y -f --data_dir ./
    # 如果内存够大，复制eggNOG至内存加速比对
    # cp eggnog.db /dev/shm
    # 手工整理COG分类注释
    wget -c wget http://210.75.224.110/share/COG.anno
    # 手工整理KO注释
    wget -c wget http://210.75.224.110/share/KO.anno

### CAZy碳水化合物

    # dbCAN2无法访问 http://cys.bios.niu.edu/dbCAN2/
    mkdir -p ${db}/dbCAN2 && cd ${db}/dbCAN2
    wget -c http://210.75.224.110/share/meta/dbcan2/CAZyDB.07312018.fa # 497 MB
    wget -c http://210.75.224.110/share/meta/dbcan2/CAZyDB.07312018.fam-activities.txt # 58 KB
    time diamond makedb --in CAZyDB.07312018.fa --db CAZyDB.07312018 # 28s
    # 提取fam对应注释
    grep -v '#' CAZyDB.07312018.fam-activities.txt|sed 's/  //'| \
      sed '1 i ID\tDescription' > fam_description.txt

### Resfam耐药基因

    # http://dantaslab.wustl.edu/resfams
    mkdir -p ${db}/resfam && cd ${db}/resfam
    # 官网的数据格式非常混乱, 推荐下载我手动整理的索引和注释
    wget http://210.75.224.110/share/Resfams-proteins.dmnd # 1.5 MB
    wget http://210.75.224.110/share/Resfams-proteins_class.tsv # 304 KB


### Metawarp分箱

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

## 1.6 其它软件和输助脚本

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


    # Bin可视化VizBin
    # sudo apt-get install libatlas3-base libopenblas-base default-jre
    # curl -L https://github.com/claczny/VizBin/blob/master/VizBin-dist.jar?raw=true > VizBin-dist.jar
    # mv VizBin-dist.jar /usr/local/bin # 或~/bin
    
    # 比对结果整理samtools
    conda install samtools
    
    ### CARD(选学) https://card.mcmaster.ca/download 
    # 方法1. 直接conda安装
    conda install --channel bioconda rgi
    # 如果报错，尝试方法2。CondaMultiError: CondaFileIOError: '/home/liuyongxin/miniconda2/pkgs/prokka-1.11-0.tar.bz2'. contains unsafe path: db/cm/READM
    # 方法2. 虚拟环境安装
    conda activate metawrap
    conda install --channel bioconda rgi
    rgi main -v # 4.0.3
    # rgi教程 https://github.com/arpcard/rgi

## 附录

### conda批量安装软件

    conda install -y fastqc multiqc kneaddata=0.6.1 humann2 graphlan export2graphlan lefse kraken2 megahit spades quast prokka cd-hit emboss salmon eggnog-mapper samtools
    
### kneaddata质控双端结果不成对
    
    # 0.7.2存在对旧格式fastq去宿主后数据极少或双端数据不对称，可指定版本安装0.6.1
    conda install kneaddata=0.6.1 # 175 MB
    
### kneaddata运行提示java版本不至
    
	# 解决思路，新建虚拟环境，安装kneaddata，再安装对应java版本
    conda install openjdk=8.0.152
    	
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

### kraken2定制数据库

    # 详见 https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
    # 下载物种注释，gb 1.8G; wgs 3.3G; 解压为9.2/21G; taxdump 50M
    kraken2-build --download-taxonomy --threads 24 --db ${db}/kraken2

    # 数据库包括 archaea bacteria plasmid viral human fungi plant protozoa nr nt env_nr env_nt UniVec
    # 下载单个库
    i=bacteria
    kraken2-build --download-library $i --threads 24 --db ${db}/kraken2
    # 批量下载数据库，除默认5种外新加植物、真菌、原生生物和质粒，下载需几小时-数天
    for i in archaea bacteria UniVec_Core viral human fungi plant protozoa plasmid; do
    kraken2-build --download-library $i --threads 24 --db ${db}/kraken2
    done
    
    # 建索引，4h, 40h
    time kraken2-build --build --threads 24 --db ${db}/kraken2
    
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

### metawrap分箱

运行metaWRAP annotate_bins时报错？

This copy of tbl2asn is more than a year old.  Please download the current version