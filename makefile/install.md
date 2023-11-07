# 宏基因组软件和数据库安装

	# 测试环境为 Ubuntu 16/18.04 LTS / CentOS 7.3
	# Ubuntu登陆时显示系统版本，可用`cat /etc/issue`命令查看；CentOS用`cat /etc/redhat-release`查看

	# 设置数据库存放位置，如`~/db`目录，务必存于自己有权限的位置
	# 本此设置为`/db`根目录下，方便多人使用，需要管理员权限
	db=/db
	mkdir -p ${db}


## 软件管理器miniconda2

	cd 
	mkdir -p software && cd software
	wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
	bash Miniconda2-latest-Linux-x86_64.sh
	# 正常默认安装即可(安装目录为`~/miniconda2`)，为方便多人使用目录更短，我安装在`/conda`目录
	# 安装说明详见：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)
	# 添加通道
	conda config --add channels conda-forge
	conda config --add channels bioconda
	# 添加清华镜像加速下载
	site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
	conda config --add channels ${site}/pkgs/free/ 
	conda config --add channels ${site}/pkgs/main/
	conda config --add channels ${site}/cloud/conda-forge/
	conda config --add channels ${site}/pkgs/r/
	conda config --add channels ${site}/cloud/bioconda/
	conda config --add channels ${site}/cloud/msys2/
	conda config --add channels ${site}/cloud/menpo/
	conda config --add channels ${site}/cloud/pytorch/


## 质控软件

	# 新建环境
	conda create -n meta
	conda activate meta
	# 序列处理
	conda install seqkit
	# 表格处理
	conda install csvtk
	# 直接下载最新版 https://github.com/shenwei356/csvtk/releases
	wget -c https://github.com/shenwei356/csvtk/releases/download/v0.22.0/csvtk_linux_amd64.tar.gz
	tar xvzf csvtk_linux_amd64.tar.gz
	cp csvtk ~/miniconda3/bin/
	# 质量评估软件fastqc
	conda install fastqc # 0.11.9
	# 多样品评估报告汇总multiqc
	conda install multiqc # v1.8
	# 质量控制流程kneaddata
	conda install kneaddata # v0.7.4
	# 相关宿主索引布署，如人类、核糖体、拟南芥和水稻，请参考附录 KneadData


## 有参分析流程MetaPhlAn2、HUMAnN2、Kraken2

	# 安装MetaPhlAn2、HUMAnN2和所有依赖关系
	conda create -n humann2 humann2 -c bioconda # 2.8.1 
	conda activate humann2
	# 测试流程是否可用
	humann2_test # Ran 185 tests in 146.291s OK
	# 数据库布置见附录：HUMAnN2

	# metaphlan2输助脚本
	conda install export2graphlan
	conda install graphlan
	# ImportError: libgfortran.so.1错误解决办法
	conda install libgfortran==1


	# 物种注释
	# 基于LCA算法的物种注释kraken2  https://ccb.jhu.edu/software/kraken/
	conda install kraken2 # 2.0.7_beta
	# 下载数据库
	kraken2-build --standard --threads 24 --db $db/kraken2
	# 此步下载数据>50GB，下载时间由网速决定，索引时间4小时33分，多线程最快35min完成
	# 标准模式下只下载5种数据库：古菌archaea、细菌bacteria、人类human、载体UniVec_Core、病毒viral
	# 添加数据库并重建索引，如添加真菌库
	kraken2-build --download-library fungi --db $db/kraken2
	# 数据库确定后，建索引
	kraken2-build --build --db $db/kraken2

	# 新位置定义不同类型数据库
	cd $db
	cp -r kraken2 kraken2mod
	# 下载物种注释
	kraken2-build --download-taxonomy --threads 24 --db ${db}/kraken2mod
	# 下载非默认数据库archaea bacteria plasmid viral human己下载，可选fungi plant protozoa，以及nr nt env_nr env_nt UniVec
	kraken2-build --download-library fungi --threads 24 --db ${db}/kraken2mod
	kraken2-build --download-library plant --threads 24 --db ${db}/kraken2mod
	kraken2-build --download-library protozoa --threads 24 --db ${db}/kraken2mod
	# 批量下载，非标准数据库
	for i in nr nt env_nr env_nt UniVec; do
		kraken2-build --download-library $i --threads 24 --db ${db}/kraken2mod
	done
	# 确定的库建索引
	kraken2-build --build --threads 24 --db ${db}/kraken2mod

	# 直载Github脚本，需要git程序，没有程序可以wget或网站上手动下载压缩包，再解压
	# LangilleLab输助脚本 https://github.com/LangilleLab/microbiome_helper
	git clone git@github.com:LangilleLab/microbiome_helper.git
	# 我写的输助脚本 https://github.com/YongxinLiu/Metagenome
	git git@github.com:YongxinLiu/Metagenome.git
	# 添加环境变量，否则需要使用脚本所在的完整路径
	# 添加目录至环境变量，注释修改位置
	echo "export PATH=`pwd`/microbiome_helper:`pwd`/Metagenome/denovo1/script/:\$PATH" >> ~/.bashrc

## 基因组拼接、注释和定量

	# megahit 快速组装
	conda install megahit # v1.1.3
	# QUEST 组装评估
	conda install quast # 5.0.1 

	# prokka 细菌基因组注释
	conda install prokka # 1.13.3
	# cd-hit 非冗余基因集
	conda install cd-hit # 4.6.8
	# emboss transeq工具翻译核酸为蛋白
	conda install emboss # 6.6.0

	# 定量工具salmon
	conda install salmon # 0.11.3


## 分箱工具

### Metawrap

	# 物种注释和分箱流程 https://github.com/bxlab/metaWRAP
	conda create -n metawrap python=2.7
	source activate metawrap
	conda config --add channels ursky
	conda install -c ursky metawrap-mg
	# 数据库见附录 metawrap


## 其它工具

	# 比对结果整理samtools
	conda install samtools


### 其它输助脚本

	# 按注释类型同类合并脚本，如KO级别合并
	cd ~/bin
	wget http://bailab.genetics.ac.cn/share/mat_gene2ko.R
	chmod +x mat_gene2ko.R

	# Bin可视化VizBin
	sudo apt-get install libatlas3-base libopenblas-base default-jre
	curl -L https://github.com/claczny/VizBin/blob/master/VizBin-dist.jar?raw=true > VizBin-dist.jar



# 附录 Appendix

  db=/db
  cd $db

## 核糖体数据库

	mkdir -p $db/rDNA

### Greengene

  # greengene13.5 下载地址：ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz

	mkdir -p $db/rDNA/gg && cd $db/rDNA/gg
	wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz
	tar xvzf gg_13_5_otus.tar.gz

### usearch

usearch提供gg13.5, rdpv16和silva123序列和物种注释文件

	mkdir -p $db/rDNA/usearch && cd $db/rDNA/usearch
	# RDP training set v16 (13k seqs)，数量少，但预测更准确
	wget -c http://www.drive5.com/sintax/rdp_16s_v16_sp.fa.gz
	# Greengenes v13.5 (1.2M seqs)
	wget -c http://www.drive5.com/sintax/gg_16s_13.5.fa.gz
	# SILVA v123 (1.6M seqs)
	wget -c http://www.drive5.com/sintax/silva_16s_v123.fa.gz
	# silva精选 SILVA v123 LTP named isolate subset (12k seq)
	wget -c http://www.drive5.com/sintax/ltp_16s_v123.fa.gz
	# 真菌ITS序列 (53k sequences in v7.1)
	wget -c https://unite.ut.ee/sh_files/utax_reference_dataset_22.08.2016.zip
	# RDP ITS training set v2 (18k sequences)
	wget -c http://www.drive5.com/sintax/rdp_its_v2.fa.tz
	# 解压全部gz文件
	gunzip *.gz


## KneadData宿主基因组HostGenome

### 人类和核糖体

	# 查看可用数据库，如宏基因组去宿主，宏转录组去核糖体
	kneaddata_database
	# 如下载人类基因组bowtie2索引，3.44G，推荐晚上过夜下载
	mkdir -p $db/host/human
	kneaddata_database --download human_genome bowtie2 $db/host/human
	# 如宏转录组用SILVA核糖体数据库，3.61G
	mkdir -p $db/host/silva
	kneaddata_database --download ribosomal_RNA bowtie2 $db/host/silva/

### 拟南芥基因组索引

	# http://plants.ensembl.org/info/website/ftp/index.html
	mkdir -p $db/host/ath
	cd $db/host/ath
	wget -c ftp://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
	gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
	mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa tair10.fa
	bowtie2-build --threads 9 tair10.fa bt2 # 28s

### 多水稻混合基因组

	# 下载建索引，以水稻籼粳混合为例，平时只有单个基因组(涉及多物种需多个基因组)
	# 设置基因组目录，改为自己有权限的目录，我设置的位置方便大家使用 Set genome download directory
	mkdir -p $db/host/rice && cd $db/host/rice
	# Ensembl plant: http://plants.ensembl.org/info/website/ftp/index.html
	# Download Oryza indica 
	wget -c ftp://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/oryza_indica/dna/Oryza_indica.ASM465v1.dna.toplevel.fa.gz
	# Download Oryza sativa (japonica)
	wget -c ftp://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
	gunzip *.gz
	rename 's/Oryza_//;s/dna.toplevel.//;' *.fa
	# combine indica and japonica
	cat <(sed 's/>/>IND/' indica.ASM465v1.fa) <(sed 's/>/>TEJ/' sativa.IRGSP-1.0.fa) > ebi41.fa
	# bowtie2 index 2.2.6
	bowtie2-build --version
	time bowtie2-build --threads 9 ebi41.fa bt2
	# 1p 37m55s, 9p 7m55s

### 小麦基因组

	#http://plants.ensembl.org/info/website/ftp/index.html
	mkdir -p $db/host/wheat
	cd $db/host/wheat
	wget -c ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
	gunzip Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
	mv Triticum_aestivum.IWGSC.dna.toplevel.fa tae.fa
	bowtie2-build --threads 9 tae.fa bt2 # 1h41m

### 苜蓿基因组
  
  # R108测序只有1.0版，目前最好的是A17 4.0 https://www.ncbi.nlm.nih.gov/assembly/GCF_000219495.3 点RefSeq download
	mkdir -p $db/host/med
	cd $db/host/med
	# 下载时间：2019/9/19，文件修改时间：2016/10/18
	wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/219/495/GCA_000219495.2_MedtrA17_4.0/GCA_000219495.2_MedtrA17_4.0_genomic.fna.gz 
	gunzip GCA_000219495.2_MedtrA17_4.0_genomic.fna.gz 
	mv GCA_000219495.2_MedtrA17_4.0_genomic.fna MedtrA17_4.fa
	bowtie2-build --threads 9 MedtrA17_4.fa bt2 # 5m


## HUMAnN2有参流程

	humann2_databases # 显示可用数据库
	wd=/db/humann2
	mkdir -p $wd # 建立下载目录
	cd $wd
	# 蛋白注释库 0.58G
	humann2_databases --download utility_mapping full $wd 
	# 微生物泛基因组数据库 5.37G
	humann2_databases --download chocophlan full $wd 
	# 非冗余功能基因diamond索引 10.3G
	humann2_databases --download uniref uniref90_diamond $wd
	
  # 	# 手动下载最新版 http://huttenhower.sph.harvard.edu/humann2_data
  # 	wd=/db/humann2v201901
  # 	mkdir -p $wd
  # 	cd $wd
  # 	# 解压前备份至 humann2v201901/ 为 humann2v201901_bak/
  # 	# 201901版泛基因组数据库 5 to 15G
  #   wget -c http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.v201901.tar.gz
  #   mkdir -p chocophlan
  #   tar xvzf full_chocophlan.v201901.tar.gz -C chocophlan
  # 	# 201901非冗余功能基因diamond索引 10 to 20G
  #   wget -c http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_v201901.tar.gz
  #   mkdir -p uniref90
  #   # 解压20G为34G，用时20m
  #   time tar xvzf uniref90_annotated_v201901.tar.gz -C uniref90
  #   # 可选小内存使用uniref50，只从2.7 to 6.9 G
  # 	wget -c http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref50_annotated_v201901.tar.gz
  #   mkdir -p uniref50
  #   tar xvzf uniref50_annotated_v201901.tar.gz -C uniref50
	# 1万条数据测试新数据库
  # head -n 4000000 temp/12concat/lyr4B3R2.fq > temp/temp.fq
  # time humann2 --input temp/temp.fq --output temp/12humann2/
  # 201901数据库不兼容，还在开发中... Error: Incompatible database versio


	# 设置数据库位置和线程数
	# 复制 cp -r ../humann2/utility_mapping ./
	humann2_config --update database_folders utility_mapping $wd/utility_mapping
	humann2_config --update database_folders protein $wd/uniref
	humann2_config --update database_folders nucleotide $wd/chocophlan
	humann2_config --update run_modes threads 24
	# 显示参数
	humann2_config --print
	# metaphlan2数据库默认位于程序所在目录的db_v20和databases下各一份，备用文件见/db/humann2/metaphlan2/
	# # The database file for MetaPhlAn does not exist at /conda/bin/db_v20/mpa_v20_m200.pkl . Please provide the location with --metaphlan-options .
	ln /db/humann2/metaphlan2/ /conda/bin/db_v20 -s 
	# ERROR: Unable to create folder for database install: /conda/bin/databases
	# --metaphlan-options '--mpa_pkl ''
	ln /db/humann2/metaphlan2/ /conda/bin/databases -s 

## 基因功能注释

	mkdir -p $db/protein

### COG/eggNOG

	# http://eggnogdb.embl.de
	# 安装eggnog比对工具
	conda install eggnog-mapper
	# 下载常用数据库，注意设置下载位置
	mkdir -p $db/eggnog && cd $db/eggnog
	download_eggnog_data.py --data_dir ./ -y -f euk bact arch viruses
	# 如果内存够大，复制eggNOG至内存加速比对
	cp $db/eggnog/eggnog.db /dev/shm/
	# 手工整理COG分类注释
	wget -c wget http://bailab.genetics.ac.cn/share/COG.anno
	# 手工整理KO注释
	wget -c wget http://bailab.genetics.ac.cn/share/KO.anno


### 碳水化合物数据库dbCAN2 

	#http://cys.bios.niu.edu/dbCAN2/
	mkdir -p $db/protein/dbCAN2 && cd $db/protein/dbCAN2
	wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fa
	wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fam-activities.txt
	diamond makedb --in CAZyDB.07312018.fa --db CAZyDB.07312018
	# 提取fam对应注释
	grep -v '#' CAZyDB.07312018.fam-activities.txt|sed 's/  //'| \
	sed '1 i ID\tDescription' > fam_description.txt


### 抗生素抗性基因Resfam

	# http://dantaslab.wustl.edu/resfams
	mkdir -p $db/protein/resfam && cd $db/protein/resfam
	# 官网的数据格式非常混乱, 推荐下载我手动整理的索引和注释
	wget http://bailab.genetics.ac.cn/share/Resfams-proteins.dmnd
	wget http://bailab.genetics.ac.cn/share/Resfams-proteins_class.tsv


## metawrap

	# 相关数据库，大小近300GB
	# 这里我们安装数据库到`db=/db`目录，其它人请保证你有权限，如改为`~/db`目录，但要保证至少有500GB的空间。
	# 请根据你的情况修改为自己有权限且空间足够的位置。多人使用，建议管理员统一安装节省空间
	# 务必启动虚拟环境再进行配置
	source activate metawrap


### CheckM

	# 用于Bin完整和污染估计和物种注释
	cd $db
	mkdir checkm && cd checkm
	# 下载文件276MB，解压后1.4GB
	wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
	tar -xvf *.tar.gz
	rm *.gz
	# 设置数据库位置
	checkm data setRoot
	# 按提示输出你数据下载的路径，按两次回车完成后退出


### KRAKEN物种注释数据库

	# 下载建索引需要 > 300GB以上空间，完成后占用192GB空间
	cd $db
	mkdir -p $db/kraken
	kraken-build --standard --threads 24 --db $db/kraken
	kraken-build --db $db/kraken --clean


### NCBI_nt核酸序列用于bin物种注释

	# 41GB，我下载大约12h；解压后99GB
	cd $db
	mkdir NCBI_nt && cd NCBI_nt
	wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
	for a in nt.*.tar.gz; do tar xzf $a; done


### NCBI物种信息

	# 压缩文件45M，解压后351M
	cd $db
	mkdir NCBI_tax
	cd NCBI_tax
	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
	tar -xvf taxdump.tar.gz


### 人类基因组去宿主

	mkdir -p $db/host/human && cd $db/host/human
	wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
	gunzip *fa.gz
	cat *fa > hg38.fa
	rm chr*.fa
	bmtool -d hg38.fa -o hg38.bitmask
	srprism mkindex -i hg38.fa -o hg38.srprism -M 100000


### 数据库位置设置

	which config-metawrap
	# 查使用vi/vim/gedit等文本编辑器来修改数据库的位置吧
	# 退出环境
	source deactivate



## 16S分析相关工具

### qiime

	conda install qiime # 1.9.1

### qiime2
  
	cd ~/soft
	wget https://data.qiime2.org/distro/core/qiime2-2018.11-py35-linux-conda.yml
	conda env create -n qiime2 --file qiime2-2018.11-py35-linux-conda.yml

	source activate qiime2
	source deactivate

### lefse

	conda install lefse
	# ImportError: libgfortran.so.3
	conda install libgfortran==3
	# ImportError: libicuuc.so.54
	conda install ICU==54 # UnsatisfiableError
	scp -r yongxin@192.168.0.110:/mnt/bai/yongxin/bin/miniconda3/pkgs/icu-54.1-0/ /conda/pkgs/ # 程序找不到
	# 设置lib的三种方法
	sys.path
	sys.path.append("/conda/pkgs/icu-54.1-0/")
	sys.path
	# 添加 ~/.bashrc PYTHONPATH=$PYTHONPATH:/conda/pkgs/icu-54.1-0

  
### picrust预测KO
  
	conda install picrust
	# 数据库下载
	mkdir -p $db/protein/picrust/ && cd $db/protein/picrust/
	download_picrust_files.py



## 其它参考代码段

### 测序数据批量提取相同行

  # 文件提取100万行测试
	mkdir -p 02seq/bak
	mv 02seq/*.gz 02seq/bak
	for i in `tail -n+2 01doc/design.txt|cut -f 1`; do
		zcat 02seq/bak/${i}_1.fq.gz | head -n 1000000 | gzip > 02seq/${i}_1.fq.gz
		zcat 02seq/bak/${i}_2.fq.gz | head -n 1000000 | gzip > 02seq/${i}_2.fq.gz
	done

### 三代工具

#### 二三代混装OPERA-MS

	# https://github.com/CSB5/OPERA-MS
	~/bin