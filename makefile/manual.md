	# 宏基因组分析流程第一版 —— 操作脚本
	# Metagenome pipeline version 1 —— Manual script

	# Biocloud中重置环境变量(可选)
	PATH=/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
	# 调置Perl模块
	PERL5LIB=/conda/envs/kraken2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/conda/envs/kraken2/lib/site_perl/5.26.2:/conda/envs/kraken2/lib/5.26.2/x86_64-linux-thread-multi:~/miniconda2/envs/kraken2/lib/5.26.2

	# 加载目标环境变量
	source /conda/bin/activate
	# 启动宏基因组通用分析环境
	conda activate meta

# 1. 有参分析流程 Reference-based pipeline

	# 0. 准备工作 Preparation

	# 设置工作目录 Set work directory
	wd=medicago/lyr4_210115b1/
	mkdir -p ~/${wd}
	cd ~/$wd
	
	# 准备流程 Prepare makefile
	ln -s /home/meta/soft/Metagenome/denovo1/parameter.md makefile
	ln -s /home/meta/soft/Metagenome/denovo1/manual.md manual.sh
	# cp /home/meta/soft/Metagenome/denovo1/manual.md manual.sh
	
	# 建立初始工作目录 Create initial working directory
	# 代码预览、执行、保存
	make -n -B init
	make init
	echo "#" `date` > result/pipeline.sh
	make -n -B init >> result/pipeline.sh

	# 准备原始数据 sequencing raw data (多样本合并和统计见附录1)
	# 链接数据至工作目录
	ln -s /mnt/m2/data/meta/$wd/*.gz seq/

	# 准备实验设计上传到result目录，至少有两列样本名和组名 Experiment design
	# 方法1. 数据来源处复制实验设计
	cp /mnt/m2/data/meta/$wd/metadata.txt result/metadata.txt
	# 方法2. 复制实验设计模板并手动填写
	cp /home/meta/soft/Metagenome/denovo1/result/design.txt result/metadata.txt
	# 方法3. 从样本名中提取，并手动补充
	ls seq/*_1.fq.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|awk '{print $1"\t"$1}'|sed '1 i SampleID\tGroupID' > result/metadata.txt

	# 检查实验设计ID唯一性，有结果则为重复
	cut -f1 result/metadata.txt|sort|uniq -d

## 1.1. 质控并移除宿主 Quality control & Remove host

	### 1.1.1 (可选)质量评估原始数据
	# 依赖软件版本：质量评估、评估报告汇总
	fastqc -v # v0.11.9
	multiqc --version # 1.8
	# 预览代码
	make -n -B qa
	# 计算：240G,72p,30m; 64p,37m,16G
	memusg -t make qa
	# 查看result/multiqc_report.html，重点检查数据量分布，重复率，接头含量
	echo "#" `date` >> result/pipeline.sh
	make -n -B qa >> result/pipeline.sh

	### 1.1.2 KneadData移除低质量和宿主
	echo -e "\n#" `date` >> result/pipeline.sh
	# 依赖软件版本：并行、移除低质量和宿主
	echo -n "kneaddata --version # " >> result/pipeline.sh
	kneaddata --version >> result/pipeline.sh 2>&1
	# 预览代码
	make -n -B qc
	# 计算：150G,24p,20h; 240G,72p,12h; ; 126G,24x3p,2.5h,99G; 800G,96p,31h,66G
	make qc
	# 异常中断处理见附录2. KneadData中断
	# 结果见 result/qc/kneaddata.txt 高质量、非宿主比例
	# 保存运行代码
	make -n -B qc >> result/pipeline.sh

	### 1.1.3 提取上传的Clean数据
	# 预览代码
	make -n -B qa2
	# 	# 240G,72p,1h，按GSA标准整理上传数据，见submit目录
	make qa2
	make -n -B qa2 >> result/pipeline.sh

	### 1.1.4 kraken2质控
	conda activate kraken2 # bicloud
	kraken2 --version # 2.1.1
	# 采用kraken完整数据库注释, 370G,24p,10h; 330G,24x3p,1h;
	make -n -B kraken2_qc
	make kraken2_qc
	make -n -B qa2 >> result/pipeline.sh

## 1.2. 物种和功能组成定量 humman2

	## 1.2.1 humman2输入文件准备：双端文件cat连接
	make humann2_concat

	# 启动humann2环境
	conda activate humann2
	humann2 --version # v2.8.1

	# 检查数据库位置
	humann2_config --print
	
	## 1.2.2 humman2计算，包括metaphlan2
	make humann2
	# 4.7 Tb水稻数据，8X12线程运行2.5 Days
	#  X Tb Data, 3 x 24，运行 xx

	## 1.2.3 功能组成整理 humman2_sum
	make humann2_sum
	# 结果见 result/12humann2目录，有功能通路及物种组成表uniref.tsv、标准化表uniref_relab.tsv，以及拆分功能表unstratified和功能物种对应表stratified
	
	# humann2转为KEGG
	humann2_regroup_table -i temp/humann2/genefamilies.tsv \
	  -g uniref90_ko -o temp/humann2/ko.tsv

## 1.3. 整理物种组成表和基本绘图 Summary metaphlan2 and plot

	# 结果见result/13metaphlan2目录

	### 1.3.1 整理物种组成表 Summary metaphlan2
	make metaphaln2_sum
	# 结果taxonomy*文件，有多级物种表.tsv、株水平表.spf用于STAMP分析，以及聚类热图_heatmap_top.pdf观察分组情况

	### 1.3.2 GraPhlAn图
	make metaphaln2_graphlan
	# 结果taxonomy_graphlan*.pdf，包括主图、图例和注释文本3个文件

	### 1.3.3 物种组成LEfSe差异分析(可选)
	# 依赖实验设计，分组比较末确定时，可选跳过此步
	make metaphaln2_lefse


## 1.4. kraken2物种组成

	### 1.4.1 基于NCBI完整基因组数据库的k-mer物种注释
	# Taxonomy assign by k-mer and based on NCBI database
	kraken2 --version # Kraken version 2.0.9-beta
	make -n -B kraken2_read
	# 500G,3x24p,3h; 100G,3x24p,11m
	make kraken2_read
	# record pipeline
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B kraken2_read >> result/pipeline.sh

	### 1.4.2 合并为矩阵 merge into matrix
	make kraken2_read_sum



# 2. 无参分析流程 De novo assemble pipeline


## 2.1. Assemble 组装

	# 启动组装环境
	conda activate meta
	megahit -v # v1.2.9

	# 方法1. 大项目推荐
	### 2.2.1 基于qc质控序列单样本拼接(大项目)
	make assemble_single
	# 66个样168G压缩数据分别装为5h 计算过程日志见 temp/22megahit/megahit.log

	# 方法2. 小项目推荐
	### 2.2.2 基于qc质控后序列拼接
	# metahit拼接所有样本, 496G,72p,3d; 38G,48p,16h; 360G,72p,58h; 297G,72p,8h; 142G,48p,7d；1.5d,243G
	make -n -B megahit_all
	make megahit_all
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B megahit_all >> result/pipeline.sh

	# 问题. 数据量大无法拼接
	# 126 - Too many vertices in the unitig graph (8403694648 >= 4294967294), you may increase the kmer size to remove tons
	# 需要增加k-mer，如21增加为29

	### 2.2.3 megahit_all_quast评估
	# 3G,15m; 2G,7m
	make -n -B megahit_all_quast
	make megahit_all_quast

	### 2.2.5 (可选)kraken2对contig物种注释
	# 10M,9G,48m; 2G,
	conda activate kraken2
	make -n -B kraken2_contig
	make kraken2_contig
	conda deactivate
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B kraken2_contig >> result/pipeline.sh

## 2.3. Genome annotation 基因组注释

	### 2.3.1 对单样品组装的每个contig文件基因注释(大数据可选)
	make prodigal_single

	### 2.3.2 对合并组装的单个contig文件基因注释(小样本可选)
	# 3G,10h; 9G,20h；2G,
	prodigal -v # V2.6.3: February, 2016
	make -n -B prodigal_all
	make prodigal_all
	
	# 并行处理，加速基因预测
	make -n -B prodigal_all_split
	# 2.5h,473M; 4h,768M
	make prodigal_all_split
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B prodigal_all_split >> result/pipeline.sh

	### 2.3.3 构建非冗余基集 Non-redundancy gene set(大数据可选)

	# 90%覆盖度，95%相似度下，拼接结果再聚类基本不变少，如7736减少为7729
	# 5,514,236聚类为4,879,276，72线程155m，8756m；
	# 15,889,664聚类为13,924,854，48线程1651m(1d)，71427m；
	# 4M,48p,2h; 11M,72p,7h
	cd-hit-est -v # 4.8.1 (built on Oct 26 2019)
	make -n -B NRgene
	make NRgene
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B NRgene >> result/pipeline.sh

	### 2.3.4 基因定量 salmon genes
	conda activate metawrap1.3 # 中salmon仅为0.13.1
	salmon -v # 1.3.0
	# 1. 建索引:-t转录本,-i索引; 4M,24p,15m; 14M,24p,83m; 22M,24p,2h,57G
	# 2. 定量: 38Gx4M,24px3,30m
	# 3. 合并样本为基因丰度表:50x22M,1h,46G;X2(count/TPM)
	# 4. 统计维度，基因count/TPM第一样本求和：1m,29M; 1m,5G; 1m,6G; 
	make -n -B salmon_gene
	make salmon_gene
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B salmon_gene >> result/pipeline.sh
	# 错误解决，详见 ### cd-hit合并多批基因salmon索引时提示ID重复

	# 统计mapping rate：比对率与统计求和一致
	# 20M genes,19m,30G
	memusg -t Rscript /db/script/table_stat.R -i result/salmon_gene/gene.count -o result/salmon_gene/gene.count
	# 以第二个样本为例求和 Lyk9CpB8R1，TPM~1000000，count为8647721.08/14684642=0.5888956
	csvtk -t summary -f 2:sum result/salmon_gene/gene.TPM
	i=`head -n1 result/salmon_gene/gene.count|cut -f2`
	# seqkit stat temp/qc/${i}_1_kneaddata_paired_1.fastq
	csvtk -t summary -f 2:sum result/salmon_gene/gene.count
	grep $i result/qc/kneaddata.txt|cut -f4
	grep 'Mapping rate' temp/salmon_gene/$i.quant/logs/salmon_quant.log|cut -f 2 -d '='|sed 's/ //g'

	### 2.3.5 基因物种注释 kraken2 annotate gene

	# 9G,24p,3m; 13G,24p,5m,55.4%; 1.7M count,1m,5G; 1m,
	make -n -B kraken2_gene
	make kraken2_gene
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B kraken2_gene >> result/pipeline.sh
	# 从gene层面去宿主(动物和植物),13G,32m,2.7G,22.06/22.18=99.45%; 9G,24m,13.5/13.92=96.9%; 
	# 从gene层面选择(细菌和古菌),13G,17m,2.5G,11.99/22.18=54.1%; 
	# 基因筛选，只筛选原核生物
	make -n -B kraken2_gene_proto
	make kraken2_gene_proto
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B kraken2_gene_proto >> result/pipeline.sh

	# 统计过滤原核注释比例，Keep: 13988920/17143886=81.6%
	# csvtk -t summary -f 2:sum result/salmon_gene/gene.count
	# csvtk -t summary -f 2:sum result/salmon_gene/gene.count.proto

## 2.4 功能数据库注释

### 2.4.2 eggNOG
	# biocloud上: conda activate biobakery
	emapper.py --version # 2.0.1
	diamond --version # 2.0.6
	# 预览代码
	make -n -B eggnog
	# eggnog注释
	# 15M 72p, 6h
	# The host system is detected to have 1082 GB of RAM. It is recommended to use this parameter for better performance: -c1
	# 默认diamond比对，采用：--more-sensitive  -e 0.001 --top 3 --query-cover 0 --subject-cover 0
	make eggnog
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B eggnog >> result/pipeline.sh

	# 22M,25m,95G
	make -n -B eggnog_sum
	make eggnog_sum
	# 注释基因数量 3781759/4879276=77.5; 9036098/13,924,854=64.9; 3025768/3886015=77.8; 7189141/9729931=73.9
	# 注释KO基因数量 2335006/4879276=47.8; 5836548/13,924,854=41.9; 1955582/3886015=50.3; 4545144/9729931=46.7
	# 合并，15Mx76,75m; 4Mx25,3m; 7Mx25,11m
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B eggnog_sum >> result/pipeline.sh

### 2.4.2 KEGG

	# 比对kegg 2021 数据库
	# 22M,64p,14.5h,1.31T
	make -n -B kegg
	make kegg
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B kegg >> result/pipeline.sh
	# 注释基因数量3853393/4879276=79%; 1919238/2550307=75%; 1413301/1780779=79%
	
	# 汇总为result/24kegg/kotab.count/tpm两个表，分别为原始count和rpm
	# 22M,10m,49G; 
	make -n -B kegg_sum
	make kegg_sum

	# 采用刘芳的脚本、和分层数据库
	# 22M,35m,59G
	make -n -B kegg_sum_liufang

	# 只筛选原核(细菌、古菌)
	make -n -B kegg_sum_proto
	make kegg_sum_proto

### 2.4.3 CAZy

	# blast比对到CAZy
	make -n -B dbcan2
	make -n -B dbcan2_sum

### 2.4.4 CARD
	
	# 启动rgi环境
	conda activate rgi
	
	make -n -B card

# 3 分箱

	conda activate metawrap1.3
	metawrap -v # 1.3.2
	megahit -v # v1.1.3
	metaspades.py -v # v3.13.0

## 3.1 混合分箱

## 3.2 分批次分箱

## 3.3 单样本分箱
	
	# 杂合低、污染低、速度快，但MAG较少
	make -n -B binning_single

## 3.3 Drep去冗余


## 3.4 物种注释和进化



# 5. 细菌基因组
	
	make init

## 5.1 质控

	conda activate meta

	### 1.1.1 (可选)质量评估原始数据
	make -n -B qa
	make qa
	echo "#" `date` > result/pipeline.sh
	make -n -B qa >> result/pipeline.sh
	# 提取样本数据量和统计总数据量
	cut -f 1,5,7 result/qc/multiqc_data/multiqc_fastqc.txt | sed 's/seq | //;s/_.\t/\t/' | uniq | awk '{print $0"\t"$2*$3*2}' | sed '1 s/0$/Total bases/' | less -S > result/qc/multiqc_data.txt
	csvtk -t summary -f 4:sum result/qc/multiqc_data.txt

	### 1.1.2 Trimmomatic质量评估原始数据(报错中需要复制代码手动运行)
	make -n -B trimmomatic
	make trimmomatic
	echo "#" `date` >> result/pipeline.sh
	make -n -B trimmomatic >> result/pipeline.sh

## 5.2 组装	
	# 务必spades>3.14.0才支持--isolate模式
	spades.py -v # 3.14.0
	make -n -B spades
	make spades
	echo "#" `date` >> result/pipeline.sh
	make -n -B spades >> result/pipeline.sh

## 5.3 评估
	conda activate metawrap1.3
	make -n -B checkm
	make checkm
	echo "#" `date` >> result/pipeline.sh
	make -n -B checkm >> result/pipeline.sh

## 5.4 分箱
	conda activate metawrap1.3
	make -n -B bin_isolate
	make bin_isolate
	echo "#" `date` >> result/pipeline.sh
	make -n -B bin_isolate >> result/pipeline.sh

# 细菌基因组


## 去冗余drep

	source /conda/bin/activate
	conda activate drep
	# 1287,6h,1.11T;1287,8.5h,1.11T;按100,0.99995,0.99,0.95聚类为1287,1056,899,586个基因组
	# 2689,
	make -n -B drep
	make drep
	echo "#" `date` >> result/pipeline.sh
	make -n -B drep >> result/pipeline.sh


## 物种注释GTDB

	make -n -B gtdb_isolate
	# 3k,32p,12h,4.6T
	make gtdb_isolate
	echo "#" `date` >> result/pipeline.sh
	make -n -B gtdb_isolate >> result/pipeline.sh
	# 苜蓿784个分离培养细菌基因组，其中538为未知种，5个为未知属

	make -n -B gtdb_isolate_stat
	# 3k,32p,12h,4.6T
	make gtdb_isolate_stat
	echo "#" `date` >> result/pipeline.sh
	make -n -B gtdb_isolate_stat >> result/pipeline.sh

	
## 物种注释Kraken2(可选)

	# 整合单菌为单条序列，然后再kraken2分类
	conda activate meta
	make -n -B kraken2_genome
	make kraken2_genome
	echo "#" `date` >> result/pipeline.sh
	make -n -B kraken2_genome >> result/pipeline.sh
	# 结果result/itol/ncbi_tax.txt，所有菌均注释到种，与软件采用阈值为0相关，尽量分给名称标签

## 基因预测
	
	# 生成MAG的列表
	ls temp/antismash/|cut -f1 -d '.'|awk '{print $0"\tMAG"}'|sed '1 i ID\tSource'|less -S>result/metadata.txt
	tail -n+2 result/metadata.txt|wc -l  # 2689
	prodigal -v # V2.6.3: February, 2016
	make -n -B prodigal_isolate
	make prodigal_isolate
	echo "#" `date` >> result/pipeline.sh
	make -n -B prodigal_isolate >> result/pipeline.sh

## KEGG注释

	diamond --version # v2.0.8
	make -n -B kegg_isolate
	make kegg_isolate
	echo "#" `date` >> result/pipeline.sh
	make -n -B kegg_isolate >> result/pipeline.sh
	# 2689个菌usearch计算距离，1h; 4K个菌,4.5h
	make -n -B kegg_isolate_stat
	make kegg_isolate_stat
	echo "#" `date` >> result/pipeline.sh
	make -n -B kegg_isolate_stat >> result/pipeline.sh

## CAZy注释
	
	diamond --version # v2.0.8
	make -n -B dbcan2_isolate
	make dbcan2_isolate
	echo "#" `date` >> result/pipeline.sh
	make -n -B dbcan2_isolate >> result/pipeline.sh

	make -n -B dbcan2_isolate_stat
	make dbcan2_isolate_stat
	echo "#" `date` >> result/pipeline.sh
	make -n -B dbcan2_isolate_stat >> result/pipeline.sh

## CARD注释

	# activate env 启动环境
	conda activate rgi
	# 使用RGI程序基于card数据库注释
	make -n -B card_isolate
	make card_isolate
	# 目前运行有warining，但不一定影响结果




# 附录

## 依赖软件版本

	# 依赖软件版本：改名、压缩、链接、md5值计算、
	parallel --version # 20201122
	csvtk version # v0.22.0
	rename --version # 0.20
	pigz --version # 2.4
	ln --version # (GNU coreutils) 8.28, 包括Linux常用命令md5sum, cat
	sed --version # 4.8
	md5sum --version # 8.28

## 附录1. 宏基因组样本多文件合并和统计

	# 水稻测试数据6个样品，实验和对照各3个，数据量109-236M，PE100, 21.8-47.2G，共198.8GB
	# 拟南芥测试数据36个样本，4个实验组，共850GB数据
	# 样有多个文件，合并各样品文件()
	# merge_sample ~ 5h, 合并单个样品, 20GB 10min;
	cd /mnt/m2/data/meta/ath/3T/
	# 按实验设计按文件夹批量合并再改名，需要输入文件每个样本一个目录
	p=3
	for i in `tail -n+2 design.txt|cut -f3`; do
		zcat `find 01.filter/${i}/ -name *.gz | grep '_1.fq'` | pigz -p ${p} > seq/${i}_1.fq.gz &
		zcat `find 01.filter/${i}/ -name *.gz | grep '_2.fq'` | pigz -p ${p} > seq/${i}_2.fq.gz &
	done
	awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_1.fq.gz seq/"$1"_1.fq.gz ");system("mv seq/"$3"_2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 design.txt)
	# 统计样本md5值
	cd seq
	md5sum *_1.fq.gz > md5sum.txt
	md5sum *_2.fq.gz >> md5sum.txt
	cat md5sum.txt

