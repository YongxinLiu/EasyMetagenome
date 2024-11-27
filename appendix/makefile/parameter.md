SHELL:=/bin/bash

	# 宏基因组第一版流程配置文件，请根据项目具体情况修改
	# Config file of Metagenome pipeline version 1, please modify according to experiment design
	# v1.0 2019/4/23 编写标准流程
	# v1.1 2020/12/28 测试完整流程
	# v1.2 2021/1/19 跨系统测试流程
	# v1.3 2021/4/8 增加分箱功能 metawrap1.3

# 1. 有参分析流程参数 Parameters of reference-based pipeline

	# 工作目录 Working directory
	# 修改wd为当前工作目录pwd
	wd=`pwd`
	# 设置检查点样本名，第一个(tail -n+2 result/metadata.txt|head -n1|cut -f1)或最后一个(tail -n1 result/metadata.txt|cut -f1)
	# i=`tail -n1 result/metadata.txt|cut -f1`
	# 设置j最大运行任务/p线程数/p1非并行任务线程数
	# 超过总线程效率会降低, htop 或 cat /proc/cpuinfo查看线程数
	j=3
	p=32
	p1=64
	# 设置数据库目录
	db=/db
	# 调置临时文件，/tmp空间不中会导致程序中断，备选 /mnt/m3/yongxin/tmp
	tmp=/mnt/m3/yongxin/tmp
	# make init # 建立分析所需子目录
	# 准备实验设计(result/design.txt)和测序数据(seq/*.fq.gz)和数据库(修改如下参数)


## 1.1. qc 质控和去宿主
	
	# 质控软件trimmomatic安装目录
	trimmomatic_path=/conda/envs/meta/share/trimmomatic/
	# 宿主基因组bowtie2索引，如人human_genome/Homo_sapiens, 拟南芥ath/水稻rice/水麦wheat/苜蓿medicago/bt2
	host_bt2=/db/kneaddata/medicago/bt2
	# 接头文件：通常cleandata已经去除接头，因无须指定接头文件。检查multiqc中Adapter Content 如有接头，则查找接头文件trimmomatic/adapters中文件，手动指定参数，如 ILLUMINACLIP:/conda/envs/meta/share/trimmomatic/adapters/TruSeq2-PE.fa:2:40:15

	# Kraken2过滤方法，去除植物和人(nohost)或选择微生物(microbe)，包括细菌、古菌、病毒、真菌和原生动物
	kraken2_qc_filter=nohost
	# 抽样方法，可选保持不变(raw)或截取头部(head)
	subsample=head
	# 输入文件结尾fq, fastq, fq.gz, fqstq.gz; 默认fq.gz
	fq=fq.gz
	# kraken2注释中筛选reads数据
	read_max=20000000
	# 默认截取6G(20M reads, 80M lines)，快速分析可选4(3G)或2(1.5G)千万行
	n=80000000

## 1.2. 物种和功能组成定量 humman2


## 1.3. 整理物种组成表和基本绘图 Summary metaphlan2 and plot
	
	# 绘制热图高丰度菌数据
	tax_top=30

## 1.4 kraken2特种注释
	
	kraken2_db=/db/kraken2/210120pfp
	kraken2_header=`tail -n+2 result/metadata.txt | cut -f 1 | head -n1`
	# 设置序列长度，通常为150，也可能为100，250或300
	bracken2_length=150
	# 设置估算的分类级别D,P,C,O,F,G,S
	# bracken2_tax=P
	# t是阈值，默认为0，越大越可靠，但可用数据越少
	bracken2_t=0


# 2. 无参分析流程 De novo assemble pipeline

## 2.1.

## 2.2. 组装 Assemble

	# megahit参数
	## 最小contig长度，默认200，推荐500
	min-contig-len=500
	# 最小kmin，默认21，越小计算量越大，土壤推荐27，25-31范围合理
	# 超过>300G数据至少29，数据量大可进一步增加39/49/59，诺禾中用K57计算
	kmin=29
	# 最长kmer，默认141适合PE150数据，PE100可设为91
	kmax=141
	# kmer步长，默认12，最小10低覆盖度组装更好，最大28，越大精度越差但计算更少
	kstep=28

	# 组装软件选择 Choose assemble method: megahit / metaspades
	# 组装模式 single/group/all
	# 单样品组装，适合大数据megahit_single；按组最合理megahit_group；全部适合小样本megahit_all
	
	# 评估长阈值，默认500，binning可达2k
	quast_len=1000

	# 单样品组装，可选 megahit metaspades、
	assemble_method=megahit


## 2.3. 基因组注释 Genome annotation

	# 基因注释软件选择 Choose assemble method: prokka / prodigal(基因存在多转录本，暂不考虑)
	annotation_method=prokka

### 2.3.3 构建非冗余基集 Non-redundancy gene set(大数据可选)

	# -c sequence identity threshold, default 0.9; -M  max available memory (Mbyte), default 400; -l  length of throw_away_sequences, default 10
	# -n 10 -d 0 -g 1 -b ${band_width} 
	cdhit_coverage=0.9
	cdhit_similarity=0.95
	cdhit_mem=900000
	# If you need highly identical matches, you may narrow the alignment band width, using options such as -b 10 or -b 5, instead of the default -b 20. This will make it faster. But with a smaller band width, similar sequences with larger gaps will not be clustered.
	# 默认为20；加速10，5
	band_width=20
	# time cd-hit-est -i temp/prodigal_all/gene.fa -o temp/NRgene/gene5.fa -aS 0.9 -c 0.95 -G 0 -g 0 -M 900000 -T 48 -b 5

## 2.3. 物种注释 kraken2
	
	# 物种注释分三个级别：reads、contig、gene
	py3=/conda/envs/humann3/bin

## 2.4 功能数据库注释

### 2.4.1 eggNOG
	
	# 拆分行2000000，1M行
	# split_line=10000
	# eggnog数据库位置 483m
	eggnog_db=/db/eggnog
	# 复制数据库至内存22G，加速检索(前提内存足够大) cp /home/meta/db/eggnog2 /dev/shm
	# eggnog_db=/dev/shm

### 2.4.2 KEGG

	# 76,4G; 2020,13G; /mnt/m1/liufang/db/KEGG_2020_DIAMOND/diamond_2.0.2.140/KEGG_diamond_2.0.2.140.dmnd
	# /db/kegg/genes/fasta/diamond2.0.5.dmnd
	# cp /mnt/m1/liufang/db/KEGG_2021_DIAMOND/kegg_2021.dmnd /home/meta/db/kegg/
	kegg_dmnd=/db/kegg/kegg_2021.dmnd
	gene_ko=/db/kegg/genes/ko/gene_ko.list
	# --sensitive，慢7倍，可填空
	kegg_mode=--sensitive
	kegg_evalue=1e-3
	
	# 标准化单位，默认为一百万tpm/rpm，可选1、100或1000
	unit=1000000
	
### 2.4.3 CAZy 碳水化合物数据库

	# http://cys.bios.niu.edu/dbCAN2/download/
	dbcan2_dmnd=/db/dbcan2/CAZyDB.07312020.dmnd
	dbcan2_anno=/db/dbcan2/fam_description.txt

### 2.4.4 CARD/ResFams 碳水化合物数据库

	# http://www.dantaslab.org/resfams
	resfams_dmnd=/mnt/zhou/yongxin/db/ResFams/Resfams-proteins
	resfams_anno=/mnt/bai/yongxin/data/db/ResFams/Resfams-proteins_class.tsv

# 3. 分箱

	# 筛选contigs长度
	l=1000
	# 定义完整度和污染率的阈值(50, 5; Finn NBT 2020;50, 10, Bowers NBT 2017)
	c=50
	x=10
	# 组装模式，megahit 或 metaspades
	mode=megahit
	# 分箱工具选择，默认为--metabat2 --maxbin2，可加concoct或为空
	bin_tools=--concoct
	# 内存(G), 默认25，设置总内存/并行数，即2000/3=667。
	mem=667

## 3.1 混合分箱 bin_all


## 3.2 分组分箱 bin_batch

	# 批次列表
	bin_list=`tail -n+2 result/metadata.txt|cut -f2|sort|uniq`
	g=Lyk9

## 3.3 单样本分箱 bin_single



# 4. 统计绘图

	# 绘图通用参数
	# 实验设计文件位置，全局，其它图默认调此变量，也可单独修改；并选择表中的组列和具体分组
	# 设置子版本目录
	sub=""
	doc=${wd}/result/${sub}
	design=${wd}/result/metadata.txt 
	g1=GroupID
	# tail -n+2 ${doc}/design.txt|cut -f 5|sort|uniq|awk '{print "\""$1"\""}'|tr "\n" ","
	# 绘图使用的实验组，顺序即图中显示顺序；为空时使用所有组和默认顺序
	#g1_list='"Col","ThasKO2","ThahKO","ThadKO","ACT2KO"'
	# 从实验设计比较组中提取组名，自动获得目录组 (推荐)
	g1_list=`cat result/${sub}/compare.txt|tr '\t' '\n'|sort|uniq|awk '{print "\""$$1"\""}'|tr "\n" ","|sed 's/,$$//'`
	# 从实验设计提取组(可选)
	# g1_list=`tail -n+2 ${doc}/design.txt|cut -f 5|sort|uniq|awk '{print "\""$$1"\""}'|tr "\n" ","|sed 's/,$$//'`

	# 组间比较列表
	compare=${doc}/compare.txt
	# 组间共有、特有比较列表
	venn=${doc}/venn.txt
	# 图片长宽，按nature全面版、半版页面设置
	# 图片长宽和字体大小，7组以下用默认，7组以上改为8x5或更大； figure size, recommend 4x2.5, 5x3(default), 8x5, 16x10, text_size 6, 7(default), 8
	width=8
	height=5
	text_size=7

	# 图中显示legend, 如taxonomy的数量，5，8(default)，10
	legend_number=10
	# 差异统计按丰度过滤 abundance filter，如丰度按万分之一过滤，减少计算量，提高OTU的FDR值，根据组数量多少可选十万5或万分之5
	abundance_thre=0.01
	# 差异比较方法wilcox/t.test/edgeR，默认wilcox
	compare_method="wilcox"
	# 显著性P值过滤 threshold of P-value，可选0.05, 0.01, 0.001。采用FDR校正，此参数意义不大，即使0.001也没有FDR < 0.2过滤严格
	pvalue=0.05
	# 统计检验方式FDR，常用0.05, 0.1, 0.2; FDR < 0.1使用9.5万次，且为菌群近期的Nature和Sciences; 0.2使用7.7万次
	FDR=0.2
	# 差异变化倍数常用1.5, 2, 4倍，对应logFC为0.585, 1, 2；菌丰度变化倍数不明显，还可用1.3和1.7倍对应0.379和0.766
	FC=1.2

	# 统计绘图和网页报告版本控制
	species="species"
	keyword="keyword"
	version=${species}_${keyword}_${sub}_v1

## 2.6 DA_compare 丰度筛选和组间差异比较
	# DA_compare2 Mearable OTU筛选和组间差异比较两步法
	Dc_input=result/kraken2_gene/proto.Genus.raw.txt
	Dc_compare=${compare}
	Dc_method=${compare_method}
	Dc_pvalue=${pvalue}
	Dc_FDR=${FDR}
	Dc_FC=${FC}
	Dc_thre=${abundance_thre}
	Dc_design=${design}
	Dc_group_name=${g1}
	Dc_group_list=${g1_list}
	Dc_output=${wd}/result/compare/



# 4. 细菌基因组

## 4.1 质控

## 4.2 组装

## 4.3 评估

## 去冗余drep(可选)

	# 0.99995 去重复, 0.99 株水平, 0.95 种水平；与MAG同一标准，nc0.3，完整50，污染10。
	sa=0.95
	# 纯菌一般都99%的完整度，nc可用90%覆盖度去冗余，但MAG只有50%完整度，仅用nc0.3
	nc=0.3

## 4.4 GTDB物种注释

	# gtdb基因组输入目录，默认为temp/antismash/，可选聚类目录，如 temp/drep/${sa}/dereplicated_genomes/
	gtdb_in=temp/antismash/
	# 输出目录，默认为空，非冗余为 ${sa}
	drep=
	


# 5. 三代宏基因组



include /home/meta/soft/Metagenome/denovo1/pipeline.md
