# 处理 taxonomy2.spf 文件，生成 taxonomy_merged.spf
cut -f 1-7,9- metaphlan4/taxonomy.spf | awk -F '\t' '
NR==1 {
    print;  # 打印表头
    header = $0;  # 保存表头
    next;  # 跳过表头
}
{
    taxonomy = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6; 
    species = $7;  # 获取物种
    for (i=8; i<=NF; i++) {
        # 统计出现次数，保持0，非0改为1
        if ($i > 0) {
            sum[species, i] += 1;  # 不为0则计为1
        } else {
            sum[species, i] += 0;  # 为0保持0
        }
    }
    species_seen[species] = taxonomy; 
}
END {
    for (sp in species_seen) {
        printf "%s\t%s", species_seen[sp], sp;
        for (i=8; i<=NF; i++) {
            printf "\t%d", sum[sp, i]; 
        }
        print "";
    }
}' > metaphlan4/taxonomy_merged.spf

# 读取已生成的 taxonomy_merged.spf 文件并进行修改，输出到 taxonomy_modified.spf
awk -F '\t' '
NR==1 {
    print;  # 打印表头
    header = $0;  # 保存表头
    next;  # 跳过表头
}
{
    # 打印前7列
    printf "%s\t%s", $1, $2;  # 打印前两列
    printf "\t%s\t%s\t%s\t%s\t%s", $3, $4, $5, $6, $7;  # 打印后五列

    # 从第8列开始进行处理
    for (i=8; i<=NF; i++) {
        if ($i > 0) {
            printf "\t1";  # 非0的值改为1
        } else {
            printf "\t0";  # 0保持为0
        }
    }
    print "";  # 换行
}' metaphlan4/taxonomy_merged.spf > metaphlan4/taxonomy_modified.spf


rm -rf metaphlan4/taxonomy_merged.spf
