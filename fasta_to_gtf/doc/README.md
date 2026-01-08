## 前言
根据定制化项目客户的需求，研发同事会设计靶向扩增目标基因片段的探针，需要针对该片段序列生成相应的 gtf 用于定量，现整理相关内容如下。

## 解决方案
step1：将目标基因片段和参考基因组比对得到bam文件，使用 bamToBed 获取 bed12 文件；

step2：使用 bedToGenePred + genePredToGtf 生成 gtf 文件。


## 脚本
/SGRNJ06/randd/USER/wangjingshen/script/fasta_to_gtf/script/fasta_to_gtf.py

#### 主要参数

--genome ： 基因组位置

--fasta ： 目标基因片段序列

--star_path ： star的位置

--prefix ： 前缀名

--force ： 是否重新生成 bam

--threads： star比对的线程数

#### 输出

{prefix}.gtf：  gtf 文件

#### 测试示例

/SGRNJ06/randd/USER/wangjingshen/script_dev/fasta_to_gtf/test/


## 特殊情况

#### 基因片段过长，超过 star 的比对长度上限（649bp）

解决方案：

1）首先把过长的 reads 调整到 fasta 文件的最后或者直接去除，避免该条 reads 之后的 reads 没有比对；
        
2）使用其他比对工具对该条 reads 进行比对，例如 blastn。根据比对得到的区域修改 bed 文件，chromStart 为比对起始点 -1，chromEnd 为比对终止点，blockSizes 为比对终止点 - 比对起始点。

#### 基因有多个 exon

解决方案：合并 bed 文件里多个 exon。具体而言：

1）chromStart： 最小的 chromStart；

2）chromEnd：最大的 chromEnd；

3）blockSizes： 按 chromStart 排序的 exon 长度；

4）blockStarts： 按 chromStart 排序的 exon 起始位置，所有 exon 的 blockStarts 均为与 chromStart 的相对位置（即原 bed 文件中对应 exon 的 chromStart - 第一条 exon 的 chromStart）。

举个例子：

| raw bed | | | | | | | | | | | |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| test | 2 | 30 | geneA_1 | 255 | - | 2 | 30 | 255，0，0 | 1 | 28 | 0 |
| test | 102 | 132 | geneA_2 | 255 | - | 102 | 132 | 255，0，0 | 1 | 30 | 0 |

| updated bed | | | | | | | | | | | |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| test | 2 | 132 | geneA | 255 | - | 2 | 132 | 255，0，0 | 2 | 28,30 | 0，100 |


上述两种情况在更新好 bed文件 之后，需添加 --update 参数再次运行脚本以生成正确的 gtf。 


## reference

1.https://github.com/alexdobin/STAR

2.https://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format