使用说明-中文版
=================

- [使用说明-中文版](#使用说明-中文版)
  - [版本更新](#版本更新)
    - [1.01版](#101版)
  - [1. 联系方式](#1-联系方式)
  - [2. 安装](#2-安装)
  - [3. 使用说明](#3功能模块)
     - [（1）参数详解](#1参数详解)
     - [（2）输入文件](#2输入文件)
     - [（3）输出文件](#3输出文件)
     - [（4）示例代码](#4示例代码)

## 版本更新

### 1.03版

时间：2023年04月03日

说明：基于用户提供的区间标记数或者每个区间包含的标记数文件运行MB-BayesAS模型

## 1. 联系方式

liwn@cau.edu.cn

## 2. 安装

下载源文件用gcc编译，或者下载编译版本直接运行
依赖：需要调用mkl库

## 3.使用说明

#### （1）参数详解

```text
 Required parameters:
  -g, --raw=FILE             plink raw file, e.g. /path/plink.raw
  -l, --bfile=FILE           plink binary file prefix, e.g. /path/plink

 Optional parameters:
  -a, --id=INT               id position at phenotype file [1]
  -A, --varOut=FILE          posteriori mean of dispersion matrix [var.txt]
  -b, --binf=FILE            file contains number of SNPs per interval [NULL]
  -B, --burnin=INT           burnin steps [100]
  -c, --varC=INT             output first c predictions in MCMC chain [NULL]
  -d, --dfvara=DOUBLE        degrees of freedom Va [4]
  -D, --dfvare=DOUBLE        degrees of freedom Ve [4]
  -e, --region=CHR           Method for dividing snp into different groups,
                             LD/fix [LD]
  -E, --effOut=FILE          posteriori mean of location effect [theta.txt]
  -f, --mapf=FILE            plink .map (.bim) use to define bins [NULL]
  -F, --fix=CHR              fixed effect position, e.g. '2 3' [NULL]
  -G, --diagG=INT            a integer large enough [100000]
  -i, --prioSa=CHR           distribution of scale of additive effect prior
  -I, --prioSe=CHR           distribution of scale of residual effect prior
  -j, --bincol=INT           position of nSNPs in each bins at binf [1]
  -L, --logf=FILE            redirect standard output to this file
  -m, --miss=DOUBLE          missing value in phenotype [-99]
  -M, --mcmcOut=FILE         output predictions in MCMC chain [MCMC_var.txt]
  -N, --nsnp=FILE            Number of SNPs per interval [100]
  -o, --outdir=CHR           output directory [.]
  -p, --phef=FILE            phenotype file
  -r, --report=INT           MCMC sampling report interval [100]
  -R, --ran=CHR              environmental random effect position [NULL]
  -s, --seed=LONG            random number seed [NULL]
  -S, --iter=INT             length of MCMC chain [1000]
  -t, --thin=INT             thinning rate [10]
  -v, --varf=FILE            SNP variance prior [varp * 0.5]
  -V, --gebvOut=FILE         breeding value output file [gebv.txt]
  -y, --phe=CHR              phenotypes position, e.g. '4 5' [NULL]

 Flag parameters:
  -C, --dontcenter           centering gene content matrix
  -H, --header               phenotype and bins files has header, [automatic
                             detection]

  -?, --help                 Give this help list
      --usage                Give a short usage message
      --version              Print program version
```

#### （2）输入文件  

基因型文件：plink raw 格式(*.raw)文件，具体参考[PLINK文件格式](https://www.cog-genomics.org/plink/1.9/formats#raw)说明

表型文件：以空格为分隔符的纯文本文件，包含个体ID、固定效应和表型值，固定效应可以用字符格式表示，字符格式的值最大长度为20字符，文件中个体顺序应与与raw中第一列ID对应

#### （3）输出文件

见参数说明

#### （4）示例代码

BayesAS --rawf test.raw --phef test.txt --nsnp 100
