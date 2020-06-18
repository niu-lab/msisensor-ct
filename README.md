# MSIsensor-blood

MSIsensor-blood is a novel method based machine learning, which can accurately
determine the MSI status of cfDNA samples with ultra-low ctDNA content.



## Advantages

1. **High sensitivity and specificity**

   In detecting MSI status in cfDNA sequencing samples, MSIsensor-blood can obtain high sensitivity and specificity in data with extremely low ctDNA content and sequencing depths.

2. **A stable threshold**

   According to the robustness test, if a panel contains 30 microsatellites overlapped with 1476 site-classifiers, the AUC for MSI calling can up to 0.99.  Besides, no matter how many sites a panel overlaps with our classifiers,  MSIsensor-blood alwasys shows a stable cutoff that MSIscore=20%.

3. **Suitable for ultra-low ctDNA frequency** 

   According to the limitation test, MSIsensor-blood can show accurate MSI calling ability on cfDNA samples with sequencing depth of 10,000X and ctDNA content of 0.05%.

4. **User-friendly installation and operation process** 

   MSIsensor-blood only requires BAM files as input, which is user-friendly and can be flexibly integrated into the routine next generation sequencing (NGS) analysis.



##Install

Currently, MSIsensor-blood is based on Linux system, and we provide binaries only. Please note your GCC version should be at least 5.0.x.

```
    git clone https://github.com/niu-lab/msisensor-blood.git
    cd msisensor-blood
    chmod +x msisensor-blood
```



## Usage

```
    Version 0.1
    Usage:  msisensor-blood <command> [options]
```

Key commands:

```
    msi            msi scoring
```

msisensor-blood msi [options]:

```
   -D   <boolean>  activate processing for ctDNA samples
   -M   <string>   models directory for tumor only data
   -t   <string>   tumor bam file
   -o   <string>   output distribution file

   -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default=20
   -b   <int>      threads number for parallel computing, default=1
   -x   <int>      output homopolymer only, 0: no; 1: yes, default=0
   -y   <int>      output microsatellite only, 0: no; 1: yes, default=0

   -h   help
```



## Example

MSI scoring:

for tumor only sequence data:

hg38 bam:

```
    msisensor-blood msi -D -M ./models_hg38 -t ./test/example.tumor.only.hg38.bam -o output.tumor.prefix
```

hg19 or GRCh37 bam:

```
    msisensor-blood msi -D -M ./models_hg19_GRCh37 -t ./test/example.tumor.only.hg19.bam -o output.tumor.prefix
```

b37 or humanG1Kv37 bam:

```
    msisensor-blood msi -D -M ./models_b37_HumanG1Kv37 -t ./test/example.tumor.only.b37.bam -o output.tumor.prefix
```

Note: bam index files are needed in the same directory as bam files



## Output

For tumor only input, the MSI scoreing step produces 3 files:

```
    output.tumor.prefix
    output.tumor.prefix_dis
    output.tumor.prefix_somatic
```

1. output.prefix: msi score output

   ```
    Total_Number_of_Sites   Number_of_Somatic_Sites %
    2     2      100.00
   ```

2. output.prefix_dis: read count distribution (T: tumor)

   ```
    chr22 29286892 AAAGC 12[T] CTCTT
    T: 0 0 0 0 0 0 0 0 25 71 4 86 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
   ```

3. output.prefix_somatic: somatic sites detected

   ```
    chromosome   location        left_flank     repeat_times    repeat_unit_bases    right_flank    discrimination_value_ML
    chr22	29286892	AAAGC	12	T	CTCTT	0.98852
   ```



## Test sample

We provided one small dataset (tumor only bam file) to test the msi scoring step:

```
    msisensor-blood msi -D -M ./models_hg38 -t ./test/example.tumor.only.hg38.bam -o output.tumor.prefix
    msisensor-blood msi -D -M ./models_hg19_GRCh37 -t ./test/example.tumor.only.hg19.bam -o output.tumor.prefix
    msisensor-blood msi -D -M ./models_b37_HumanG1Kv37 -t ./test/example.tumor.only.b37.bam -o output.tumor.prefix
```

We also provided a R script to visualize MSI score distribution of MSIsensor-blood output. ( msi score list only or msi score list accompanied with known msi status). 

For msi score list only as input:

```
    R CMD BATCH "--args msi_score_only_list msi_score_only_distribution.pdf" plot.r
```

For msi score list accompanied with known msi status as input:

```
    R CMD BATCH "--args msi_score_and_status_list msi_score_and_status_distribution.pdf" plot.r
```



## Contact

If you have any questions, please contact one or more of the following folks: Beifang Niu [bniu@sccas.cn](mailto:bniu@sccas.cn); Xinyin Han [hanxinyin@cnic.cn ](mailto:hanxinyin@cnic.cn); Shuying Zhang zhangshuying@cnic.cn.
