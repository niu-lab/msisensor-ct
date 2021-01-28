# MSIsensor-ct

MSIsensor-ct is a novel method based machine learning, which can accurately determine the MSI status of cfDNA sequencing data with ultra-low ctDNA content. In practice, we recommend aiming for a minimum of 30 valid microsatellites to classify the MSI status using MSIsensor-ct. MSIsensor (https://github.com/niu-lab/msisensor) is a viable option if paired tumor and normal sequencing data are available, and MSIsensor2 (https://github.com/niu-lab/msisensor2) is effective for tumor-only samples.

## Citation
Xinyin Han, Shuying Zhang, Daniel Cui Zhou, Dongliang Wang, Xiaoyu He, Danyang Yuan, Ruilin Li, Jiayin He, Xiaohong Duan, Michael C Wendl, Li Ding, Beifang Niu, MSIsensor-ct: microsatellite instability detection using cfDNA sequencing data, Briefings in Bioinformatics, 2021. 




## Advantages

1. **High sensitivity and specificity**

   For detecting MSI status in cfDNA sequencing data, MSIsensor-ct reached 100% sensitivity and specificity in our genuine samples and 17 simulation datasets.

2. **A stable threshold**

   According to the robustness test, if a panel contains 30 microsatellites overlapped with 1,476 site-classifiers in MSIsensor-ct, the AUC for MSI calling can attain up to 0.99. In addition, no matter how many sites overlapped with our classifiers, MSIscore = 20% can always be a potential stable threshold to distinguish MSI and MSS samples.

3. **Suitable for extra-low ctDNA content** 

   The limitation test demonstrated that MSIsensor-ct accurately detected the MSI status on cfDNA samples with the sequencing depth of 3,000× and ctDNA content at the level of 0.05%.

4. **User-friendly installation and operation process** 

   MSIsensor-ct requires only BAM files as input and is free from additional baseline establishment, which is user-friendly and can be flexibly integrated into the routine next generation sequencing analysis.



## Install

Currently, MSIsensor-ct is based on Linux system, and we provide binaries only. Please note your GCC version should be at least 5.0.x.

```
    git clone https://github.com/niu-lab/msisensor-ct.git
    cd msisensor-ct
    chmod +x msisensor-ct
```



## Usage

```
    Version 0.1
    Usage:  msisensor-ct <command> [options]
```

Key commands:

```
    msi            msi scoring
```

msisensor-ct msi [options]:

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


hg38 bam:

```
    msisensor-ct msi -D -M ./models_hg38 -t ./test/example.cfdna.hg38.bam -o output.prefix
```

hg19 or GRCh37 bam:

```
    msisensor-ct msi -D -M ./models_hg19_GRCh37 -t ./test/example.cfdna.hg19.bam -o output.prefix
```

b37 or HumanG1Kv37 bam:

```
    msisensor-ct msi -D -M ./models_b37_HumanG1Kv37 -t ./test/example.cfdna.b37.bam -o output.prefix
```

Note: bam index files are needed in the same directory as bam files



## Output

The MSI scoring step produces 3 files:

```
    output.prefix
    output.prefix_dis
    output.prefix_somatic
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

We provided one small dataset to test the msi scoring step:

```
    msisensor-ct msi -D -M ./models_hg38 -t ./test/example.cfdna.hg38.bam -o output.prefix
    msisensor-ct msi -D -M ./models_hg19_GRCh37 -t ./test/example.cfdna.hg19.bam -o output.prefix
    msisensor-ct msi -D -M ./models_b37_HumanG1Kv37 -t ./test/example.cfdna.b37.bam -o output.prefix
```

We also provided a R script to visualize MSI score distribution of MSIsensor-ct output. ( msi score list only or msi score list accompanied with known msi status). 

For msi score list only as input:

```
    R CMD BATCH "--args msi_score_only_list msi_score_only_distribution.pdf" plot.r
```

For msi score list accompanied with known msi status as input:

```
    R CMD BATCH "--args msi_score_and_status_list msi_score_and_status_distribution.pdf" plot.r
```



## Contact

If you have any questions, please contact one or more of the following folks: Beifang Niu [bniu@sccas.cn](mailto:bniu@sccas.cn); Xinyin Han [hanxinyin@cnic.cn](mailto:hanxinyin@cnic.cn); Shuying Zhang [zhangshuying@cnic.cn](mailto:zhangshuying@cnic.cn).
