# SnpCountCU
Count common and unique SNPs among several populations from a VCF format file.

## _Brief description_
This script is used to count the number of SNPs that are common among populations and unique within populations from VCF formart file. 

The populations and individuals included in the population are specified by the `-l/--list` parameter. The number of population-shared and population-specific SNPs are calculated based on the intersection and difference of the snp contained in each population. By default, SNPs contained in a population are defined as all non-"0/0" type sites of the individual corresponding to the population in the vcf file. You can set the `-f/--freq-threshold` parameter to filter SNPs in each population according to the frequency of ALT alleles. SNPs shared among all population and unique to each population determined in the list file are counted by default. You can use the `--common-pop` parameter to determine additional population combinations to calculate the number of snp they shared. 

Additionally, you can uses multiple threads to perform parallel calculations on each population to speed up the program. Therefore, the number of threads should be less than the set total number of populations, otherwise it will not provide higher calculation efficiency.

Note: This script is based on python3. Before running this script, make sure that your python version meets the requirements, and the  modules required in the script have been installed, especially [cyvcf2](https://github.com/brentp/cyvcf2).

Please don't hesitate to open an [`Issue`](https://github.com/JingfangSI/SnpCountCU/issues) if you find any problem or suggestions for a new feature.

## _Usage_
Just type `python3 SnpCountCU.py -h` or `./SnpCountCU.py -h` to show the help of the program:
```
usage: SnpCountCU.py [-h] -v VCF -l LIST -r REGION -o OUT [-f FREQ_THRESHOLD]
                     [--common-pop COMMON_POP] [-nt NUM_THREADS]

This script is used to count the number of SNPs that are common among population and unique within population from VCF formart file.

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Name of the input vcf file, must be gziped and indexed
                        by 'bcftools index'
  -l LIST, --list LIST  Name of the input population list, two columns: pop_id
                        sample_id
  -r REGION, --region REGION
                        Name the input list file, one chromosome name per line
  -o OUT, --out OUT     Name of the output file
  -f FREQ_THRESHOLD, --freq-threshold FREQ_THRESHOLD
                        A frequence threshold value, SNPs with non-REF
                        genotype frequence greater than this value are
                        regarded as existing in a population.(default:
                        0.00001)
  --common-pop COMMON_POP
                        Add additional population combinations to count common
                        snp, e.g. 'GroupA:pop1,pop2,pop3;GroupB:pop4,pop5'
  -nt NUM_THREADS, --num-threads NUM_THREADS
                        Number of threads. This value should be lower than the
                        number of population in the list, otherwise it will
                        not provide additional efficiency
```

## _Examples_
In the following examples you can omit `python3` if you change the permissions of `vcf2phylip.py` to executable.

_Example 1:_ Use default parameters and 4 threads:
```bash
python3 SnpCountCU.py -v file.vcf.gz -l pop.list  -r chrom.list -o snpcout_common_uniq.out -nt 4
```
_Example 2:_ Add two new population combinations(G1 and G2) and use frequence threshold 0.8:
```bash
python3 SnpCountCU.py -v file.vcf.gz -l pop.list  -r chrom.list -f 0.8 --commom-pop G1:pop1,pop2;G2:pop3,pop4 -o snpcout_common_uniq.out -nt 4
```
