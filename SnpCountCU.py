#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
"This script is used to count the number of SNPs that are common among 
population and unique within population from VCF formart file."
'''

__author__      = "Jingfang SI"
__version__     = "1.0"
__email__       = "sijingfang@foxmail.com"
__date__        = "2020-06-06"


import sys
import os
import time
import psutil
import argparse
import numpy as np
from multiprocessing import Pool
from functools import wraps
from cyvcf2 import VCF

def timethis(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        r = func(*args, **kwargs)
        end = time.perf_counter()
        print('[{}] INFO: Done {}.{}, time used:{:.4f}'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),func.__module__, func.__name__, end - start))
        return r
    return wrapper

def get_region_list(region_file):
    region_list = []
    with open(region_file, "r") as fp:
        for line in fp:
            region_list.append(line.strip())
    return region_list

## get the list of populations and link the samples to their populations
def get_list_check(in_file, in_list):
    vcf = VCF(in_file)
    ## get all the sample ids in the vcf
    samples_list = vcf.samples
    dict_check = {}
    for sample_id in samples_list:
        dict_check[sample_id] = True

    dict_id_pop = {}
    pop_list = []

    with open(in_list, "r") as fp:
        for line in fp:
            pop = line.strip().split()[0]
            sample = line.strip().split()[1]
            ## check whether the sample id in list exists in vcf file
            if sample in dict_check:
                dict_id_pop[sample] = pop
            else:
                sys.exit("ERROR: Sample ID {0} is not exsit in vcf.".format(sample))

            if pop not in pop_list:
                pop_list.append(pop)
    return [pop_list, dict_id_pop]


def count(allele_list):
    '''gt_types is a array of 0,1,2,3, represent HOM_REF, HET, UNKNOWN, HOM_ALT'''
    allele_array = np.array(allele_list)
    count_size = np.size(allele_array)
    # count_HOM_REF = np.size(allele_array[allele_array == 0])
    count_HET = np.size(allele_array[allele_array == 1])
    # count_UNKNOWN = np.size(allele_array[allele_array == 2])
    count_HOM_ALT = np.size(allele_array[allele_array == 3])
    
    count_ALL = 2 * count_size
    count_ALT = 2 * count_HOM_ALT + count_HET
    freq_ALT = count_ALT / count_ALL
    return freq_ALT

def count_snp_index(count_out, pop, pop_index, feq_threshold):
    '''
    using number 0,1,2,3,4 ... repersent snp id in count_out, extract a list of snps of which ALT allele frequence greater than the threshold
    '''
    threshold = feq_threshold
    i = 0
    index_list = []
    for ferq_1 in count_out[pop_index]:
        if ferq_1 >= threshold:
            index_list.append(i) 
        else:
            pass
        i += 1
    return index_list


@timethis
def vcf_io(vcf, chrom, pop, dict_id_pop):
    print('[{0}] INFO: Start count population: {1}'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), pop))
    vcf = VCF(in_file)
    samples_list = vcf.samples

    ## select samples in the specific population and get the index of those samples
    samples_list_index = []
    for sample in samples_list:
        if sample in dict_id_pop and dict_id_pop[sample] == pop:
            samples_list_index.append(samples_list.index(sample))
    
    if(len(samples_list_index) == 0):
        sys.exit("ERROR: Population {0} is empty.".format(pop))
    else:
        ## extract genotype recode of each sample depned on sample index
        output = []
        for line in vcf(chrom):
            record_GT = [line.gt_types[i] for i in samples_list_index]
            ## count
            output.append(count(record_GT))
        
        process = psutil.Process(os.getpid())
        print("[{}] INFO: PID:{}  MEM:{:.2f}MB".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), process.pid, process.memory_info().rss / 1000000))
        return output


@timethis
def count_snp_common(pop_target, list_pop_all, count_snp_index_box):
    '''
    Calculate intersection of all populations as the common snp
    '''
    print('[{0}] INFO: Find common snps among populations.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    set_index = set(count_snp_index_box[list_pop_all.index(pop_target[0])])
    for i in pop_target:
        set_index = set_index & set(count_snp_index_box[list_pop_all.index(i)])
    return len(set_index)

@timethis
def count_snp_unique(pop_target, list_pop_all, count_snp_index_box):
    '''
    Calculate the difference between the target group and other groups as the population unique snp
    '''
    print('[{0}] INFO: Find unique snp in population: {1}'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), pop_target))
    set_index_others = set()
    for pop_x in list_pop_all:
        if pop_x != pop_target:
            set_index_others = set_index_others | set(count_snp_index_box[list_pop_all.index(pop_x)])
        else:
            set_index_target = set(count_snp_index_box[list_pop_all.index(pop_target)])
    return len(set_index_target - set_index_others)


def use_pool(func, args_list, num_threads):
    with Pool(num_threads) as pool:
        out = pool.starmap(func, args_list)
        pool.close()
        pool.join()
    return out

def final_run(in_file, chrom, list_pop, dict_id_pop, common_pop, out_file):
    process = psutil.Process(os.getpid())
    print("[{}] INFO: ========================================================".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    print("[{}] INFO: PID:{}  MEM:{:.2f}MB".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), process.pid, process.memory_info().rss / 1000000))
    print('[{0}] INFO: Chromosome {1}'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), chrom))
    print("[{}] INFO: ========================================================".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    ## calculate the count and the ferquence of the alt snp
    runlist_vcf_io = [(in_file, chrom, pop, dict_id_pop) for pop in list_pop]
    count_out = use_pool(vcf_io, runlist_vcf_io, num_threads)
    
    runlist_count_snp_index = [(count_out, pop, list_pop.index(pop), feq_threshold) for pop in list_pop]
    count_snp_index_box = use_pool(count_snp_index, runlist_count_snp_index, num_threads)
    count_out = None
    
    ## get common snp count
    '''
    Finding the common snp among several poplutions. 
    The command specify the combined group ID and poplutions including in the group. 
    e.g. GroupA:pop1,pop2,pop3;GroupB:pop4,pop5
    '''
    list_pop_common = [list_pop]
    if common_pop != None:
        for group_pop in common_pop.split(";"):
            # group = group_pop.split(":")[0]
            pops = group_pop.split(":")[1].split(",")
            list_pop_common.append(pops)

    runlist_count_snp_common = [(pop_target, list_pop, count_snp_index_box) for pop_target in list_pop_common]
    count_snp_common_out = use_pool(count_snp_common, runlist_count_snp_common, num_threads)
    
    ## get unique snp count
    runlist_count_snp_unique = [(pop_target, list_pop, count_snp_index_box) for pop_target in list_pop]
    count_snp_unique_out = use_pool(count_snp_unique, runlist_count_snp_unique, num_threads)
    
    ## write out the result of common/unique snp count
    with open(out_file, "a") as wp:
        out = [chrom] + [str(i) for i in count_snp_common_out] + [str(i) for i in count_snp_unique_out]
        wp.write("\t".join(out) + "\n")
    out = None
    count_snp_index_box = None
    count_snp_common_out = None
    count_snp_unique_out = None
    print("[{}] INFO: ========================================================".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    print("[{}] INFO: PID:{}  MEM:{:.2f}MB".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), process.pid, process.memory_info().rss / 1000000))
    print("[{}] INFO: ========================================================".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    print('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, description = __doc__)
    parser.add_argument("-v", "--vcf", required = True, 
                        help = "Name of the input vcf file, must be gziped and indexed by 'bcftools index'")
    parser.add_argument("-l", "--list", required = True, 
                        help = "Name of the input population list, two columns: pop_id sample_id")
    parser.add_argument("-r","--region", required = True, 
                        help = "Name the input list file, one chromosome name per line")
    parser.add_argument("-o", "--out", required = True, 
                        help = "Name of the output file")
    parser.add_argument("-f", "--freq-threshold", required = False, type = float, default = 0.00001,
                        help="A frequence threshold value, SNPs with ALT allele frequence greater than this value are regarded as existing in a population.(default: 0.00001)")
    parser.add_argument("--common-pop", required = False, type = str, default = None,
                        help="Add additional population combinations to count common snp, e.g. 'GroupA:pop1,pop2,pop3;GroupB:pop4,pop5'")
    parser.add_argument("-nt", "--num-threads", required = False, type = int, default = 1,
                        help = "Number of threads. This value should be lower than the number of populations in the list, otherwise it will not provide additional efficiency")
    args = parser.parse_args()

    in_file = args.vcf
    in_list = args.list

    region_file = args.region
    num_threads = args.num_threads
    feq_threshold = args.freq_threshold

    ## get the list of populations
    get_list_check_out = get_list_check(in_file, in_list)
    list_pop = get_list_check_out[0]
    dict_id_pop = get_list_check_out[1]
    
    chrom_list = get_region_list(region_file)
    
    print('[{0}] INFO: Start SnpCount. Total populations : {1}.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), len(list_pop)))
    print('[{0}] INFO: Start SnpCount. Total chromosome  : {1}.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), len(chrom_list)))
    ## write the header of output file
    with open(args.out, "w") as wp:
        wp.write("\t".join(["CHROM", "COMMOM"]) + "\t")
        if args.common_pop != None:
            wp.write("\t".join([i.split(":")[0] for i in args.common_pop.split(";")]) + "\t")
        wp.write("\t".join([i for i in list_pop]) + "\n")
    
    ## start the loop of chromosomes
    for chrom in chrom_list:
        final_run(in_file, chrom, list_pop, dict_id_pop, args.common_pop, args.out)

    print( '[{0}] INFO: Done'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())) )

