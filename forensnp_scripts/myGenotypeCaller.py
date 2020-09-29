
"""

## Create: Zilin Ren
## Date:   2020-09-02
## Email:  Zilin.Ren@outlook.com

## Workflow:
    1. output: genotype and coverage

"""

import numpy as np
import pandas as pd 


def func_cal_coverage(x):
    coverage = sum([x["A"], x["G"], x["C"], x["T"]])
    return coverage


def func_annotate_rs_into_ntfreq(ntfreq_file, rs_ref_path):
    
    ## search for each line
    def search_rs_num(x, rs_ref):
        rs_arr = rs_ref.loc[(rs_ref['#chr'] == x['#chr']) 
                            & (rs_ref['start'] <= int(x['pos']))
                            & (rs_ref['end'] >= int(x['pos'])), 'rs'].values
        return rs_arr[0]
    
    ## 
    ntfreq_dat = pd.read_csv(ntfreq_file, sep="\t")
    rs_ref = pd.read_csv(rs_ref_path, sep="\t", header=None)
    rs_ref.columns = ['#chr', 'start', 'end', 'rs', 'allele']
    ntfreq_dat['rs'] = ntfreq_dat.apply(lambda x: search_rs_num(x, rs_ref), axis=1)
    ntfreq_dat['ref'] = ntfreq_dat['ref'].apply(lambda x: x.upper())
    return ntfreq_dat


def get_allele_each_row(x, homo_cutoff, cov_cutoff=10):

    coverage= x['cov']

    if coverage == 0:
        return [pd.NA, pd.NA, coverage, 'no-coverage'] + [0 for i in range(8)]

    forware_reverse = x["A5'/3'"].split('/') + x['G5/3'].split('/') + x['C5/3'].split('/') + x['T5/3'].split('/')

    ratio_A = float(x['A']) / coverage
    ratio_G = float(x['G']) / coverage
    ratio_C = float(x['C']) / coverage
    ratio_T = float(x['T']) / coverage

    freq_df = pd.DataFrame([['A', ratio_A], ['G', ratio_G], ['C', ratio_C], ['T', ratio_T]], columns=['allele', 'ratio'])
    INFO    = "".join(["{}:{:.2%};".format(i[0], i[1]) for i in freq_df.values.tolist()])
    freq_df.sort_values(by = 'ratio', ascending = False, inplace=True)

    if coverage < cov_cutoff:
        return [pd.NA, pd.NA, coverage, INFO] + forware_reverse

    if freq_df.loc[freq_df['ratio'] > homo_cutoff, :].shape[0] == 1:
        allele1 = freq_df['allele'].values[0]
        allele2 = allele1

    elif freq_df.loc[freq_df['ratio'] > homo_cutoff, :].shape[0] == 2:
        allele1 = freq_df['allele'].values[0]
        allele2 = freq_df['allele'].values[1]

    else:
        allele1 = '?' + freq_df['allele'].values[0]
        allele2 = '?' + freq_df['allele'].values[1]

    return [allele1, allele2, coverage, INFO] + forware_reverse



def get_alleles(sampleID, output_ntfreq, snp_bed_path, output_genotype, homo_cutoff):



    
    ntfreq_df = func_annotate_rs_into_ntfreq(output_ntfreq, snp_bed_path)
    ntfreq_df['cov'] = ntfreq_df.apply(lambda x: func_cal_coverage(x), axis=1)
    
    
    alleles_series = ntfreq_df.apply(lambda x: get_allele_each_row(x, homo_cutoff), axis=1)
    alleles_df = pd.DataFrame(alleles_series.tolist(), columns=['allele1', 'allele2', 'totalCoverage', 'INFOs', 'A+', 'A-', 'G+', 'G-', 'C+', 'C-', 'T+', 'T-'])

    raw_output_df = pd.concat([ntfreq_df, alleles_df], axis=1)
    output_df = raw_output_df.loc[:, ['#chr', 'pos', 'rs', 'ref', 'allele1', 'allele2', 'totalCoverage', 'INFOs', 'A+', 'A-', 'G+', 'G-', 'C+', 'C-', 'T+', 'T-']]
    output_df['SampleID'] = sampleID
    output_df.to_csv(output_genotype, index=None, na_rep="NA")
    return 