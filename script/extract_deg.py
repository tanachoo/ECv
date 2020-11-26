# Author: yoshi
# Date: 12/13/2019
# Updated: 12/17/2019
# Project: ECv
# Script: extract differential express genes following PDK4 papers

import pandas as pd
import math

def extract_deg(filename):
    ## Set DEG cutoff
    FC_cutoff = 2
    q_value_cutoff = 0.00001    

    print(f'load: {filename}\n'
          f'FC cutoff: {FC_cutoff}\n'
          f'q_value_cutoff: {q_value_cutoff}')

    ## Load data
    degdata = pd.read_table(filename, sep='\t', header=0, index_col=0)

    ## extract gene using cutoff ##
    ## using FC
    degdata_FC_cutoff = degdata[abs(degdata['log2FC'].astype(float)) > math.log2(FC_cutoff)]
    FC_genes = list(degdata_FC_cutoff.index) 

    ## using q-value
    degdata_qcutoff = degdata[degdata['q_value'].astype(float) < q_value_cutoff]
    qvalue_genes = list(degdata_qcutoff.index)

    ## acquire shared gene between FC and q-value
    FC_qvalue_shared_genes = list(set(FC_genes) & set(qvalue_genes))

    return FC_genes, qvalue_genes, FC_qvalue_shared_genes


def ECv_gene():
    """ aquire ΔECv(th=1) gene """
    ecv_gene = []

    with open('../../../../ownCloud2/MicroarrayAnalysis_BN/microarray_data/GSE49644/main_analysis/ECv_GSE49644_deltaECv/gene_forVenn/3shared_genes_post_Benn_pairs/shared_3celllines_120pair_150genes_genelist.txt', 'r') as f:
        for i in f:
            gene = i.strip()        
            ecv_gene.append(gene)

    return ecv_gene


def writeout(list_, filename):
    """ Write list as txt file """
    with open(filename, 'wt') as f:
        for i in list_:
            f.write(i+'\n')


if __name__ == '__main__':

    ## DEG for HCC827 
    FC_h, qvalue_h, FC_qvalue_shared_h  = extract_deg('./GSE49644_DEG_HCC827.txt')
    print(f'#FC_genes: {len(FC_h)}\n'
          f'#qvalue_gene/Users/yoshi/Desktop/microarray_data/GSE49644/data : {len(qvalue_h)}\n'
          f'#FC_qvalue_shared_gene: {len(FC_qvalue_shared_h)}\n')

    ## DEG for A549
    FC_a, qvalue_a, FC_qvalue_shared_a  = extract_deg('./GSE49644_DEG_A549.txt')
    print(f'#FC_genes: {len(FC_a)}\n'
          f'#qvalue_gene: {len(qvalue_a)}\n'
          f'#FC_qvalue_shared_gene: {len(FC_qvalue_shared_a)}\n')

    ## DEG for NCI-H358
    FC_n, qvalue_n, FC_qvalue_shared_n  = extract_deg('./GSE49644_DEG_NCIH358.txt')
    print(f'#FC_genes: {len(FC_n)}\n'
          f'#qvalue_gene: {len(qvalue_n)}\n'
          f'#FC_qvalue_shared_gene: {len(FC_qvalue_shared_n)}\n')

    ## extract genes
    share_h_a_n = list(set(FC_qvalue_shared_h) & set(FC_qvalue_shared_a) & set(FC_qvalue_shared_n))
    share_h_a = list(set(FC_qvalue_shared_h) & set(FC_qvalue_shared_a))
    share_h_or_a = list(set(FC_qvalue_shared_h) | set(FC_qvalue_shared_a))
    share_h_n = list(set(FC_qvalue_shared_h) & set(FC_qvalue_shared_n))
    share_h_or_n = list(set(FC_qvalue_shared_h) | set(FC_qvalue_shared_n))
    share_a_n = list(set(FC_qvalue_shared_a) & set(FC_qvalue_shared_n))
    share_a_or_n = list(set(FC_qvalue_shared_a) | set(FC_qvalue_shared_n))
    share_han_or = list(set(set(share_h_a) | set(share_h_n) | set(share_a_n)))
    print(f'#share_h_a_n: {len(share_h_a_n)}\n'
          f'#share_h_a: {len(share_h_a)}\n'
          f'#share_h_or_a: {len(share_h_or_a)}\n'
          f'#share_h_n: {len(share_h_n)}\n'
          f'#share_h_or_n: {len(share_h_or_n)}\n'
          f'#share_a_n: {len(share_a_n)}\n'
          f'#share_a_or_n: {len(share_a_or_n)}\n'
          f'#share_han_or: {len(share_han_or)}\n')
    #print(share_han_or)

    ecv_gene = ECv_gene()
    print(f'#ECv_gene: {len(ecv_gene)}\n')

    share_deg_ecv = list(set(set(ecv_gene) & set(share_h_a_n)))
    print(f'#share_deg_ecv: {len(share_deg_ecv)}\n')


    writeout(share_deg_ecv, './DEG125_3cellshareFC2q000001_deltaECv150_share_genes.txt')
