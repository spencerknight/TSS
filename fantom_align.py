import numpy as np
import pandas as pd
from tqdm import tqdm, tqdm_pandas

chromosomes = ['chr{}'.format(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
chromosomes = sorted(chromosomes)
columns = ['PMA Control', 'PMA 1hr', 'PMA 4hr', 'PMA 12hr', 'PMA 24hr', 'PMA 96hr']

def is_subset(promoter_index, start, end):
    '''
    Determines if promoter is subset of index range
    '''
    if promoter_index > start and promoter_index < end:
        return True
    else:
        return False

def promoter_constructor(column, gene_dataframe, cage_dataframe):
    '''
    Aligns promoters from cage reads to gene regions
    #param column: the sample of interest (e.g., "PMA 1hr")
    #param gene_dataframe: the gene region dataframe from UCSC
    #param cage_dataframe: the CAGE alignment dataframe from UCSC
    '''
    #select the dataframe according to the appropriate column
    cg2 = cage_dataframe[cage_dataframe['rna'] == column]
    #sort the data frame
    gene_dataframe = gene_dataframe.sort_values(by = 'seqname', ascending = True)
    #initialize promoter dictionary list that has index and read info for each promoter
    promoter_dict_ls = []
    #iterate through each chromosome
    for chrom in tqdm(chromosomes):
        #isolate starts and ends of gene regions as well as strand
        gfx = gene_dataframe[gene_dataframe['seqname'] == chrom]
        gene_starts = gfx['start'].astype(int).tolist()
        gene_ends = gfx['end'].astype(int).tolist()
        gene_signs = gfx['strand'].tolist()
        cfx = cg2[cg2['chrom'] == chrom]
        cage_starts = cfx['start'].astype(int).tolist()
        cage_ends = cfx['end'].astype(int).tolist()
        cage_strands = cfx['strand'].tolist()
        cage_reads = cfx['rna_count'].astype(int).tolist()
        #Loop through all of the gene starts
        for i in tqdm(range(len(gene_starts))):
            #initialize start and end positions as well as the strand
            s = gene_starts[i]
            e = gene_ends[i]
            st = gene_signs[i]
            #initialize the tss and read values as well as promoter dictionary
            indices = []
            reads = []
            promoter_dict = {}
            #isolate the cage_starts that are in range
            in_range_cage = list(set([c for c in cage_starts if is_subset(c, s, e) == True]))
            #isolate the in-range indices
            in_range_indices = []
            for c in in_range_cage:
                in_range_indices += [k for k,x in enumerate(cage_starts) if x == c]
            #use the indices to map the in range starts, reads, and strands
            in_range_start = [cage_starts[k] for k in in_range_indices]
            in_range_end = [cage_ends[k] for k in in_range_indices]
            in_range_read = [cage_reads[k] for k in in_range_indices]
            in_range_strands = [cage_strands[k] for k in in_range_indices]
            #Loop through candidates and extract the ones that are in frame
            for j in range(len(in_range_start)):
                if in_range_strands[j] == st:
                    if st == '+':
                        indices.append(in_range_start[j]-s-1)
                        reads.append(in_range_read[j])
                    elif st == '-':
                        indices.append(e-in_range_end[j])
                        reads.append(in_range_read[j])
                else:
                    continue
            for k,ind in enumerate(indices):
                if ind in promoter_dict:
                    promoter_dict[ind] += reads[k]
                else:
                    promoter_dict[ind] = reads[k]
            promoter_dict_ls.append(promoter_dict)
    gene_dataframe['{}'.format(column)] = promoter_dict_ls
    return gene_dataframe

if __name__ == "__main__":
    print 'LOADING DATA...\n'
    genes = pd.read_pickle('ucsc_gene_regions.pkl')
    fantom = pd.read_pickle('i03_mapping.pkl')
    for column in columns:
        print 'NOW PROCESSING DATA FOR COLUMN: {}'.format(column)
        genes = promoter_constructor(column, genes, fantom)
    genes.to_pickle('ucsc_genes_with_reads.pkl')
    print 'DONE!'
