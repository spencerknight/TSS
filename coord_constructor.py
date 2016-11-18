import pandas as pd
from Bio import SeqIO
from tqdm import tqdm, tqdm_pandas

tqdm.pandas()

genes = pd.read_pickle('ucsc_genes.pkl')
genes = genes.dropna(subset = ['start', 'end'])

def new_coord(old_start, old_end, sign):
    '''
    Calculates the new region upstream of the initiation codon to scan for promoters
    #param old_start: old start sequence
    #param old_end: old end sequence
    #sign: + or - strand
    '''
    if sign == '+':
        return [old_start - 2001, old_end]
    #note: indexing is different for start site in database (reason for -2001 instead of 2000)
    elif sign == '-':
        return [old_start, old_end + 2000]

def reverse_complement(seq):
    '''
    #param seq: sequence to be reverse complemented
    '''
    seq = seq.upper()[::-1]
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N': 'N'}
    return ''.join([complement[base] for base in list(seq)])

def sequencer(seq, sign):
    '''
    Determines if you need to apply the (+) or reverse complement of the sequence
    #param seq: sequence of interest
    #param sign: strand sign (+) or (-)
    '''
    if sign == '+':
        return seq
    elif sign == '-':
        return reverse_complement(seq)

genes['new_coordinates'] = genes.progress_apply(lambda row: new_coord(row['start'], row['end'], row['strand']), axis = 1)

if __name__ == "__main__":
    #generate your dataframe:
    df = pd.DataFrame()
    seqs = []
    df['gene_id'] = genes['gene_id']
    df['seqname'] = genes['seqname']
    df['start'] = genes['new_coordinates'].apply(lambda x: x[0])
    df['end'] = genes['new_coordinates'].apply(lambda x: x[1])
    df['strand'] = genes['strand']

    loci = df['seqname'].unique()
    for locus in tqdm(loci):
        print 'NOW PROCESSING DATA FOR CHROMOSOME: {}'.format(locus)
        dfx = df[df['seqname'] == locus]
        handle = open('genomes/{}.fa'.format(locus))
        for record in SeqIO.parse(handle, 'fasta'):
            locus_seq = str(record.seq).upper()
            dfx['sequence'] = dfx.progress_apply(lambda row: locus_seq[row['start']: row['end']], axis = 1)
            dfx['sequence'] = dfx.progress_apply(lambda row: sequencer(row['sequence'], row['strand']), axis = 1)
        seqs += dfx['sequence'].tolist()
    df['sequence'] = seqs
    df.to_pickle('ucsc_gene_regions.pkl')
    print 'DONE!'
