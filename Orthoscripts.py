import pandas as pd
import numpy as np
import scipy.stats as stats
import pingouin as pg
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
from pathlib import PosixPath

# from GFF2BED package: https://pypi.org/project/gff2bed/#files --------------------

# parse an entry from the "attr" column of a GFF3 file and return it as a dict
def parse_gff_attributes(attr: str, graceful: bool = False):
    """
    inputs:
    attr : str
        feature attribute string
    graceful : bool
        if False, throw an error when a malformed tag-value pair is encountered
        if True, ignore malformed pairs gracefully

    returns:
    dict
        attr entries as a dict
    """

    return dict(pair.split('=') for pair in attr.split(';')
                if ('=' in pair) or not graceful)

# parse a GFF3 file and yield its lines as tuples
def readGFF(gff_file, type: str = 'gene', parse_attr: bool = True,
          graceful: bool = False):
    """
    inputs:
    gff_file
        String, PosixPath, or file-like object representing a GFF3 file
    type
        string indicating feature type to include, or None to include all
        features [gene]
    parse_attr : bool
        if False, do not parse attributes [True]
    graceful : bool
        if False, throw an error when a malformed tag-value pair is encountered
        if True, ignore malformed pairs gracefully [False]
    
    returns:
    seqid, start, end, strand, attr
        coordinates of a feature
    """

    with (
        gff_file if not isinstance(gff_file, (str, PosixPath))
        else gzip.open(gff_file, 'rt') if str(gff_file).endswith('.gz')
        else open(gff_file, 'r')
    ) as f:
        for line in f:
             if not line.startswith('#'):
                seqid, _, t, start, end, _, strand, _, attr = line.rstrip().split('\t')
                if ((t == type) or (type is None)):
                    if parse_attr:
                        yield (seqid, int(start), int(end), strand,
                            parse_gff_attributes(attr, graceful=graceful))
                    else:
                        yield seqid, int(start), int(end), strand, '.'

# converts rows of GFF file into BED data                
def convert(gff_data, tag: str = 'protein_id'):
    """
    input:
    gff_data
        iterable of data from gff2bed.parse()
    tag : str
        GFF3 attribute tag to parse [ID]

    returns:
    tuple
        a row of BED data
    """

    for seqid, start, end, strand, attr in gff_data:
        yield seqid, start - 1, end, attr[tag], 0, strand
        
# -----------------------------------------------------------------------------------

# reads a BED file into a dataframe - make it non-specific to row number
def readBED(file):
    rawdata = []
    with open(file) as f:
        for line in f:
            rawdata.append(line.strip().split())
    data = pd.DataFrame(rawdata,
                      columns = ['Chromosome', 'Start', 'End', 'Name', 'Dot'])
    
    return data

# removes suffix/prefix from ortholog data
def orthFix(orthology, col, string, position):
    
    """
    inputs:
    orthology: 
        orthology input
    col: 
        column with suffix: A or B
    string: 
        string to remove
    position: 
        0 for suffix, 1 for prefix
    
    returns:
        orthology input without the suffix/prefix
    """
    
    orthology = pd.DataFrame(orthology, columns = ['Code', 'A', 'B'])
    orthology[col] = orthology[col].str.rsplit(string).str.get(position)
    orthology = orthology.to_numpy()
    
    return orthology

# selects and removes all genelist entries from non-chromosome scaffolds
def unscaff(data, scope = 100):
    
    """
    inputs
    data: 
        dataframe
    scope: 
        level at which to filter scaffolds, default is 100
    """
    
    scaffs = data.groupby('Chromosome').size()
    scaffs = scaffs.reset_index()

    scaffs.columns = ['Chromosome', 'Count']
    scaffs = scaffs.loc[scaffs['Count'] >= scope]

    scaffolds = scaffs.Chromosome.tolist() # Remove all values from non-chromosome scaffolds
    data = data.loc[data['Chromosome'].isin(scaffolds)]
    
    return data

# converts genelist, orthology file into a df with significant ortholog combinations
def orthofind(genelistA, genelistB, orthologies):
    
    """
    inputs:
    genelistA: 
        gene list for species A
    genelistB: 
        gene list for species B
    orthologies: 
        orthology dataset 
    
    outputs: dataframe with significant ortholog combinations 
             and their location in species A and B and p-Values
    """
    
    # make ortholog dictionaries (ortholog : gene name)
    orthdictA = dict(zip(orthologies[:, 1], orthologies[:, 0]))
    orthdictB = dict(zip(orthologies[:, 2], orthologies[:, 0]))

    # replace genelist values with ortholog dictionary keys
    A_data = genelistA.copy()
    B_data = genelistB.copy()
    A_data['Name'] = A_data['Name'].map(lambda x: orthdictA.get(x, x))
    B_data['Name'] = B_data['Name'].map(lambda x: orthdictB.get(x, x))
    
    # make orthology location dictionaries (ortholog : chromosome)
    dictA = dict(zip(A_data.loc[A_data['Name'].str.contains('ortholog')].Name, 
                     A_data.loc[A_data['Name'].str.contains('ortholog')].Chromosome))
    dictB = dict(zip(B_data.loc[B_data['Name'].str.contains('ortholog')].Name, 
                     B_data.loc[B_data['Name'].str.contains('ortholog')].Chromosome))
    
    # seperate all orthology entries into new dataframe
    AB_data = pd.DataFrame({'Orthologs': orthologies[:, 0],
                            'A' : orthologies[:, 0],
                            'B' : orthologies[:, 0]})
    
    # replace location in A and B with ortholog location dictionary keys
    AB_data['A'] = AB_data['A'].map(dictA)
    AB_data['B'] = AB_data['B'].map(dictB)
    
    # calculate number of orthologs for each pair of chromosomes
    AB_data = AB_data.groupby(['A', 'B']).count().reset_index()
    
    M = len(list(set(A_data.Name.values.tolist()) & set(B_data.Name.values.tolist())))
    
    # define inner function for hypergeometric testing
    def hypertest(chrA, chrB):
        nA = AB_data.loc[(AB_data['A'] == chrA), 'Orthologs'].sum()
        nB = AB_data.loc[(AB_data['B'] == chrB), 'Orthologs'].sum()
        x = AB_data.loc[(AB_data['A'] == chrA) & 
                        (AB_data['B'] == chrB), 'Orthologs'].sum()
    
        p = stats.hypergeom.sf(x - 1, M, nA, nB)
        
        return p

    # conduct hypergeometric testing
    AB_data['p-Values'] = AB_data.apply(lambda x : hypertest(x['A'], x['B']), axis = 1)
    
    # apply BH testing correction
    AB_data['Results'], AB_data['p-Values'] = pg.multicomp(AB_data['p-Values'], method = 'fdr_bh')
    
    # remove all rows that have been rejected in BH correction
    AB_data = AB_data.loc[AB_data['Results'] == True]
    
    return AB_data

# plots an oxford dot plot 
def orthoplot(data, titleA, titleB, x, y):

    """
    input: 
    dataset:
        species A chromosome | species B chromosome | n. orthologs
    titleA: 
        x-axis title
    titleB: 
        y-axis title
    x: 
        species on x-axis
    y: 
        species on y-axis
    """

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [8, 8]
    sns.set_style("whitegrid")
    sns.scatterplot(data = data, x = x, y = y, 
                size = 'Orthologs', sizes = (50, 200),
                hue = 'Orthologs', palette = "crest")

    plt.xlabel(titleA, fontsize = 13, 
           labelpad = 15, style = 'italic')
    plt.ylabel(titleB, fontsize = 13, 
           labelpad = 15, style = 'italic')

    plt.xticks(rotation='vertical', fontsize = 9)
    plt.yticks(fontsize = 8)

    plt.legend(bbox_to_anchor=(1, 1), loc='upper left', 
           fontsize=10, title = 'Orthologous genes', frameon = False)

    plt.show()
    
def rearrangements(data):
    fusions = data.pivot(index = 'B', columns='A', values = 'Orthologs')
    fusions = fusions.loc[(fusions.where(fusions.isnull(), 1).sum(axis=1) > 1) | (fusions.sum(axis=0) > 1)]
    fusions = fusions.stack(dropna = True).reset_index().groupby('B')['A'].apply(list).reset_index(name = 'A')

    fissions = data.pivot(index = 'A', columns='B', values = 'Orthologs')
    fissions = fissions.loc[(fissions.where(fissions.isnull(), 1).sum(axis=1) > 1) | (fissions.sum(axis=0) > 1)]
    fissions = fissions.stack(dropna = True).reset_index().groupby('A')['B'].apply(list).reset_index(name = 'B')
    
    f = open("rearrangements.txt", "w+")
    for index, row in fissions.iterrows():
        print('Fission of ancestral chromosome', row['A'], 'into', row['B'])
        f.write('{0} {1} {2} {3}\n'.format('Fission of', row['B'], 'into', row['A']))
    for index, row in fusions.iterrows():
        print('Fusion of ancestral chromosomes', row['A'], 'into', row['B'])
        f.write('{0} {1} {2} {3}\n'.format('Fusion of', row['A'], 'into', row['B']))
    f.close()
    