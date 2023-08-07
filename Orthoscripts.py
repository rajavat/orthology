import pandas as pd
import numpy as np
import scipy.stats as stats
import pingouin as pg
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
from pathlib import PosixPath

# reads a BED file
def readBED(file, str = 't'):
    cols = ['Chromosome', 'Start', 'End', 'Name']
    if str == 't':
        delim = '\t'
    if str == 's':
        delim = '\s+'
    df = pd.read_csv(file, sep = delim, header = None, 
                     names = cols, usecols = [0,1,2,3])
    return df

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

# returns a df with the number of orthologs for each pair of chromosomes
def orthologies(genelistA, genelistB, orthologies):
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
    
    return AB_data 

# returns a df with significant chromosome pairs a
def sigorthologies(genelistA, genelistB, orthologies):
    
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

# plots
def orthoplot(data, titleA, titleB, x = 'A', y = 'B'):

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

# counts the rearrangements
def rearrangements(data, outfile):
    # Converts table into dotplot
    fissions = data.pivot(index = 'A', columns='B', values = 'Orthologs')
    
    # Picks out all rows and columns with more than one dot
    fissions = fissions.loc[(fissions.where(fissions.isnull(), 1).sum(axis=1) > 1) | (fissions.sum(axis=0) > 1)]
    fissions = fissions.stack(dropna = True).reset_index().groupby('A')['B'].apply(list).reset_index(name = 'B')
    fissions['B'] = [', '.join(map(str, l)) for l in fissions['B']] # Convert list to str

    # Identify all translocations
    translocations = fissions.groupby('B').filter(lambda g: len(g) > 1)
    # Remove all translocations from list of fissions
    fissions = fissions[~ fissions.isin(translocations)]
    fissions.dropna(inplace = True)

    # Converts table into dotplot
    fusions = data.pivot(index = 'B', columns = 'A', values = 'Orthologs')
    fusions = fusions.loc[(fusions.where(fusions.isnull(), 1).sum(axis=1) > 1) | (fusions.sum(axis=0) > 1)]
    fusions = fusions.stack(dropna = True).reset_index()
    
    # Picks out all rows and columns with more than one dot
    translocations = fusions.groupby('A').filter(lambda g: len(g) > 1)
    
    # Remove all translocations from list of fissions
    fusions = fusions[~ fusions.isin(translocations)]

    # Identify and isolate the translocations
    translocations = translocations.groupby('B')['A'].apply(list).reset_index()
    translocations['A'] = [', '.join(map(str, l)) for l in translocations['A']]
    translocations = translocations.groupby('A')['B'].apply(list).reset_index()
    translocations['B'] = [', '.join(map(str, l)) for l in translocations['B']]

    fusions = fusions.groupby('B')['A'].apply(list).reset_index(name = 'A')
    fusions['A'] = [', '.join(map(str, l)) for l in fusions['A']]

    events = []
    for index, row in fissions.iterrows():
        events.append(''.join(('Fission of ancestral chromosome', row['A'], 'into', row['B'])))

    for index, row in fusions.iterrows():
        events.append(''.join(('Fusion of ancestral chromosome', row['A'], 'into', row['B'])))
        
    for index, row in translocations.iterrows():
        events.append(''.join(('Translocation of ancestral chromosome', row['A'], 'into', row['B'])))
        
    with open(outfile, 'w') as out:
        out.write('\n'.join(events))
    
    return events