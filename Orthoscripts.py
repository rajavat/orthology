import pandas as pd
import numpy as np
import scipy.stats as stats
import pingouin as pg
import matplotlib.pyplot as plt
import seaborn as sns

def readBED(file):
    rawdata = []
    with open(file) as f:
        for line in f:
            rawdata.append(line.strip().split())
    data = pd.DataFrame(rawdata,
                      columns = ['Chromosome', 'Start', 'End', 'Name', 'Dot'])
    
    return data

def orthofy(genelistA, genelistB, orthologies):
    
    """
    inputs:
    genelistA: gene list for species A
    genelistB: gene list for species B
    orthologies: orthology dataset
    
    outputs: dataframe with significant ortholog combinations 
             and their location in species A and B and p-Values
    """
    
    # Make ortholog dictionaries (ortholog : gene name)
    orthdictA = dict(zip(orthologies[:, 1], orthologies[:, 0]))
    orthdictB = dict(zip(orthologies[:, 2], orthologies[:, 0]))

    # Replace genelist values with ortholog dictionary keys
    A_data = genelistA.copy()
    B_data = genelistB.copy()
    A_data['Name'] = A_data['Name'].map(lambda x: orthdictA.get(x, x))
    B_data['Name'] = B_data['Name'].map(lambda x: orthdictB.get(x, x))
    
    # Make orthology location dictionaries (ortholog: chromosome)
    dictA = dict(zip(A_data.loc[A_data['Name'].str.contains('ortholog')].Name, 
                     A_data.loc[A_data['Name'].str.contains('ortholog')].Chromosome))
    dictB = dict(zip(B_data.loc[B_data['Name'].str.contains('ortholog')].Name, 
                     B_data.loc[B_data['Name'].str.contains('ortholog')].Chromosome))
    
    # Seperate all orthology entries into new dataframe
    AB_data = pd.DataFrame({'Orthologs': orthologies[:, 0],
                            'A' : orthologies[:, 0],
                            'B' : orthologies[:, 0]})
    
    # Replace location in A and B with ortholog location dictionary keys
    AB_data['A'] = AB_data['A'].map(dictA)
    AB_data['B'] = AB_data['B'].map(dictB)
    
    # Calculate number of orthologs for each pair of chromosomes
    AB_data = AB_data.groupby(['A', 'B']).count().reset_index()
    
    M = len(list(set(A_data.Name.values.tolist()) & set(B_data.Name.values.tolist())))
    
    # Define inner function for hypergeometric testing
    def hypertest(chrA, chrB):
        nA = AB_data.loc[(AB_data['A'] == chrA), 'Orthologs'].sum()
        nB = AB_data.loc[(AB_data['B'] == chrB), 'Orthologs'].sum()
        x = AB_data.loc[(AB_data['A'] == chrA) & 
                        (AB_data['B'] == chrB), 'Orthologs'].sum()
    
        p = stats.hypergeom.sf(x - 1, M, nA, nB)
        
        return p

    # Conduct hypergeometric testing
    AB_data['p-Values'] = AB_data.apply(lambda x : hypertest(x['A'], x['B']), axis = 1)
    
    # Apply BH testing correction
    AB_data['Results'], AB_data['p-Values'] = pg.multicomp(AB_data['p-Values'], method = 'fdr_bh')
    
    # Remove all rows that have been rejected in BH correction
    AB_data = AB_data.loc[AB_data['Results'] == True]
    
    return AB_data

def orthoplot(data, titleA, titleB, x, y):
    """
    input: 
    dataset
    titleA: x-axis title
    titleB: y-axis title
    x: species on x-axis
    y: species on y-axis
    """
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [8, 8]
    sns.set_style("whitegrid")
    sns.scatterplot(data = data, x = x, y = y, 
                size = 'Orthologs', sizes = (10, 200),
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
    """
    inputs:
    dataset in format: speciesA chromosomes | speciesB chromosomes | Orthologs
    
    outputs: 
    """
    matrix = data.pivot(index = 'B', columns='A', values = 'Orthologs')
    matrix = matrix.where(matrix.isnull(), 1).fillna(0).astype(int)
    
    fusions = matrix.sum(axis = 1).loc[lambda x : x >= 2] - 1
    fissions = matrix.sum().loc[lambda x : x >= 2] - 1
    
    print("Fusions:", fusions.sum())
    print(fusions.index.tolist())
    print("Fissions:", fissions.sum())
    print(fissions.index.tolist())
    
def orthFix(orthology, col, string, position):
    """
    inputs:
    orthology: orthology input
    col: column with suffix: A or B
    string: string to remove
    position: 0 for suffix, 1 for prefix
    """
    orthology = pd.DataFrame(orthology, columns = ['Code', 'A', 'B'])
    orthology[col] = orthology[col].str.rsplit(string).str.get(position)
    orthology = orthology.to_numpy()
    
    return orthology
    
def unscaff(data, scope):
    """
    inputs
    data: dataframes
    scope: level at which to filter scaffolds
    """
    scaffs = data.groupby('Chromosome').size()
    scaffs = scaffs.reset_index()

    scaffs.columns = ['Chromosome', 'Count']
    scaffs = scaffs.loc[scaffs['Count'] >= scope]

    scaffolds = scaffs.Chromosome.tolist() # Remove all values from non-chromosome scaffolds
    data = data.loc[data['Chromosome'].isin(scaffolds)]
    
    return data