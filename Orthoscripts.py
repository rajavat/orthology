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
    
    # Make ortholog dictionaries
    A_orthdict = dict(zip(orthologies[:, 1], orthologies[:, 0]))
    B_orthdict = dict(zip(orthologies[:, 2], orthologies[:, 0]))

    # Replace genelist values with ortholog dictionaries
    A_data = genelistA.replace({'Name': A_orthdict})
    B_data = genelistB.replace({'Name' : B_orthdict})
    
    # Add column for orthologs: 1 if ortholog, 0 if not
    B_data['Ortholog'] = B_data['Name'].apply(lambda x:1 if 'ortholog' in x.lower() else 0)
    A_data['Ortholog'] = A_data['Name'].apply(lambda x:1 if 'ortholog' in x.lower() else 0)
    
    # Isolate orthologies
    A_ortho = A_data.loc[A_data['Ortholog'] == 1]
    A_dict = dict(zip(A_ortho.Name, A_ortho.Chromosome))

    B_ortho = B_data.loc[B_data['Ortholog'] == 1]
    B_dict = dict(zip(B_ortho.Name, B_ortho.Chromosome))
    
    # Seperate all orthology entries into new dataframe
    AB_data = pd.DataFrame({'Orthologs' : orthologies[:, 0],
                            'speciesA' : orthologies[:, 0],
                            'speciesB' : orthologies[:, 0]})
    
    # Replace location in A and B with orthology dictionary keys
    AB_data['speciesA'] = AB_data['speciesB'].map(A_dict)
    AB_data['speciesB'] = AB_data['speciesB'].map(B_dict)
    
    # Calculate number of orthologs for each pair of chromosomes
    AB_data = AB_data.groupby(['speciesA', 'speciesB']).count().reset_index()
    
    A = A_data.Name.values.tolist()
    B = B_data.Name.values.tolist()
    M = len(list(set(A) & set(B)))
    
    # Define inner function for hypergeometric testing
    def hypertest(chrA, chrB):
        nA = AB_data.loc[(AB_data['speciesA'] == chrA), 'Orthologs'].sum()
        nB = AB_data.loc[(AB_data['speciesB'] == chrB), 'Orthologs'].sum()
        x = AB_data.loc[(AB_data['speciesA'] == chrA) & 
                        (AB_data['speciesB'] == chrB), 'Orthologs'].sum()
    
        p = stats.hypergeom.sf(x - 1, M, nA, nB)
        
        return p

    # Conduct hypergeometric testing
    AB_data['p-Values'] = AB_data.apply(lambda x : hypertest(x['speciesA'], x['speciesB']), axis = 1)
    
    # Apply BH testing correction
    AB_data['Results'], AB_data['p-Values'] = pg.multicomp(AB_data['p-Values'], method = 'fdr_bh')
    
    # Remove all rows that have been rejected in BH correction
    AB_data = AB_data.loc[AB_data['Results'] == True]
    
    return AB_data

def orthoplot(data, titleA, titleB):
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [8, 8]
    sns.set_style("whitegrid")
    sns.scatterplot(data = data, x = 'speciesA', y = 'speciesB', 
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
    
def orthoChr(data):
    """
    inputs:
    dataset in format: speciesA chromosomes | speciesB chromosomes | Orthologs
    
    outputs: 
    """
    matrix = data.pivot(index = 'speciesB', columns='speciesA', values = 'Orthologs')
    matrix = matrix.where(matrix.isnull(), 1).fillna(0).astype(int)
    
    fusions = matrix.sum(axis = 1).loc[lambda x : x >= 2] - 1
    fissions = matrix.sum().loc[lambda x : x >= 2] - 1
    
    print("Fusions:", fusions.sum())
    print(fusions.index.tolist())
    print("Fissions:", fissions.sum())
    print(fissions.index.tolist())