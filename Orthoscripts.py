import pandas as pd
import numpy as np
import io

import scipy.stats as stats
import pingouin as pg

import random
from random import randrange

import matplotlib.pyplot as plt
import seaborn as sns

import argparse

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
def rearrangements(data, outfile = 'F'):
    # Converts table into dotplot
    fissions = data.pivot(index = 'A', columns='B', values = 'Orthologs')
    fusions = data.pivot(index = 'B', columns = 'A', values = 'Orthologs')

    # Picks out all rows and columns with more than one dot
    fissions = fissions.loc[(fissions.where(fissions.isnull(), 1).sum(axis = 1) > 1) | (fissions.sum(axis = 0) > 1)]
    fissions = fissions.stack(dropna = True).reset_index()

    fusions = fusions.loc[(fusions.where(fusions.isnull(), 1).sum(axis = 1) > 1) | (fusions.sum(axis = 0) > 1)]
    fusions = fusions.stack(dropna = True).reset_index()

    # Identify all translocations
    translocations = (fusions[fusions.B.isin(fissions.B)])

    # Remove all translocations from list of fissions
    fissions = (fissions[~ fissions.A.isin(translocations.A)]).groupby('A')['B'].apply(list).reset_index(name = 'B')
    fissions['B'] = [', '.join(map(str, l)) for l in fissions['B']] # Convert list to str

    fusions = (fusions[~ fusions.B.isin(translocations.B)]).groupby('B')['A'].apply(list).reset_index(name = 'A')
    fusions['A'] = [', '.join(map(str, l)) for l in fusions['A']]  # Convert list to str

    # Identify and isolate the translocations
    translocations = translocations.groupby('B')['A'].apply(list).reset_index()
    translocations['A'] = [', '.join(map(str, l)) for l in translocations['A']]
    translocations = translocations.groupby('A')['B'].apply(list).reset_index()
    translocations['B'] = [', '.join(map(str, l)) for l in translocations['B']]

    events = []
    for index, row in fissions.iterrows():
        logstr = (''.join(('Fission of ancestral chromosome ', str(row['A']), ' into ', str(row['B']))))
        print(logstr)
        events.append(logstr)

    for index, row in fusions.iterrows():
        logstr = (''.join(('Fusion of ancestral chromosome ', str(row['A']), ' into ', str(row['B']))))
        print(logstr)
        events.append(logstr)
        
    for index, row in translocations.iterrows():
        logstr = (''.join(('Translocation between ancestral chromosome ', str(row['A']), ' into ', str(row['B']))))
        print(logstr)
        events.append(logstr)
    
    if outfile != 'F':
        with open(outfile, 'w') as out:
            out.write('\n'.join(events))

    fusions['Events'] = 'FUS'
    fissions['Events'] = 'FIS'
    translocations['Events'] = 'TRA'

    log = pd.concat([fusions, fissions, translocations])
    log = log[['Events', 'A', 'B']]
    
    return log

def simulator(Nchr = 20, Ngene = 100, Nevents = 10):

    # Make ancestor genome
    def makeancestor(Nchr, Ngene):
        ancestor = pd.DataFrame(columns = ['Chromosome'])
        for i in range(Nchr):
            row = {'Chromosome' : (i + 1)}
            for i in range(Ngene):
                    ancestor = pd.concat([ancestor, pd.DataFrame([row])], ignore_index = True)
        ancestor['Name'] = (ancestor.reset_index().index + 1)
        # ancestor['Name'] = 'g_' + ancestor['Name'].astype(str)

        return ancestor

    # Dummy BED files :: type 'anc' for ancestor, 'des' for descendant
    def dummyBED(genome, type):
        if type == 'anc':
            genome['Chromosome'] = 'AncChr' + genome['Chromosome'].astype(str)
            genome['Name'] = 'ancg_' + genome['Name'].astype(str)
            
        if type == 'des':
            genome['Chromosome'] = 'Chr' + genome['Chromosome'].astype(str)
            genome['Name'] = 'g_' + genome['Name'].astype(str)
        
        genome['Start'] = np.arange(len(genome))
        genome['End'] = np.arange(len(genome)) + 5
        
        genome = genome[['Chromosome', 'Start', 'End', 'Name']]
            
        return genome

    # Dummy ortholog file
    def dummyOrthologs(genome):
        orthologs = pd.DataFrame()
        
        orthologs['Orthologs'] = np.arange(len(genome)) + 1
        orthologs['speciesA'] = np.arange(len(genome)) + 1
        orthologs['speciesB'] = np.arange(len(genome)) + 1
        
        orthologs['Orthologs'] = 'orthologs_' + orthologs['Orthologs'].astype(str)
        orthologs['speciesA'] = 'ancg_' + orthologs['speciesA'].astype(str)
        orthologs['speciesB'] = 'g_' + orthologs['speciesB'].astype(str)
        
        orthologs = orthologs.to_numpy()
        
        return orthologs

    def fusion(genome, chr, mixing = 0):
        '''
        inputs: 
        ancestor : 
            df with chromosome name | gene name
        mixing : 
            float between 0 and 1, where 1 implies 
            extreme mixing and 0 implies no mixing
        '''
        
        # Randomly select two chromosomes to fuse
        A = random.choice(chr)
        B = random.choice(chr)
        
        chr = [x for x in chr if x not in (A, B)]
        
        if A == B: # Just so the same chromosome isn't selected twice
            B = random.choice(chr)
            if A == B:
                B = random.choice(chr)

        fusion = ancestor.loc[ancestor['Chromosome'].isin([A, B])]
        
        # Apply mixing if required
        if mixing > 0:
            genes = fusion['Name'].to_numpy()
            n = len(genes)
            for i in range(int(mixing * n)):
                g1, g2 = randrange(n), randrange(n)
                genes[g2], genes[g1] = genes[g1], genes[g2]

            fusion['Name'] = genes
            chrom = f'{A}x{B}'
            fusion['Chromosome'] = chrom
            
        else:
            chrom = f'{A}+{B}'
            fusion['Chromosome'] = chrom
            
        log = ['FUS', 
            ('AncChr' + str(A) + ', ' + 'AncChr' + str(B)), 
            ('Chr' + str(chrom))]
        
        print(f'Fusion of AncChr{A} and AncChr{B} into {chrom}')
        
        # Remove the unfused chromosomes
        genome.drop(genome[genome['Chromosome'].isin([A, B])].index, inplace = True)
        genome = pd.concat([genome, fusion])
        
        return genome, log, chr

    def fission(genome, chr):
        # Randomly select a chromosome for fission
        A = random.choice(chr)
        fission = genome.loc[genome['Chromosome'] == A]
        chr.remove(A)
        
        pos = random.choice(range(1, Ngene))

        # Add the new chromosomes back into the genome
        chr1 = fission.iloc[: pos]
        chr1['Chromosome'] = f'{A}_1'
        
        chr2 = fission.iloc[pos :]
        chr2['Chromosome'] = f'{A}_2'
        
        # Remove the fission chromosome from the genome
        genome = pd.concat([genome, chr1, chr2])
        genome = genome[genome.Chromosome != A]
        log = ['FIS', 
            ('AncChr' + str(A)), 
            ('Chr' + str(A) + '_1' + ', Chr' + str(A) + '_2')]
        
        print(f'Fission of AncChr{A} into Chr{A}_1, Chr{A}_2')
        
        return genome, log, chr

    def translocation(genome, chr):
        # Randomly select two chromosomes for translocation
        A = random.choice(chr)
        B = random.choice(chr)
        
        if A == B: # Just so the same chromosome isn't selected twice
            B = random.choice(chr)
        
        chr = [x for x in chr if x not in (A, B)]
        
        chrA = genome.loc[genome['Chromosome'] == A]
        chrB = genome.loc[genome['Chromosome'] == B]
        
        # Randomly select two break point positions
        posA = random.choice(range(1, Ngene))
        posB = random.choice(range(1, Ngene))
        
        r = np.random.uniform()
        if r <= 0.30:
            # Join the fragments to form recombinant chromosomes
            chr1 = pd.concat([chrA.iloc[: posA], chrB.iloc[posB :]])
            chr1['Chromosome'] = f'{A};{B}'
            chr2 = pd.concat([chrB.iloc[: posB], chrA.iloc[posA :]])
            chr2['Chromosome'] = f'{B};{A}'
            log = ['TR', 
                ('AncChr' + str(A) + ', AncChr' + str(B)), 
                (('Chr' + str(A) + ';' + str(B) + ', Chr' + str(B) + ';' + str(A)))]
            
            print(f'Translocation between AncChr{A} and AncChr{B} into Chr{A};{B} and Chr{B};{A}')
        
        else:
            chr1 = pd.concat([chrA.iloc[: posA]])
            chr1['Chromosom'] = f'{A};'
            chr2 = pd.concat([chrB, chrA.iloc[posA :]])
            chr2['Chromosome'] = f'{B};{A}'
        
            log = ['TR', 
                ('AncChr' + str(A) + ', AncChr' + str(B)),
                ('Chr' + str(B) + ';' + str(A))]
            
            print(f'Translocation between AncChr{A} and AncChr{B} into Chr{B};{A}')

        # Remove the original chromosomes from the genome
        genome = pd.concat([genome, chr1, chr2]).drop(genome[(genome['Chromosome'] == A) & (genome['Chromosome'] == B)].index)
        return genome, log, chr


    # Apply macro-rearrangements to the ancestor
    ancestor = makeancestor(Nchr, Ngene)
    chr = ancestor.Chromosome.unique().tolist()
    speciesA = ancestor.copy()

    events = []
    for event in range(Nevents):
        r = np.random.uniform()
        
        if r <= 0.30:
            if len(ancestor) < 2: continue
            speciesA, log, chr = fission(speciesA, chr)
            events.append(log)
        
        elif r <= 0.45:
            speciesA, log, chr = translocation(speciesA, chr)
            events.append(log)
        
        elif r <= 0.70:
            speciesA, log, chr = fusion(speciesA, chr)
            events.append(log)
        
        elif r <= 0.99:
            speciesA, log, chr = fusion(speciesA, chr, mixing = 0.5)
            events.append(log)
            
        else:
            # speciesA, log, chr = syntenyloss(speciesA, chr)
            # events.append(log)
            # print(log)
            continue
        
    ancestor = dummyBED(ancestor, 'anc')
    speciesA = dummyBED(speciesA, 'des')
    orthologs = dummyOrthologs(ancestor)
    
    events =  pd.DataFrame(events, columns = ['Events', 'A', 'B'])
    
    return ancestor, speciesA, orthologs, events