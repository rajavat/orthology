import pandas as pd
import numpy as np
import Orthoscripts

# Import genelists
Braflo = Orthoscripts.readBED("Data/Genelists/Branchiostoma.floridae.genelist.bed")
Pecmax = Orthoscripts.readBED("Data/Genelists/Pecten.maximus.genelist.bed")
Holleu = Orthoscripts.readBED("Data/Genelists/Holothuria.leucospilota.genelist.bed")

# Import orthologies
Pecmax_Braflo = np.loadtxt("Orthology pipeline/orthologs/Pecmax+Braflo_sensitive.txt", dtype = "str")
Pecmax_Holleu = np.loadtxt("Orthology pipeline/orthologs/Pecmax+Holleu_sensitive.txt", dtype = "str")
Holleu_Braflo = np.loadtxt("Orthology pipeline/orthologs/Holleu+Braflo_sensitive.txt", dtype = "str")

# Modified version of ortholog function - outputs just the orthologies in df
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
    genelistA['Name'] = genelistA['Name'].map(lambda x: orthdictA.get(x, x))
    genelistB['Name'] = genelistB['Name'].map(lambda x: orthdictB.get(x, x))
    
    # Make orthology location dictionaries (ortholog : chromosome)
    dictA = dict(zip(genelistA.loc[genelistA['Name'].str.contains('ortholog')].Name, 
                     genelistA.loc[genelistA['Name'].str.contains('ortholog')].Chromosome))
    dictB = dict(zip(genelistB.loc[genelistB['Name'].str.contains('ortholog')].Name, 
                     genelistB.loc[genelistB['Name'].str.contains('ortholog')].Chromosome))
    
    # Seperate all orthology entries into new dataframe
    AB_data = pd.DataFrame({'Orthologs': orthologies[:, 0],
                            'A' : orthologies[:, 0],
                            'B' : orthologies[:, 0]})
    
    # Replace location in A and B with ortholog location dictionary keys
    AB_data['A'] = AB_data['A'].map(dictA)
    AB_data['B'] = AB_data['B'].map(dictB)
    
    return AB_data

PB = orthofy(Pecmax, Braflo, Pecmax_Braflo)
PB = PB.dropna()

# Make matrix with corresponding chromosomes
Amp = ['BFL_11', 'BFL_10', 'BFL_16', 'BFL_8', 'BFL_3', 'BFL_1', 'BFL_18', 'BFL_14', 'BFL_15', 'BFL_5', 'BFL_7', 'BFL_3', 'BFL_17', 'BFL_3', 'BFL_19', 'BFL_12', 'BFL_1', 'BFL_13', 'BFL_2', 'BFL_2', 'BFL_6', 'BFL_9', 'BFL_4', 'BFL_4']
Sca = ['PYE_10', 'PYE_13', 'PYE_1', 'PYE_1', 'PYE_17', 'PYE_5', 'PYE_19', 'PYE_15', 'PYE_4', 'PYE_6', 'PYE_7', 'PYE_2', 'PYE_18', 'PYE_2', 'PYE_3', 'PYE_14', 'PYE_16', 'PYE_2', 'PYE_4', 'PYE_9', 'PYE_8', 'PYE_3', 'PYE_11', 'PYE_12']
Anc = ['G', 'B1', 'B2', 'M', 'C2', 'A1aA1b', 'B3', 'P', 'L', 'EaEb', 'F', 'QbQa', 'J1', 'QcQd', 'O2', 'N', 'A2', 'H', 'J2', 'C1', 'D', 'K', 'I', 'O1']
ChrCorr = np.column_stack((Sca, Amp, Anc))

# Make dataframe with corresponding chromosomes
PBgenes = pd.DataFrame()
for i in range (0, 24): 
    PBorthologs = PB.loc[(PB['A'] == ChrCorr[i, 0]) & (PB['B'] == ChrCorr[i, 1])]
    PBorthologs['Chr'] = ChrCorr[i, 2]

    PBgenes = pd.concat([PBgenes, PBorthologs])

# Manually add PYE_12
PBorthologs = PB.loc[(PB['A'] == 'PYE_12') & (PB['B'] != 'BFL_4')]
PBorthologs['Chr'] = 'R'
PBgenes = pd.concat([PBgenes, PBorthologs])

PBgenes['BGenes'] = PBgenes.loc[:, 'Orthologs']
PBgenes = PBgenes.rename(columns = {'Orthologs' : 'PGenes'})
PBgenes = PBgenes[['Chr', 'A', 'PGenes', 'B', 'BGenes']]

# Make reverse ortholog dictionaries (ortholog : gene name)
orthdictA = dict(zip(Pecmax_Braflo[:, 0], Pecmax_Braflo[:, 1]))
orthdictB = dict(zip(Pecmax_Braflo[:, 0], Pecmax_Braflo[:, 2]))

# Replace values
PBgenes['PGenes'] = PBgenes['PGenes'].map(lambda x: orthdictA.get(x, x))
PBgenes['BGenes'] = PBgenes['BGenes'].map(lambda x: orthdictB.get(x, x))

# Make dictionaries (H gene name : P/B gene name)
orthdictP = dict(zip(Pecmax_Holleu[:, 1], Pecmax_Holleu[:, 2]))
orthdictB = dict(zip(Holleu_Braflo[:, 2], Holleu_Braflo[:, 1]))

# Replace values
PBgenes['PGenes'] = PBgenes['PGenes'].map(lambda x: orthdictP.get(x, x))
PBgenes['BGenes'] = PBgenes['BGenes'].map(lambda x: orthdictB.get(x, x))

# Select all values orthologous in both columns
Ancestor = PBgenes.loc[(PBgenes['PGenes'].str.contains('gene-HOLleu_')) & 
                       (PBgenes['BGenes'].str.contains('gene-HOLleu_'))]

Ancestor = Ancestor.rename(columns = {'Chr' : 'Chromosome',
                                      'PGenes' : 'Name', 
                                      'A' : 'Pchr',
                                      'B' : 'Bchr'})
Ancestor = Ancestor[['Chromosome', 'Name', 'Pchr', 'Bchr']]
Ancestor['Hchr'] = Ancestor.loc[:, 'Name']

# Add column with native sea cucumber chromosome
Hol = Holleu.to_numpy()
orthdictHchr = dict(zip(Hol[:, 3], Hol[:, 0]))
Ancestor['Hchr'] = Ancestor['Hchr'].map(lambda x: orthdictHchr.get(x, x))

# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  
#   print(Ancestor)
    
# Make ancestor BED file
Ancestorgenelist = pd.DataFrame()
Ancestorgenelist['Chromosome'] = Ancestor.loc[:, 'Chromosome']
Ancestorgenelist['Start'] = Holleu.loc[:, 'Start']
Ancestorgenelist['End'] = Holleu.loc[:, 'End']
Ancestorgenelist['Name'] = Ancestor.loc[:, 'Name']
Ancestorgenelist['Dot'] = Holleu.loc[:, 'Dot']

np.savetxt(r'Data/Genelists/Ancestor.genelist.bed', Ancestorgenelist.values, fmt = '%s')