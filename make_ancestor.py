import pandas as pd
import numpy as np
import syntiny

# Disable chained assignments
pd.options.mode.chained_assignment = None

# Import genelists
Braflo = syntiny.readBED("Data/Genelists/Branchiostoma.floridae.genelist.bed")
Pecmax = syntiny.readBED("Data/Genelists/Pecmax.genelist.bed")
Holleu = syntiny.readBED("Data/Genelists/Holothuria.leucospilota.genelist.bed")
Ephmue = syntiny.readBED("Data/Genelists/Ephmue.genelist.bed")

# Import orthologies
Pecmax_Braflo = np.loadtxt("Orthology pipeline/orthologs/Pecmax+Braflo_sensitive.txt", dtype = "str")
Pecmax_Holleu = np.loadtxt("Orthology pipeline/orthologs/Pecmax+Holleu_sensitive.txt", dtype = "str")
Holleu_Braflo = np.loadtxt("Orthology pipeline/orthologs/Holleu+Braflo_sensitive.txt", dtype = "str")
Holleu_Ephmue = np.loadtxt("Orthology pipeline/orthologs/Holleu+Ephmue_sensitive.txt", dtype = "str")
Braflo_Ephmue = np.loadtxt("Orthology pipeline/orthologs/Braflo+Ephmue_sensitive.txt", dtype = "str")
Pecmax_Ephmue = np.loadtxt("Orthology pipeline/orthologs/Pecmax+Ephmue_sensitive.txt", dtype = "str")

Pecmax = Pecmax.loc[Pecmax['Chromosome'].str.contains('PYE_')]
Ephmue = Ephmue.loc[Ephmue['Chromosome'].str.contains('EMU_')]
Braflo = Braflo.loc[Braflo['Chromosome'].str.contains('BFL_')]

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

# BILATERIAN ANCESTOR -----------------------------------------------------------------------
PB = orthofy(Pecmax, Braflo, Pecmax_Braflo).dropna()

# Make matrix with corresponding chromosomes
BFL = ['BFL_11', 'BFL_10', 'BFL_16', 'BFL_8', 'BFL_3', 
       'BFL_1', 'BFL_18', 'BFL_14', 'BFL_15', 'BFL_5', 
       'BFL_7', 'BFL_17', 'BFL_19', 'BFL_12', 'BFL_1', 
       'BFL_13', 'BFL_2', 'BFL_2', 'BFL_6', 'BFL_9', 
       'BFL_4', 'BFL_4']

PYE = ['PYE_10', 'PYE_13', 'PYE_1', 'PYE_1', 'PYE_17', 
       'PYE_5', 'PYE_19', 'PYE_15', 'PYE_4', 'PYE_6', 
       'PYE_7', 'PYE_18', 'PYE_3', 'PYE_14', 'PYE_16', 
       'PYE_2', 'PYE_4', 'PYE_9', 'PYE_8', 'PYE_3', 
       'PYE_11', 'PYE_12']

Anc = ['G', 'B1', 'B2', 'M', 'C2', 'A1aA1b', 'B3', 
       'P', 'L', 'EaEb', 'F', 'J1', 'O2', 'N', 'A2', 
       'H', 'J2', 'C1', 'D', 'K', 'I', 'O1']

ChrCorr = np.column_stack((PYE, BFL, Anc))

# Make dataframe with corresponding chromosomes
PBgenes = pd.DataFrame()
for i in range (0, 22): 
    PBorthologs = PB.loc[(PB['A'] == ChrCorr[i, 0]) & (PB['B'] == ChrCorr[i, 1])]
    PBorthologs['Chr'] = ChrCorr[i, 2]

    PBgenes = pd.concat([PBgenes, PBorthologs])

# Manually add R
R = PB.loc[(PB['A'] == 'PYE_12') & (PB['B'] != 'BFL_4')]
R['Chr'] = 'R'
PBgenes = pd.concat([PBgenes, R])

# Manually add QbQa and QdQc
PBorthologs = PB.loc[(PB['A'] == 'PYE_2') & (PB['B'] == 'BFL_3')]
QbQa = PBorthologs.iloc[:101,:]
QbQa['Chr'] = 'QbQa'
QcQd = PBorthologs.iloc[101:,:]
QcQd['Chr'] = 'QcQd'

PBgenes = pd.concat([PBgenes, QbQa, QcQd])

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

BilAnc = Ancestorgenelist

np.savetxt(r'Data/Genelists/BilAnc.genelist.bed', Ancestorgenelist.values, fmt = '%s')

# ANIMAL ANCESTOR ----------------------------------------------------------------------
PE = orthofy(Pecmax, Ephmue, Pecmax_Ephmue).dropna()

# Make matrix with corresponding chromosomes
PYE = ['PYE_5', 'PYE_5', 'PYE_6', 'PYE_6', 'PYE_2', 
       'PYE_2', 'PYE_2', 'PYE_2']
EMU = ['EMU_19', 'EMU_14', 'EMU_01', 'EMU_14', 'EMU_2', 
       'EMU_10', 'EMU_7', 'EMU_10']
Anc = ['A1a', 'A1b', 'Ea', 'Eb', 'Qa', 'Qb', 'Qc', 'Qd']
ChrCorr = np.column_stack((PYE, EMU, Anc))

# Make dataframe with corresponding chromosomes
PEgenes = pd.DataFrame()
for i in range (0, 8): 
    PEorthologs = PE.loc[(PE['A'] == ChrCorr[i, 0]) & (PE['B'] == ChrCorr[i, 1])]
    PEorthologs['Chr'] = ChrCorr[i, 2]

    PEgenes = pd.concat([PEgenes, PEorthologs])
    
PEgenes['EGenes'] = PEgenes.loc[:, 'Orthologs']
PEgenes = PEgenes.rename(columns = {'Orthologs' : 'PGenes'})

# Make reverse ortholog dictionaries (ortholog : gene name)
orthdictA = dict(zip(Pecmax_Ephmue[:, 0], Pecmax_Ephmue[:, 1]))
orthdictB = dict(zip(Pecmax_Ephmue[:, 0], Pecmax_Ephmue[:, 2]))

# Replace values
PEgenes['PGenes'] = PEgenes['PGenes'].map(lambda x: orthdictA.get(x, x))
PEgenes['EGenes'] = PEgenes['EGenes'].map(lambda x: orthdictB.get(x, x))

# Make dictionaries (H gene name : P/B gene name)
orthdictE = dict(zip(Holleu_Ephmue[:, 2], Holleu_Ephmue[:, 1]))
orthdictP = dict(zip(Pecmax_Holleu[:, 1], Pecmax_Holleu[:, 2]))

# Replace values
PEgenes['PGenes'] = PEgenes['PGenes'].map(lambda x: orthdictP.get(x, x))
PEgenes['EGenes'] = PEgenes['EGenes'].map(lambda x: orthdictE.get(x, x))

# Select all values orthologous in both columns
Ancestor = PEgenes.loc[(PEgenes['PGenes'].str.contains('gene-HOLleu_')) & 
                       (PEgenes['EGenes'].str.contains('gene-HOLleu_'))]

Ancestor = Ancestor.rename(columns = {'Chr' : 'Chromosome',
                                      'PGenes' : 'Name', 
                                      'A' : 'Pchr',
                                      'B' : 'Echr'})
Ancestor = Ancestor[['Chromosome', 'Name', 'Pchr', 'Echr']]
Ancestor['Hchr'] = Ancestor.loc[:, 'Name']

# Add column with native sea cucumber chromosome
Hol = Holleu.to_numpy()
orthdictHchr = dict(zip(Hol[:, 3], Hol[:, 0]))
Ancestor['Hchr'] = Ancestor['Hchr'].map(lambda x: orthdictHchr.get(x, x))

# Make ancestor BED file
Ancestorgenelist = pd.DataFrame()
Ancestorgenelist['Chromosome'] = Ancestor.loc[:, 'Chromosome']
Ancestorgenelist['Start'] = Holleu.loc[:, 'Start']
Ancestorgenelist['End'] = Holleu.loc[:, 'End']
Ancestorgenelist['Name'] = Ancestor.loc[:, 'Name']
Ancestorgenelist['Dot'] = Holleu.loc[:, 'Dot']

np.savetxt(r'Data/Genelists/Ancestor.genelist.bed', Ancestorgenelist.values, fmt = '%s')

BilChr = BilAnc.loc[BilAnc['Chromosome'].isin(['G', 'B1', 'B2', 'M', 'C2', 'B3','P',
                                       'L', 'F', 'J1', 'R', 'O2', 'N', 'A2',
                                       'H', 'J2', 'C1', 'D', 'K', 'I', 'O1'])]

Ancestorgenelist = pd.concat([BilChr, Ancestorgenelist])

np.savetxt(r'Data/Genelists/AniAnc.genelist.bed', Ancestorgenelist.values, fmt = '%s')