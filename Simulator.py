# Imports
import numpy as np
import pandas as pd
import argparse
import os

import random
from random import randrange

# Disable chained assignments
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None 

# Inputs
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', '--Nchr', type = int, 
                        required = False, default = 20,
                        help = 'Number of chromosome in the ancestor')
    parser.add_argument('-g', '--Ngene', type = int, 
                        required = False, default = 100,
                        help = "Number of genes on each chromosome")
    parser.add_argument('-n', '--Nevents', type = int, 
                        required = False, default = 10,
                        help = "Number of macro-syntetic rearrangement events")
    parser.add_argument('-p', "--out_prefix", type=str, 
                        required=False, default='Simulations/')
    args = vars(parser.parse_args())
    
    Nchr = args['Nchr']
    Ngene = args['Ngene']

# Create output folder
if '/' in args['out_prefix']:
    folder = '/'.join(args['out_prefix'].split('/')[:-1])
    os.makedirs(folder, exist_ok = True)
    
# Make ancestor genome
def makeancestor(Nchr, Ngene):
    ancestor = pd.DataFrame(columns = ['Chr'])
    for i in range(Nchr):
        row = {'Chr' : (i + 1)}
        for i in range(Ngene):
                ancestor = ancestor.append(row, ignore_index = True)
    ancestor['Genes'] = (ancestor.reset_index().index + 1)
    # ancestor['Genes'] = 'g_' + ancestor['Genes'].astype(str)

    return ancestor

# Dummy BED files :: type 'anc' for ancestor, 'des' for descendant
def dummyBED(genome, type, outfile):
    if type == 'anc':
        genome['Chr'] = 'AncChr' + genome['Chr'].astype(str)
        genome['Genes'] = 'ancg_' + genome['Genes'].astype(str)
        
    if type == 'des':
        genome['Chr'] = 'Chr' + genome['Chr'].astype(str)
        genome['Genes'] = 'g_' + genome['Genes'].astype(str)
    
    genome['Start'] = np.arange(len(genome))
    genome['End'] = np.arange(len(genome)) + 5
    
    genome = genome[['Chr', 'Start', 'End', 'Genes']]
    
    with open(outfile, 'w') as out:
        out.write(genome.to_string(header = False, index = False))
        
    return genome

# Dummy ortholog file
def dummyOrthologs(genome, outfile):
    
    orthologs = pd.DataFrame()
    
    orthologs['Orthologs'] = np.arange(len(genome)) + 1
    orthologs['speciesA'] = np.arange(len(genome)) + 1
    orthologs['speciesB'] = np.arange(len(genome)) + 1
    
    orthologs['Orthologs'] = 'orthologs_' + orthologs['Orthologs'].astype(str)
    orthologs['speciesA'] = 'ancg_' + orthologs['speciesA'].astype(str)
    orthologs['speciesB'] = 'g_' + orthologs['speciesB'].astype(str)
    
    with open(outfile, 'w') as out:
        out.write(orthologs.to_string(header = False, index = False))

def mixing(genome, mixing):
    genes = genome['Genes'].to_numpy()
    n = len(genes)
    for i in range(int(mixing * n)):
        g1, g2 = randrange(n), randrange(n)
        genes[g2], genes[g1] = genes[g1], genes[g2]

        genome['Genes'] = genes
        # genome['Chr'] = f'{fuse1}x{fuse2}'
        
def fusion(ancestor, mixing = 0):
    '''
    inputs: 
    ancestor : df with chromosome name | gene name
    mixing : float between 0 and 1, where 1 implies extreme mixing and 0 implies no mixing
    '''
    
    # Randomly select two chromosomes to fuse
    fuse1 = random.choice(range(len(ancestor.Chr.unique())))
    fuse2 = random.choice(range(len(ancestor.Chr.unique())))

    # Fuse the chromosomes
    fusion = ancestor.loc[ancestor['Chr'].isin([fuse1, fuse2])]
    fusion['Chr'] = f'{fuse1}+{fuse2}'
    
    # Apply mixing if required
    if mixing > 0:
        genes = fusion['Genes'].to_numpy()
        n = len(genes)
        for i in range(int(mixing * n)):
            g1, g2 = randrange(n), randrange(n)
            genes[g2], genes[g1] = genes[g1], genes[g2]

        fusion['Genes'] = genes
        fusion['Chr'] = f'{fuse1}x{fuse2}'

    # Add the fused chromosome back into the genome
    speciesA = ancestor.append([ancestor, fusion])
    
    # Remove the unfused chromosomes
    speciesA.drop(speciesA[speciesA['Chr'].isin([fuse1, fuse2])].index, inplace = True)
    
    log = f'Fusion of AncChr{fuse1} and AncChr{fuse2} into Chr{fuse1}+{fuse2}'
    
    return speciesA, log

def fission(ancestor):
    # Randomly select a chromosome for fission
    fiss = random.choice(range(len(ancestor.Chr.unique())))
    fission = ancestor.loc[ancestor['Chr'] == fiss]

    pos = random.choice(range(1, Ngene))

    # Add the new chromosomes back into the genome
    chr1 = fission.iloc[: pos]
    chr1['Chr'] = f'{fiss}_1'
    chr2 = fission.iloc[pos :]
    chr2['Chr'] = f'{fiss}_2'
    
    # Remove the fission chromosome from the genome
    speciesA = ancestor.append([chr1, chr2])
    speciesA = speciesA[speciesA.Chr != fiss]
    
    log = f'Fission of AncChr{fiss} into Chr{fiss}_1 and Chr{fiss}_2'
    
    return speciesA, log

def translocation(ancestor):
    # Randomly select two chromosomes for translocation
    cA = random.choice(range(len(ancestor.Chr.unique())))
    cB = random.choice(range(len(ancestor.Chr.unique())))
    
    chrA = ancestor.loc[ancestor['Chr'] == cA]
    chrB = ancestor.loc[ancestor['Chr'] == cB]
    
    # Randomly select two break point positions
    posA = random.choice(range(1, Ngene))
    posB = random.choice(range(1, Ngene))
    
    # Break up the two chromosomes into four fragments
    fragmentA1 = chrA.iloc[: posA]
    fragmentA2 = chrA.iloc[posA :]
    fragmentB1 = chrB.iloc[: posB]
    fragmentB2 = chrB.iloc[posB :]
    
    # Join the fragments to form recombinant chromosomes
    chr1 = fragmentA1.append(fragmentB2)
    chr1['Chr'] = f'{cA};{cB}'
    chr2 = fragmentA2.append(fragmentB1)
    chr2['Chr'] = f'{cA};{cB}'
    
    # Remove the original chromosomes from the genome
    speciesA = ancestor.append([chr1, chr2])
    
    speciesA = speciesA[speciesA.Chr != cA]
    speciesA = speciesA[speciesA.Chr != cB]
    
    log = f'Translocation between AncChr{cA} and AncChr{cB}'
    
    return speciesA, log

def synteny_loss(ancestor):
    syn = random.choice(range(len(ancestor.Chr.unique())))
    synchr = ancestor.loc[ancestor['Chr'] == syn]

    speciesA = ancestor.copy() # Create species A
    speciesA = speciesA[speciesA.Chr != syn] # Remove lost chromosome
    
    # Assign all elements to a random chromosome
    synchr['Chr'] = random.choices(ancestor.Chr.unique(), k = len(synchr))
    
    # Add back into the genome
    speciesA = speciesA.append(synchr)
    
    log = f'Synteny loss of AncChr{syn}'

    return speciesA

# Apply macro-rearrangements to the ancestor
ancestor = makeancestor(Nchr, Ngene)
speciesA = ancestor.copy(deep = True)
events = {}
for event in range(args['Nevents']):
    r = np.random.uniform()
    
    if r <= 0.40:
        if len(ancestor) < 2: continue
        speciesA, event_log = fission(speciesA)
        events['EVENT_' + str(event + 1)] = event_log
        print(event_log)
    
    elif r <= 0.60:
        speciesA, event_log = translocation(speciesA)
        events['EVENT_' + str(event + 1)] = event_log
        print(event_log)
    
    elif r <= 0.80:
        speciesA, event_log = fusion(speciesA)
        events['EVENT_' + str(event + 1)] = event_log
        print(event_log)
    
    elif r <= 0.90:
        speciesA, event_log = fusion(ancestor, mixing = 0.5)
        events['EVENT_' + str(event + 1)] = event_log
        print(event_log)
        
    else:
        print('Synteny loss: skipping for now...')
        # Need to implement the synteny loss (and we could also have a WGD event)
        continue
    
# Create BED files and orthology file
outfile = args['out_prefix'] + 'Ancestor.genelist.bed'
dummyBED(ancestor, 'anc', outfile)
outfile = args['out_prefix'] + 'Ancestor+SpeciesA.txt'
dummyOrthologs(ancestor, outfile)
outfile = args['out_prefix'] + 'SpeciesA.genelist.bed'
dummyBED(speciesA, 'des', outfile)

# Create list of rearrangements
outfile = args['out_prefix'] + 'rearrangements.txt'
with open(outfile, 'w') as out:
    for event in events:
        out.write(f'{event}: {events[event]}\n')