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
                ancestor = pd.concat([ancestor, pd.DataFrame([row])], ignore_index = True)
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
        
def fusion(genome, mixing = 0):
    '''
    inputs: 
    ancestor : df with chromosome name | gene name
    mixing : float between 0 and 1, where 1 implies extreme mixing and 0 implies no mixing
    '''
    
    # Randomly select two chromosomes to fuse
    fuse1 = random.choice(genome.Chr.unique())
    fuse2 = random.choice(genome.Chr.unique())
    
    if fuse1 == fuse2: # Just so the same chromosome isn't selected twice
        fuse2 = random.choice(range(1, len(genome.Chr.unique())))

    fusion = ancestor.loc[ancestor['Chr'].isin([fuse1, fuse2])]
    
    # Apply mixing if required
    if mixing > 0:
        genes = fusion['Genes'].to_numpy()
        n = len(genes)
        for i in range(int(mixing * n)):
            g1, g2 = randrange(n), randrange(n)
            genes[g2], genes[g1] = genes[g1], genes[g2]

        fusion['Genes'] = genes
        fusion['Chr'] = f'{fuse1}x{fuse2}'
        
    else:
         fusion['Chr'] = f'{fuse1}+{fuse2}'
    
    # Remove the unfused chromosomes
    genome.drop(genome[genome['Chr'].isin([fuse1, fuse2])].index, inplace = True)
    genome = pd.concat([genome, fusion])
    
    log = f'Fusion of AncChr{fuse1} and AncChr{fuse2} into Chr{fuse1}+{fuse2}'
    
    return genome, log

def fission(genome):
    # Randomly select a chromosome for fission
    fiss = random.choice(genome.Chr.unique())
    fission = genome.loc[genome['Chr'] == fiss]

    pos = random.choice(range(1, Ngene))

    # Add the new chromosomes back into the genome
    chr1 = fission.iloc[: pos]
    chr1['Chr'] = f'{fiss}_1'
    chr2 = fission.iloc[pos :]
    chr2['Chr'] = f'{fiss}_2'
    
    # Remove the fission chromosome from the genome
    genome = pd.concat([genome, chr1, chr2])
    genome = genome[genome.Chr != fiss]
    
    log = f'Fission of AncChr{fiss} into Chr{fiss}_1 and Chr{fiss}_2'
    
    return genome, log

def translocation(genome):
    # Randomly select two chromosomes for translocation
    cA = random.choice(genome.Chr.unique())
    cB = random.choice(genome.Chr.unique())
    
    chrA = genome.loc[genome['Chr'] == cA]
    chrB = genome.loc[genome['Chr'] == cB]
    
    # Randomly select two break point positions
    posA = random.choice(range(1, Ngene))
    posB = random.choice(range(1, Ngene))
    
    # Join the fragments to form recombinant chromosomes
    chr1 = pd.concat([chrA.iloc[: posA], chrB.iloc[posB :]])
    chr1['Chr'] = f'{cA};{cB}'
    chr2 = pd.concat([chrB.iloc[: posB], chrA.iloc[posA :]])
    chr2['Chr'] = f'{cA};{cB}'
    
    # Remove the original chromosomes from the genome
    genome = pd.concat([genome, chr1, chr2]).drop(genome[(genome['Chr'] == cA) & (genome['Chr'] == cB)].index)
    
    log = f'Translocation between AncChr{cA} and AncChr{cB}'
    
    return genome, log


def synteny_loss(genome):
    syn = random.choice(genome.Chr.unique())
    synchr = genome.loc[genome['Chr'] == syn]
    genome = genome[genome.Chr != syn]
    
    # Assign all elements to a random chromosome
    synchr['Chr'] = random.choices(genome.Chr.unique(), k = len(synchr))
    
    # Add back into the genome
    genome = genome.append(synchr)
    
    log = f'Synteny loss of AncChr{syn}'

    return genome, log

# Apply macro-rearrangements to the ancestor
ancestor = makeancestor(Nchr, Ngene)
speciesA = ancestor.copy()

events = {}
for event in range(args['Nevents']):
    r = np.random.uniform()
    
    if r <= 0.40:
        if len(ancestor) < 2: continue
        speciesA, log = fission(speciesA)
        events['EVENT_' + str(event + 1)] = log
        print(log)
    
    elif r <= 0.70:
        speciesA, log = fusion(speciesA)
        events['EVENT_' + str(event + 1)] = log
        print(log)
    
    else:
        speciesA, log = fusion(speciesA, mixing = 0.5)
        events['EVENT_' + str(event + 1)] = log
        print(log)
        continue
    
# Create BED files and orthology file
outfile = args['out_prefix'] + 'SpeciesA.genelist.bed'
dummyBED(speciesA, 'des', outfile)
outfile = args['out_prefix'] + 'Ancestor.genelist.bed'
dummyBED(ancestor, 'anc', outfile)
outfile = args['out_prefix'] + 'Ancestor+SpeciesA.txt'
dummyOrthologs(ancestor, outfile)

# Create list of rearrangements
outfile = args['out_prefix'] + 'rearrangements.txt'
with open(outfile, 'w') as out:
    for event in events:
        out.write(f'{event}: {events[event]}\n')