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
    chrNames = ['AncChr' + str (i + 1) for i in range(Nchr)]
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

    # Randomly select a fission position and apply fission
    n = len(fission['Chr'])

    pos = random.choice(range(1, Ngene))
    genes = fission['Genes'].to_numpy()

    # Add the new chromosomes back into the genome
    chr1 = pd.DataFrame({'Chr': [f'{fiss}_1'] * len(genes[:pos]),
                        'Genes': genes[:pos]})
    chr2 = pd.DataFrame({'Chr': [f'{fiss}_2'] * len(genes[pos:]),
                       'Genes': genes[pos:]})

    speciesA = ancestor.append([chr1, chr2])
    
    # Remove the fission chromosome from the genome
    speciesA.drop(ancestor[ancestor['Chr'] == fiss].index, inplace = True)
    
    log = f'Fission of AnchChr{fiss} into Chr{fiss}_1 and Chr{fiss}_2'
    
    return speciesA, log

# Apply macro-rearrangements to the ancestor
ancestor = makeancestor(Nchr, Ngene)
speciesA = ancestor.copy(deep = True)
events = {}
for event in range(args['Nevents']):
    r = np.random.uniform()
    
    if r <= 0.48:
        if len(ancestor) < 2: continue
        speciesA, event_log = fission(speciesA)
        events['EVENT_' + str(event + 1)] = event_log
        print(event_log)
    
    elif r <= 0.81:
        speciesA, event_log = fission(speciesA)
        events['EVENT_' + str(event + 1)] = event_log
        print(event_log)
    
    elif r <= 0.94:
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