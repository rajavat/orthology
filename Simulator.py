# Imports
import numpy as np
import pandas as pd
import argparse
import os

import random
from random import randrange

# Disable chained assignments
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
        
def fusion(genome, chr, mixing = 0):
    '''
    inputs: 
    ancestor : df with chromosome name | gene name
    mixing : float between 0 and 1, where 1 implies extreme mixing and 0 implies no mixing
    '''
    
    # Randomly select two chromosomes to fuse
    A = random.choice(chr)
    B = random.choice(chr)
    
    chr = [x for x in chr if x not in (A, B)]
    
    if A == B: # Just so the same chromosome isn't selected twice
        B = random.choice(chr)

    fusion = ancestor.loc[ancestor['Chr'].isin([A, B])]
    
    # Apply mixing if required
    if mixing > 0:
        genes = fusion['Genes'].to_numpy()
        n = len(genes)
        for i in range(int(mixing * n)):
            g1, g2 = randrange(n), randrange(n)
            genes[g2], genes[g1] = genes[g1], genes[g2]

        fusion['Genes'] = genes
        fusion['Chr'] = f'{A}x{B}'
        log = f'Fusion of AncChr{A} and AncChr{B} into Chr{A}x{B}'
        
    else:
         fusion['Chr'] = f'{A}+{B}'
         log = f'Fusion of AncChr{A} and AncChr{B} into Chr{A}+{B}'
    
    # Remove the unfused chromosomes
    genome.drop(genome[genome['Chr'].isin([A, B])].index, inplace = True)
    genome = pd.concat([genome, fusion])
    
    return genome, log, chr

def fission(genome, chr):
    # Randomly select a chromosome for fission
    A = random.choice(chr)
    fission = genome.loc[genome['Chr'] == A]
    chr.remove(A)
    
    pos = random.choice(range(1, Ngene))

    # Add the new chromosomes back into the genome
    chr1 = fission.iloc[: pos]
    chr1['Chr'] = f'{A}_1'
    
    chr2 = fission.iloc[pos :]
    chr2['Chr'] = f'{A}_2'
    
    # Remove the fission chromosome from the genome
    genome = pd.concat([genome, chr1, chr2])
    genome = genome[genome.Chr != A]
    
    log = f'Fission of AncChr{A} into Chr{A}_1 and Chr{A}_2'
    
    return genome, log, chr

def translocation(genome, chr):
    # Randomly select two chromosomes for translocation
    A = random.choice(chr)
    B = random.choice(chr)
    
    chr = [x for x in chr if x not in (A, B)]
    
    chrA = genome.loc[genome['Chr'] == A]
    chrB = genome.loc[genome['Chr'] == B]
    
    # Randomly select two break point positions
    posA = random.choice(range(1, Ngene))
    posB = random.choice(range(1, Ngene))
    
    # Join the fragments to form recombinant chromosomes
    chr1 = pd.concat([chrA.iloc[: posA], chrB.iloc[posB :]])
    chr1['Chr'] = f'{A};{B}'
    chr2 = pd.concat([chrB.iloc[: posB], chrA.iloc[posA :]])
    chr2['Chr'] = f'{B};{A}'
    
    # Remove the original chromosomes from the genome
    genome = pd.concat([genome, chr1, chr2]).drop(genome[(genome['Chr'] == A) & (genome['Chr'] == B)].index)
    
    log = f'Translocation between AncChr{A} and AncChr{B}'
    
    return genome, log, chr

def syntenyloss(genome, chr):
    A = random.choice(chr)
    syn = genome.loc[genome['Chr'] == A]
    genome = genome[genome.Chr != syn]
    
    chr.remove(A)
    
    # Assign all elements to a random chromosome
    syn['Chr'] = random.choices(genome.Chr.unique(), k = len(syn))
    
    # Add back into the genome
    genome = pd.concat([genome, syn])
    
    log = f'Synteny loss of AncChr{A}'

    return genome, log, chr

# Apply macro-rearrangements to the ancestor
ancestor = makeancestor(Nchr, Ngene)
chr = ancestor.Chr.unique().tolist()
speciesA = ancestor.copy()

events = {}
for event in range(args['Nevents']):
    r = np.random.uniform()
    
    if r <= 0.30:
        if len(ancestor) < 2: continue
        speciesA, log, chr = fission(speciesA, chr)
        events['EVENT_' + str(event + 1)] = log
        print(log)
    
    elif r <= 0.45:
        speciesA, log, chr = translocation(speciesA, chr)
        events['EVENT_' + str(event + 1)] = log
        print(log)
    
    elif r <= 0.70:
        speciesA, log, chr = fusion(speciesA, chr)
        events['EVENT_' + str(event + 1)] = log
        print(log)
    
    elif r <= 0.95:
        speciesA, log, chr = fusion(speciesA, chr, mixing = 0.5)
        events['EVENT_' + str(event + 1)] = log
        print(log)
        
    else:
        speciesA, log, chr = syntenyloss(speciesA, chr)
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