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
    parser.add_argument('-r', '--Nruns', type = int, 
                        required = False, default = 1,
                        help = "Number of runs")
    args = vars(parser.parse_args())
    
    Nchr = args['Nchr']
    Ngene = args['Ngene']
    Nevents = args['Nevents']
    Nruns = args['Nruns']

# Create output folder
os.makedirs(os.path.dirname('Simulations/'), exist_ok=True)

def filename(path):
    i = 1
    while os.path.exists(path % i):
        i += 1
    return path % i
    
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
def dummyBED(genome, type, outfile):
    if type == 'anc':
        genome['Chromosome'] = 'AncChr' + genome['Chromosome'].astype(str)
        genome['Name'] = 'ancg_' + genome['Name'].astype(str)
        
    if type == 'des':
        genome['Chromosome'] = 'Chr' + genome['Chromosome'].astype(str)
        genome['Name'] = 'g_' + genome['Name'].astype(str)
    
    genome['Start'] = np.arange(len(genome))
    genome['End'] = np.arange(len(genome)) + 5
    
    genome = genome[['Chromosome', 'Start', 'End', 'Name']]
    
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
    genes = genome['Name'].to_numpy()
    n = len(genes)
    for i in range(int(mixing * n)):
        g1, g2 = randrange(n), randrange(n)
        genes[g2], genes[g1] = genes[g1], genes[g2]

        genome['Name'] = genes
        # genome['Chromosome'] = f'{fuse1}x{fuse2}'
        
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
    
    else:
        chr1 = pd.concat([chrA.iloc[: posA]])
        chr1['Chromosom'] = f'{A};'
        chr2 = pd.concat([chrB, chrA.iloc[posA :]])
        chr2['Chromosome'] = f'{B};{A}'
    
        log = ['TR', 
               ('AncChr' + str(A) + ', AncChr' + str(B)),
               ('Chr' + str(B) + ';' + str(A))]
        
    # Remove the original chromosomes from the genome
    genome = pd.concat([genome, chr1, chr2]).drop(genome[(genome['Chromosome'] == A) & (genome['Chromosome'] == B)].index)
    return genome, log, chr

def syntenyloss(genome, chr):
    A = random.choice(chr)
    syn = genome.loc[genome['Chromosome'] == A]
    genome = genome[genome.Chromosome != syn]
    
    chr.remove(A)
    
    # Assign all elements to a random chromosome
    syn['Chromosome'] = random.choices(genome.Chromosome.unique(), k = len(syn))
    
    # Add back into the genome
    genome = pd.concat([genome, syn])
    
    log = f'Synteny loss of AncChr{A}'

    return genome, log, chr

# Apply macro-rearrangements to the ancestor
for i in range(Nruns):
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
            print(log)
        
        elif r <= 0.45:
            speciesA, log, chr = translocation(speciesA, chr)
            events.append(log)
            print(log)
        
        elif r <= 0.70:
            speciesA, log, chr = fusion(speciesA, chr)
            events.append(log)
            print(log)
        
        elif r <= 0.99:
            speciesA, log, chr = fusion(speciesA, chr, mixing = 0.5)
            events.append(log)
            print(log)
            
        else:
            # speciesA, log, chr = syntenyloss(speciesA, chr)
            # events.append(log)
            # print(log)
            continue
        
    # Create BED files and orthology file
    outfile = 'Simulations/SpeciesA_' + str(i) + '.bed'
    dummyBED(speciesA, 'des', outfile)
    outfile = 'Simulations/Ancestor_' + str(i) + '.bed'
    dummyBED(ancestor, 'anc', outfile)
    outfile = 'Simulations/Ancestor+SpeciesA_' + str(i) + '.txt'
    dummyOrthologs(ancestor, outfile)

    # Create list of rearrangements
    outfile = 'Simulations/Events_' + str(i) + '.txt'
    with open(outfile, 'w') as out:
        for event in events:
            # write each item on a new line
            out.write("%s\n" % event)