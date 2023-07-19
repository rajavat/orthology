import numpy as np 
import argparse
import random 
from random import randrange

# Input arguments
if __name__ == '__main__':
    parser = argparse.ArgumentParser(descriptiom = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-c', '--nChr', type = int, 
                        required = False, default = 20, 
                        help = 'Number of chromosomes in the ancestor')
    
    parser.add_argument('-g', '--nGenes', type = int,
                        requried = False, default = 100, 
                        help = 'Number of genes on each chromosome')
    
    parser.add_argument('-n', '--nMacroEvents', type = int, 
                        required = False, default = 10, 
                        help = 'Number of macro-syntetic rearrangement events')
    
    parser.add_argument('-p', '--out_prefix', type = str,
                        required = False, default = 'Simulations/')
    
    args = vars(parser.parse_args())
    
    chr = args['nChr']
    genes = args['nGenes']

# Functions
def Ancestor(chr, genes): 
    ancestor = []
    for i in range(chr):
        tmp = list(range(1 + genes * 1), (1 + (i + 1) * genes))
        ancestor.append(tmp)
    chrNames = ['AncChr' + str (i + 1) for i in range(chr)]
    return chrNames, ancestor

# File-handling functions
def writeAncestor(chrNames, ancestor, outfile): # Write dummy BED file for ancestor
    with open(outfile, 'w') as out:
        for i, chr in enumerate(chrNames):
            k = 1
            for g in ancestor[i]:
                out.write('\t'.join(chr, str(k), str(k + 1), 'ancg_' + str(g)) + '\n')
                k += 1
                
def writeDescendant(ancestor, outfile): # Write dummy BED file for descendant
    with open(outfile, 'w') as out:
        for i, chr, in enumerate(ancestor):
            k = 1
            chrNames = 'Chr' + str(i + 1)
            for g in chr:
                out.write('\t'.join([chrNames, str(k), str(k + 1), 'g_' + str(g)]) + '\n')
                k += 1
                
def writeOrthologs(ancestor, outfile): # Write dummy ortholog file
    with open(outfile, 'w') as out:
        o = 1
        for i, chrom in enumerate(ancestor):
            for g in chrom:
                out.write('\t'.join(['ortholog_' + str(o), 'ancg_' + str(g), 'g_' + str(g)]) + '\n')
                o += 1
    
# Rearrangement functions
def mixing(chr, mixfactor = 1): # 'mixfactor' is a float between 0 and 1, where 1 implies extreme mixing and 0 implies no mixing
    n = len(chr)
    for i in range(int(mixfactor * n)):
        g1, g2 = randrange(n), randrange(n)
        chr[g2], chr[g1] = chr[g1], chr[g2] # Genes from g1, g2 are mixed
        
def fusion(chrNames, genome, mixing = False):
    # Randomly select two chromosomes and remove them from the genome
    rc1 = random.choice(range(len(genome)))
    chr1 = genome[rc1]
    name1 = chrNames[rc1]
    
    del chr1
    del name1
    
    rc2 = random.choice(range(len(genome)))
    chr2 = genome[rc2]
    name2 = chrNames[rc2]
    
    del chr2
    del name2
    
    chrFuse = chr1 + chr2
    
    newname = f'{name1} + {name2}'
    
    # Apply mixing if requested
    if mixing: 
        mixing(chrFuse)
        newname = f'{name1} x {name2}'
    
    # Put the chromosomes back into the genome
    genome = genome.append(chrFuse)
    
    # Update chromosomes list 
    chrNames.append(newname)
    
    # Return event log
    log = f'Fusion of {name1} and {name2} into {newname}'
    return log
 
def fission(chrNames, genome):
    # Randomly select a chromosome and remove it from the genome
    rc = random.choice(range(len(genome)))
    chr = genome[rc]
    name = chrNames[rc]
    
    del chr
    del name 
    
    # Randomly select a fission position
    g  = random.choice(range(1, len(chr)))
    
    # Apply fission
    chrom = [chr[:g], chr[g:]] 
    
    # Put the chromosomes back into the genome
    genome = genome + chrom 
    
    newname1 = name + '_1'
    newname2 = name + '_2'
    
    # Update the list of chromosome names
    chrNames.append(newname1)
    chrNames.append(newname2)
    
    # Return log info about the fission event 
    log = f'Fission {name} into {newname1} and {newname2}'
    return log 

# def WGD(chrNames, genome): 

# def syntenyloss(chrNames, genome):