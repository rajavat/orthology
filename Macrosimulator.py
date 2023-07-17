import argparse
import os

import numpy as np

import random
from random import randrange


def make_ancestor(c, g):
	ancestor = []
	for i in range(c):
		tmp = list(range(1+g*i, 1+(i+1)*g))
		ancestor.append(tmp)
	chrom_names = ['AncChr'+str(i+1) for i in range(c)]
	return chrom_names, ancestor

def write_ancestor(chrom_names, ancestor, outfile):
	with open(outfile, 'w') as out:
		for i, chrom in enumerate(chrom_names):
			k = 1
			for g in ancestor[i]:
				out.write('\t'.join([chrom, str(k), str(k+1), 'ancg_'+str(g)])+'\n')
				k += 1

def write_dummy_orthologs_file(ancestor, outfile):
	with open(outfile, 'w') as out:
		o = 1
		for i, chrom in enumerate(ancestor):
			for g in chrom:
				out.write('\t'.join(['ortholog_'+str(o), 'ancg_'+str(g), 'g_'+str(g)])+'\n')
				o += 1

def write_evolved_ancestor(ancestor, outfile):
	with open(outfile, 'w') as out:
		for i, chrom in enumerate(ancestor):
			k = 1
			chrom_name = 'Chr'+str(i+1)
			for g in chrom:
				out.write('\t'.join([chrom_name, str(k), str(k+1),'g_'+str(g)])+'\n')
				k += 1

def apply_fusion(chrom_names, genome, mixing=False, mixfactor=1):

	#randomly select two chromosomes and remove them from the genome
	idx1 = random.choice(range(len(genome)))
	chrom1 = genome[idx1]
	name1 = chrom_names[idx1]

	del genome[idx1]
	del chrom_names[idx1]

	idx2 = random.choice(range(len(genome)))
	chrom2 = genome[idx2]
	name2 = chrom_names[idx2]

	del genome[idx2]
	del chrom_names[idx2]

	chrom_fuse = chrom1 + chrom2

	newname = f'{name1}+{name2}'
	#apply mixing if requested
	if mixing:
		mix(chrom_fuse, mixfactor)
		newname = f'{name1}x{name2}'

	#Put the chromosomes back into the genome
	genome = genome.append(chrom_fuse)

	#Update chromosome list

	chrom_names.append(newname)

	#Return log info about the fusion event
	log_info = f'Fusion of {name1} and {name2} into {newname}'
	return log_info

def mix(chrom, mixfactor=1): #mixfactor should be a float between 0 and 1, where 1 implies extreme mixing and 0 no mixing
	n = len(chrom)
	for i in range(int(mixfactor*n)):
		g1, g2 = randrange(n), randrange(n)
		chrom[g2], chrom[g1] = chrom[g1], chrom[g2]

def apply_fission(chrom_names, genome):

	#randomly select a chromosome and remove it from the genome and chrom names
	idx = random.choice(range(len(genome)))
	chrom = genome[idx]
	name = chrom_names[idx]

	del genome[idx]
	del chrom_names[idx]

	#randomly select a fission position
	g = random.choice(range(1, len(chrom)))

	#apply fission
	chroms = [chrom[:g], chrom[g:]]

	#Put the chromsomoes back into the genome
	genome = genome + chroms

	newname1 = name + '_1'
	newname2 = name + '_2'

	#Update the list of chromosome names
	chrom_names.append(newname1)
	chrom_names.append(newname2)

	#Return log info about the fission event
	log_info = f'Fission {name} into {newname1} and {newname2}'
	return log_info

# def synteny_loss(chrom, all_chromsomes):
# 	return

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-c', '--nb_chr', type=int, required=False, default=20, help="Number of chromosomes in the ancestor")

	parser.add_argument('-g', '--nb_genes', type=int, required=False, default=100, help="Number of genes on each chromosomes")

	parser.add_argument('-n', '--nb_macro_events', type=int, required=False, default=10,
						help="Number of macro-syntenic rearrangement events, selected amongst [fission, fusion, fusion_with_mixing, synteny_loss] with proba [0.48, 0.33, 0.15, 0.04]")

	# parser.add_argument('-r', '--mixing-events-ratio', type=str, required=False, default=5, help="Number of intra-chromosomal")

	parser.add_argument('-p', "--out_prefix", type=str, required=False, default='simulations/sim_test/')

	args = vars(parser.parse_args())

	c = args['nb_chr']
	g = args['nb_genes']

	#Create output folders
	if '/' in args['out_prefix']:
		folder = '/'.join(args['out_prefix'].split('/')[:-1])
		os.makedirs(folder, exist_ok=True)

	#Create and write .bed file and ortholog files for the ancestor
	anc_chr, ancestor = make_ancestor(c, g)

	outfile = args['out_prefix'] + 'ancestor.geneslist.bed'
	write_ancestor(anc_chr, ancestor, outfile)

	outfile = args['out_prefix'] + 'ancestor+speciesA.txt'
	write_dummy_orthologs_file(ancestor, outfile)

	#Apply macro-rearrangments to the ancestor
	list_of_events = {}
	for event in range(args['nb_macro_events']):
		r = np.random.uniform()

		if r <= 0.48:
			if len(ancestor) < 2: continue
			event_log = apply_fission(anc_chr, ancestor)
			list_of_events['EVENT_'+str(event+1)] = event_log
			print(event_log)

		elif r <= 0.81:
			event_log = apply_fusion(anc_chr, ancestor)
			list_of_events['EVENT_'+str(event+1)] = event_log
			print(event_log)

		elif r <= 0.94:
			event_log = apply_fusion(anc_chr, ancestor, mixing=True)
			list_of_events['EVENT_'+str(event+1)] = event_log
			print(event_log)

		else:
			print('Synteny loss: skipping for now...')
			# need to implement the synteny loss (and we could also have a WGD event)
			continue

	#Write .bed file for Species A evolved from ancestor
	outfile = args['out_prefix'] + 'speciesA.geneslist.bed'
	write_evolved_ancestor(ancestor, outfile)

	outfile = args['out_prefix'] + 'list_of_rearrangements.txt'
	with open(outfile, 'w') as out:
		for event in list_of_events:
			out.write(f'{event}: {list_of_events[event]}\n')
	# write_event_log(list_of_events, outfile)