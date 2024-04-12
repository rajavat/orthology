#!/usr/bin/env python


"""
This script parses blast outputs to extract reciprocal best hits.
Results are written in ancgene format.

Example:
    $ python get_recblast_ortho.py -c1 blast/Salmo_trutta+Salvelinus_namaycush.blp \
-c2 blast/Salvelinus_namaycush+Salmo_trutta.blp -o orthologs_Salmo_trutta-Salvelinus_namaycush.txt

"""

import sys
import argparse

def parse_best_hit(input_file):

    """
    Parses a blast output file '.blb' to extract all best macthes from g1 to g2.

    Args:
        input_file (str): path to the blast outfile

    Returns:
        dict : for each protein in g1 (key) its best hit in g2 (value)
    """

    best_hits = {}
    with open(input_file, 'r') as infile:
        for line in infile:

            line = line.strip().split('\t')
            gene1, gene2, pident = line[:3]
            pident = float(pident)
            evalue = float(line[-2])

            if gene1 not in best_hits:
                best_hits[gene1] = (gene2, pident)

            # #usually blast hits are sorted by evalue, but let's check anyway
            elif pident > best_hits[gene1][1]: 
                best_hits[gene1] = (gene2, pident)

    return best_hits


def write_reciprocal_best_hits(bbh_1, bbh_2, output_file):

    """
    Compares best blast hits g1 vs g2 in `bbh_1` and best blast hits g2 vs g1 in `bbh_2`, and writes
    reciprocal best hits to file, in ancGene format.

    Args:
        bbh_1, bbh_2 (dicts): best blasts with key = prot in g_ref and
                              value = tuple (best prot in g_query, pident).
        output_fil (str): path to the output file.
    """
    i = 1
    with open(output_file, 'w') as out:
        for prot1 in bbh_1:
            prot2, _ = bbh_1[prot1]
            if prot2 not in bbh_2:
                continue
            if bbh_2[prot2][0] == prot1:
                out.write(f'ortholog_{i} {prot1} {prot2}\n')
                i += 1

#            else:
#                print(prot2, prot1, bbh_2[prot2], bbh_1[prot1])


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    REQUIRED = PARSER.add_argument_group('required named arguments')

    REQUIRED.add_argument('-c1', '--comp1', type=str, help="First blast comparison", required=True)

    REQUIRED.add_argument('-c2', '--comp2', type=str, help="Second blast comparison", required=True)

    REQUIRED.add_argument('-o', '--output', type=str, help='Output file', required=True)

    ARGS = vars(PARSER.parse_args())

    sys.stderr.write(f"Parsing {ARGS['comp1']}...")
    sys.stderr.flush()
    BBH1 = parse_best_hit(ARGS["comp1"])
    sys.stderr.write(f"done\n")
    sys.stderr.write(f"Parsing {ARGS['comp2']}...")
    sys.stderr.flush()
    BBH2 = parse_best_hit(ARGS["comp2"])
    sys.stderr.write(f"done\n")

    sys.stderr.write(f"Writing reciprocal best hits...")
    sys.stderr.flush()
    write_reciprocal_best_hits(BBH1, BBH2, ARGS["output"])
    sys.stderr.write(f"done\n")



