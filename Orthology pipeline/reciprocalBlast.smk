

# genomes, = glob_wildcards("processed_data/prot/{genomes}.fa")

import os
import itertools

# genomes = ['Afil_fr2py', 'Mar_glac', 'Sti_chlo', 'Parliv', 'Antmed', 'Pecmax', 'Braflo', 'Holleu']

# genomes = ['Strpur', 'Helery']

# genomes = ['Parhyale', 'Ambmex']

# genomes = [ 'Pecmax', 'Braflo', 'Holleu', 'Parliv', 'Afil_fr2py']

# genomes = ['Afil_fr2py', 'Spu3']

genomes = [ 'Pecmax', 'Holleu', 'Braflo', 'Ephmue']


# pairs = [('+').join(pair) for pair in itertools.combinations(genomes, 2) if 'Holleu' in pair]

# print(pairs)
pairs = [('+').join(pair) for pair in itertools.combinations(genomes, 2)]

rule all:
    input:  expand("orthologs/{pair}_sensitive.txt", pair=pairs)

wildcard_constraints:
    spec1 = '[A-Za-z0-9_]+',
    spec2 = '[A-Za-z0-9_]+'

#Afil_lgest is ferdi previous annotation, where I took only the longest trancripts 38,461 genes
#Afil_fr2py is my new annotation, where I took only the longest trancripts 33,725 genes

# rule all:
#     input: "orthologs/Afil_fr2py+Afil_lgest.txt"

rule makedb:
    input: fa = "processed_data/prot/{spec2}.fa"
    output: out = "processed_data/prot/{spec2}.dmnd"
    params: out = lambda w: f"processed_data/prot/{w.spec2}"
    shell: "diamond makedb --in {input.fa} -d {params.out}"


rule blast:
    input:
       sp1 =  "processed_data/prot/{spec1}.fa",
       sp2 =  "processed_data/prot/{spec2}.dmnd"
    output: "blastp/{spec1}"+"+"+"{spec2}.blp"
    threads: 10
    shell:"""
    diamond blastp --query {input.sp1} --db {input.sp2} --ultra-sensitive --evalue 1e-5 --outfmt 6 --out {output} -p {threads}
    """

rule reciprocal_best_blast_hits:
    input:
        comp1 = "blastp/{spec1}+{spec2}.blp",
        comp2 = "blastp/{spec2}+{spec1}.blp"
    output: "orthologs/{spec1}+{spec2}_sensitive.txt"
    shell:"""
    python get_recblast_ortho.py -c1 {input.comp1} -c2 {input.comp2} -o {output}
    """
