import sys # use to access arguments
import os # use in order to call commands from terminal script is called in
import re # regular expressions used to 


# python genome2bed.py [genome fasta file] [output file name]

#read in file name i want to change
filename = sys.argv[1] 

#name of the output file we wish to create, then activly create it
output = sys.argv[2] 
os.system('touch '+output) 

file = open(filename, "r")
out = open(output, "w")



length_of_sequence = 0 # a variable to ultimately hold the length of the sequence
iteration = 0 #a counter that keeps certain logical things on track
name = '' # initialize the name variable to hold the chromosome name


#loop through every line in the genome fasta
for line in file:
    #record the name variableon the first line
    if line[0] == '>' and iteration == 0:
        name = line[1:].rstrip()
    #everytime we encounter a new label, write out to output file
    #reset length of sequence
    elif line[0] == '>' and iteration > 0:
        out.write(name+'\t0\t'+str(length_of_sequence)+'\n')
        #remove any trailing white spaces from the line.
        name = line[1:].rstrip()
        length_of_sequence = 0
    #keep adding the length of each line otherwise    
    else:
        length_of_sequence = len(line) + length_of_sequence
    #increment iteration by 1    
    iteration += 1
    
    
def genbed(filename, outname):
    file = open(filename, "r")
    out = open(outname, "w")
    
    os.system('touch '+outname)  # create outfile
    
    seqlen = 0 # a variable to ultimately hold the length of the sequence
    iteration = 0 #a counter that keeps certain logical things on track
    name = '' # initialize the name variable to hold the chromosome name
    
    #loop through every line in the genome fasta
    for line in file:
        #record the name variableon the first line
        if line[0] == '>' and iteration == 0:
            name = line[1:].rstrip()
        #everytime we encounter a new label, write out to output file
        #reset length of sequence
        elif line[0] == '>' and iteration > 0:
            out.write(name+'\t0\t'+str(seqlen)+'\n')
            #remove any trailing white spaces from the line.
            name = line[1:].rstrip()
            seqlen = 0
        #keep adding the length of each line otherwise    
        else:
            seqlen = len(line) + seqlen
        #increment iteration by 1    
        iteration += 1
        
        return out
        