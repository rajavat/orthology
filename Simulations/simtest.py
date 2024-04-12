import numpy as np
import syntiny

for i in range(100+1):
    input = 'Simulations/Ancestor_' + str(i) + '.bed'
    ancestor = syntiny.readBED(input, 's')
    input = 'Simulations/SpeciesA_' + str(i) + '.bed'
    speciesA = syntiny.readBED(input, 's')
    input = 'Simulations/Ancestor+SpeciesA_' + str(i) + '.txt'
    orthos = np.loadtxt(input, dtype = "str")

    data = syntiny.orthologies(ancestor, speciesA, orthos)
    
    outfile = 'Simulations/Rearrangements_' + str(i) + '.txt'
    syntiny.rearrangements(data, outfile)