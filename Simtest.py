import numpy as np
import Orthoscripts

for i in range(100+1):
    input = 'Simulations/Ancestor_' + str(i) + '.bed'
    ancestor = Orthoscripts.readBED(input, 's')
    input = 'Simulations/SpeciesA_' + str(i) + '.bed'
    speciesA = Orthoscripts.readBED(input, 's')
    input = 'Simulations/Ancestor+SpeciesA_' + str(i) + '.txt'
    orthos = np.loadtxt(input, dtype = "str")

    data = Orthoscripts.orthologies(ancestor, speciesA, orthos)
    
    outfile = 'Simulations/Rearrangements_' + str(i) + '.txt'
    Orthoscripts.rearrangements(data, outfile)