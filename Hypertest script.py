def hypertest(dataset, speciesA, speciesB, orthologies):
    """
    inputs:
    chrA: species A chromosome name 
    chrB: species B chromosome name
    speciesA: column name in dataset for species A
    speciesB: column name in dataset for species B
        orthologies: orthology dataset
    """
    
    A = A_AH_data.Name.values.tolist()
    B = H_AH_data.Name.values.tolist()
    M = len(list(set(A) & set(B)))
    
    nA = dataset.loc[(dataset[speciesA] == chrA), 'Orthologs'].sum()
    nB = dataset.loc[(dataset[speciesB] == chrB), 'Orthologs'].sum()
    x = dataset.loc[(dataset[speciesA] == chrA) & (dataset[speciesB] == chrB), 'Orthologs'].sum()
    
    p = stats.hypergeom.sf(x - 1, M, nA, nB)
     
