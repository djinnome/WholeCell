function genes = genes_of_operon( sim, tus )
  rna = sim.state('Rna');
  for (i = 1:length(tus))
    tu = char(tus(i));
    genes.(tu) =  sim.gene.wholeCellModelIDs(find(rna.nascentRNAGeneComposition(:, strcmp(rna.wholeCellModelIDs(rna.nascentIndexs), tu))))';
  end