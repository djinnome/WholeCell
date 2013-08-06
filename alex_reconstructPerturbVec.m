function perturbVec2 = alex_reconstructPerturbVec(simParams)

% Parameter IDs
genesTusRxns = {
    'MG_006' 'TU_003' 'Tmk'
    'MG_023' 'TU_011' 'Fba'
    'MG_047' 'TU_027' 'MetK'
    'MG_111' 'TU_069' 'Pgi'
    'MG_272' 'TU_180' 'AceE'
    'MG_299' 'TU_203' 'Pta'
    'MG_330' 'TU_233' 'CmkA2'
    'MG_357' 'TU_260' 'AckA'
    'MG_407' 'TU_294' 'Eno'
    'MG_431' 'TU_307' 'TpiA'
    };
  
perturbVec2 = nan(30,1); % reconstructed perturbVec

sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();

rnaPolTuBindingProbs_unmod = sim.getRnaPolTuBindingProbs();
rxnKinetics_unmod = sim.getMetabolicReactionKinetics();
rnaHalfLives_unmod = sim.getRnaHalfLives();

sim.applyAllParameters(simParams);

rnaPolTuBindingProbs_mod = sim.getRnaPolTuBindingProbs();
rxnKinetics_mod = sim.getMetabolicReactionKinetics();
rnaHalfLives_mod = sim.getRnaHalfLives();


a = 1;
for i_gene = 1:10
  tuId = genesTusRxns{i_gene, 2};
  rxnId = genesTusRxns{i_gene, 3};
  % Promotor Affinity Parameter
  perturbVec2(a) = rnaPolTuBindingProbs_mod.(tuId)/rnaPolTuBindingProbs_unmod.(tuId);
  % Half Life Parameter
  perturbVec2(a+1) = rnaHalfLives_mod.(tuId)/rnaHalfLives_unmod.(tuId);
  % Reaction kCat Parameter
  perturbVec2(a+2) = rxnKinetics_mod.(rxnId).for/rxnKinetics_unmod.(rxnId).for;
  a = a+3;
end