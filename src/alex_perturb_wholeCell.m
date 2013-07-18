function filename = alex_perturb_wholeCell(perturbVec)
assert(length(perturbVec) == 30);
rng('shuffle')

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

setPath();
setWarnings();
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rnaHalfLives = sim.getRnaHalfLives();
rxnKinetics = sim.getMetabolicReactionKinetics();

a = 1;
for i_gene = 1:10
  tuId = genesTusRxns{i_gene, 2};
  rxnId = genesTusRxns{i_gene, 3};
  % Promotor Affinity Parameter
  sim.applyRnaPolTuBindingProbs(struct(tuId, perturbVec(a) * rnaPolTuBindingProbs.(tuId)));
  % Half Life Parameter
  sim.applyRnaHalfLives(struct(tuId, perturbVec(a+1) * rnaHalfLives.(tuId)));
  % Reaction kCat Parameter
  sim.applyMetabolicReactionKinetics(struct(rxnId, struct('for', perturbVec(a+2) * rxnKinetics.(rxnId).for)));
  a = a+3;
end

seed = randi(220516568);
sim.applyOptions('lengthSec', 65000, 'seed', seed);

parameterVals = sim.getAllParameters();
filename = sprintf('whole-cell-sim-%06i-%09i.mat',length(dir('./output'))-2,seed);
simulateHighthroughputExperiments(...
    'seed', seed, ...
    'parameterVals', parameterVals, ...
    'simPath', ['output/' filename] ...
    );