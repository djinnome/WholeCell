% This was the script I used to replicate the "in silico" experiments
% described on the competition website. I perturbed each of the parameters
% by 0.5x and 2.0x individually... Resulting in 60 different sets of
% parameters. I ran this script many times until I had enough replicates to
% my satisfaction.

clear; close all; clc;
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
parameterTypes = {
    'PromAffinity'
    'HalfLife'
    'RxnKcat'
    };
parameterVals = {
    '05X' 0.5
    '2X'  2.0
    };

setPath();
setWarnings();
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();

rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rnaHalfLives = sim.getRnaHalfLives();
rxnKinetics = sim.getMetabolicReactionKinetics();

index = randi(60);

i = 1; j = 1; k = 0;
for a = 1:index
  k = k+1;
  if k > size(parameterVals, 1)
    k = 1; j = j+1;
  end
  if j > numel(parameterTypes)
    j = 1; i = i+1;
  end
  if i > size(genesTusRxns, 1)
    error('All simulations have been done!')
  end
end
           
switch parameterTypes{j}
  case 'PromAffinity'
    tuId = genesTusRxns{i, 2};
    sim.applyRnaPolTuBindingProbs(struct(tuId, parameterVals{k, 2} * rnaPolTuBindingProbs.(tuId)));
  case 'HalfLife'
    tuId = genesTusRxns{i, 2};
    sim.applyRnaHalfLives(struct(tuId, parameterVals{k, 2} * rnaHalfLives.(tuId)));
  case 'RxnKcat'
    rxnId = genesTusRxns{i, 3};
    sim.applyMetabolicReactionKinetics(struct(rxnId, struct('for', parameterVals{k, 2} * rxnKinetics.(rxnId).for)));
end

seed = randi(220516568);
sim.applyOptions('lengthSec', 65000, 'seed', seed); % test code by simulating 2 seconds

parameterVals = sim.getAllParameters();
simulateHighthroughputExperiments(...
    'seed', seed, ...
    'parameterVals', parameterVals, ...
    'simPath', ['output/sim-' num2str(index) '_' num2str(seed) '.mat'] ...
    );
averageHighthroughputExperiments(...
   'simPathPattern', ['output/sim-' num2str(index) '_*.mat'], ...
   'avgValsPath', ['averaged_output/averaged_sim-' num2str(index) '.mat'] ...
   );