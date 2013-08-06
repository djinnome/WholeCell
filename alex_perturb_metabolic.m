function growth = alex_perturb_metabolic(perturbVec)
% This is an un-finished function. I think for this to be useful, we would
% want to include something a little more complicated. Either have the
% protein copy numbers increase exponentially over time, or run the
% translation model with the metabolic model.

rng('shuffle')

tuID_enzymeId_rxnId = {
    'TU_003' 'MG_006_DIMER' 'Tmk'
    'TU_011' 'MG_023_DIMER' 'Fba'
    'TU_027' 'MG_047_TETRAMER' 'MetK'
    'TU_069' 'MG_111_DIMER' 'Pgi'
    'TU_180' 'MG_271_272_273_274_192MER' 'AceE'
    'TU_203' 'MG_299_DIMER' 'Pta'
    'TU_233' 'MG_330_MONOMER' 'CmkA2'
    'TU_260' 'MG_357_DIMER' 'AckA'
    'TU_294' 'MG_407_DIMER' 'Eno'
    'TU_307' 'MG_431_DIMER' 'TpiA'
    };

% Get stuff:
setPath();
setWarnings();
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rnaHalfLives = sim.getRnaHalfLives();
rxnKinetics = sim.getMetabolicReactionKinetics();
%get handle to metabolism sub-model and metabolic reaction state
met = sim.process('Metabolism');
mr = sim.state('MetabolicReaction');

a = 1;
for i_gene = 1:10
  tuId = tuID_enzymeId_rxnId{i_gene, 1};
  enzymeId = tuID_enzymeId_rxnId{i_gene, 2};
  rxnId = tuID_enzymeId_rxnId{i_gene, 3};
  % Calculate proportional change in promotor affinity times RNA half life
  init = rnaPolTuBindingProbs.(tuId) * rnaHalfLives.(tuId);
  final = (perturbVec(a)*rnaPolTuBindingProbs.(tuId)) * (perturbVec(a+1)*rnaHalfLives.(tuId));
  propChng = (final-init)/init;
  
  % Perturb enzyme copy numbers by proportional change
  enzyme_index = strcmp(met.enzymeWholeCellModelIDs, enzymeId);
  met.enzymes(enzyme_index) = propChng * met.enzymes(enzyme_index);
  % Reaction kCat Parameter
  sim.applyMetabolicReactionKinetics(struct(rxnId, struct('for', perturbVec(a+2) * rxnKinetics.(rxnId).for)));
  a = a+3;
end

seed = randi(220516568);
sim.applyOptions('seed', seed);

parameterVals = sim.getAllParameters();
filename = sprintf('metabolic-sim-%06i-%09i.mat',length(dir('./metabolic-output'))-2,seed);

    
%simulate dynamics for 10s - It shouldn't change over time.
lengthSec = 10;
growth = nan(lengthSec, 1);
for i = 1:lengthSec
    met.evolveState();
    growth(i) = mr.growth;
end
