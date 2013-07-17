setWarnings();
setPath();


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

% fluxRatio = struct('TpiA', 0.10950078509522165, 'Tmk', 0.67000000000000004, ...
%                    'AckA', 0.10907225186548129, 'CmkA2', 0.66999999999999071, ...
%                    'AceE', 0.10908082824118392, 'Eno', 0.10908084168275271, ...
%                    'Fba', 0.10946107633831005, 'Pta', 0.10907225186548129, ...
%                    'MetK', 0.67000000000000015, 'Pgi', 0.11105379245342617); ...
    
fluxRatio = struct('TpiA', 0.12262036324436351, 'Tmk', 0.75000000000000011, ...
                   'AckA', 0.12214087397029504, 'CmkA2', 0.74999999999997413, ...
                   'AceE', 0.12215047370175966, 'Eno', 0.12215105103105892, ...
                   'Fba', 0.12257591632129075, 'Pta',  0.12214087397029504,   ...
                   'MetK', 0.75000000000000011, 'Pgi', 0.12435882605630455);      

parameterTypes = {
    'PromAffinity'
    'HalfLife'
    'RxnKcat'
    };
baseDir = '/msc/neurospora/src/WholeCell/input';
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();

originalParams = sim.getAllParameters();

rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rnaHalfLives = sim.getRnaHalfLives();
rxnKinetics = sim.getMetabolicReactionKinetics();

newRnaPolTuBindingProbs = struct();
newRnaHalfLives = struct();
newRxnKinetics = struct();
rxns = fieldnames( fluxRatio );


sim.applyAllParameters(originalParams);
for i = 1:size(genesTusRxns, 1)
  rxnId = genesTusRxns{i, 3};
  tuId = genesTusRxns{i, 2 };
  setfield(newRnaPolTuBindingProbs,tuId,getfield(fluxRatio,rxnId)* ...
                         rnaPolTuBindingProbs.(tuId));
  setfield(newRnaHalfLives,tuId,getfield(fluxRatio,rxnId)* ...
                        rnaHalfLives.(tuId));
  setfield(newRxnKinetics,rxnId,struct('for',getfield(fluxRatio,rxnId)*rxnKinetics.(rxnId).for));
end

for j = 1:numel(parameterTypes)
  sim.applyAllParameters(originalParams);
  switch parameterTypes{j}
   case 'PromAffinity'
    sim.applyRnaPolTuBindingProbs( newRnaPolTuBindingProbs );
   case 'HalfLife'
    sim.applyRnaHalfLives( newRnaHalfLives );
   case 'RxnKcat'
    sim.applyMetabolicReactionKinetics( newRxnKinetics );
  end
  parameters = sim.getAllParameters();
  paramFileName = fullfile(baseDir,sprintf('parameters_fluxRatio_%s.mat',parameterTypes{j}));
  save(paramFileName, '-struct', 'parameters')
  
end

