function paramFileName = perturbBest( sim, bestFile, baseDir)

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

sim = CachedSimulationObjectUtil.load();

bestParameters = load(bestFile);
sim.applyAllParameters(bestParameters);

rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rnaHalfLives = sim.getRnaHalfLives();
rxnKinetics = sim.getMetabolicReactionKinetics();

for i = 1:size(genesTusRxns, 1)
    for j = 1:numel(parameterTypes)
        for k = 1:size(parameterVals, 1)
            sim.applyAllParameters(bestParameters);
            
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
            
            parameters = sim.getAllParameters();
            paramFileName{i,j,k} = fullfile(baseDir, sprintf('parameters_%s_%s_%s.mat', genesTusRxns{i, 1}, parameterTypes{j}, parameterVals{k, 1}));
            save(paramFileName{i,j,k}, '-struct', 'parameters');
        end
    end
end
