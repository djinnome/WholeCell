function alex_correlate_metabolic_model
rng('shuffle')
setPath();
setWarnings();

index = 1;
met_growth = nan(60,1); % metabolic growth prediction [mmol/gDCW/s]
cell_growth = nan(60,1); % whole cell growth actual [g/s]
avg_fldr = 'C:\Users\Alex\Documents\Projects\DREAM8\Analysis\averaged_output\';
for i = 1:10
  for j = 1:3
    for k = 1:2
      % simulate metabolic model
      met_growth(index) = run_metabolic_model(i,j,k);
      % load and estimate initial growth in whole-cell model
      S = load([avg_fldr 'averaged_sim-' num2str(index) '.mat']);
      if size(S.singleCell.mass,1)>1
        y = mean(S.singleCell.mass(:,1:1e4));
      else
        y = S.singleCell.mass(:,1:1e4);
      end
      x = 1:length(y);
      P = polyfit(x,y,1);
      cell_growth(index) = P(1);
      % go to next simulation
      index = index+1;
    end
  end
end
figure;
plot(met_growth,cell_growth,'ob')
h = lsline;
set(h,'Color','r')
ylabel('Actual cell growth')
xlabel('Metabolic Predicted Growth')

function predicted_growth = run_metabolic_model(i,j,k)

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
parameterTypes = {
    'PromAffinity'
    'HalfLife'
    'RxnKcat'
    };
parameterVals = {
    '05X' 0.5
    '2X'  2.0
    };

sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
met = sim.process('Metabolism');
mr = sim.state('MetabolicReaction');
rxnKinetics = sim.getMetabolicReactionKinetics();

switch parameterTypes{j}
  case 'PromAffinity'
    % Multiply enzyme copy number instead of promotor affinity
    enzymeId = tuID_enzymeId_rxnId{i, 2};
    enzyme_index = strcmp(met.enzymeWholeCellModelIDs, enzymeId);
    met.enzymes(enzyme_index) = parameterVals{k, 2} * met.enzymes(enzyme_index);
  case 'HalfLife'
    % Multiply enzyme copy number instead of half life
    enzymeId = tuID_enzymeId_rxnId{i, 2};
    enzyme_index = strcmp(met.enzymeWholeCellModelIDs, enzymeId);
    met.enzymes(enzyme_index) = parameterVals{k, 2} * met.enzymes(enzyme_index);
  case 'RxnKcat'
    rxnId = tuID_enzymeId_rxnId{i, 3};
    sim.applyMetabolicReactionKinetics(struct(rxnId, struct('for', parameterVals{k, 2} * rxnKinetics.(rxnId).for)));
end

lengthSec = 10;
growth = zeros(lengthSec, 1);
for i = 1:lengthSec
    met.evolveState();
    growth(i) = mr.growth;
end
predicted_growth = mean(growth);