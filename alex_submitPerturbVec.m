function alex_submitPerturbVec(perturbVec)

% Check inputs
assert(length(perturbVec) == 30);
assert(all(perturbVec > 0));

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


% Set up the model
setPath();
setWarnings();
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();

rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rxnKinetics = sim.getMetabolicReactionKinetics();
rnaHalfLives = sim.getRnaHalfLives();


% The index "a" keeps track of the parameter (within a gene) that is being
% perturbed {order = 'promoter affinity','rna half life','kcat'}. This
% index is incremented by three every iteration of the loop.
a = 1;
for i_gene = 1:10
  tuId = genesTusRxns{i_gene, 2};
  rxnId = genesTusRxns{i_gene, 3};
  % Promotor Affinity Parameter
  sim.applyRnaPolTuBindingProbs(struct(tuId,...
                               perturbVec(a)*rnaPolTuBindingProbs.(tuId)));
  % Half Life Parameter
  sim.applyRnaHalfLives(struct(tuId,...
                               perturbVec(a+1)*rnaHalfLives.(tuId)));
  % Reaction kCat Parameter
  sim.applyMetabolicReactionKinetics(struct(rxnId,...
                  struct('for', perturbVec(a+2)*rxnKinetics.(rxnId).for)));
  a = a+3;
end

% Generate simName
c = clock;
month = num2str(c(2));
day = num2str(c(3));
if c(3)<10
  day = ['0' day];
end
hour = num2str(c(4));
simName = ['perturbVec_' month '_' day '_' hour];
parameters = sim.getAllParameters();

% Submit to the cloud
bucketUrl = 's3://alexander.williams.1967383.wcpe.sagebase.org';
[jobId, status, errMsg] = postCloudSimulation(...
    'simName', simName, ...
    'bucketUrl', bucketUrl, ...
    'parameterVals', parameters ...
    );

% Save parameters locally for later reference 
save(['job-13-0' month '-' day '-' hour],'jobId','status','errMsg','perturbVec','simName')