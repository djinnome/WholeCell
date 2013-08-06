clear; close all; clc;

% Set up the model
setPath();
setWarnings();
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();

rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rxnKinetics = sim.getMetabolicReactionKinetics();

% Alter Promoter Affinities
%sim.applyRnaPolTuBindingProbs(struct('TU_011', 0.5 * rnaPolTuBindingProbs.('TU_011')));

% Alter kcats
sim.applyMetabolicReactionKinetics(struct('Pgi', struct('for', 0.5 * rxnKinetics.('Pgi').for)));

% Save the parameters locally
simName = 'Pgi_kcat_x05';
parameters = sim.getAllParameters();

% Submit to the cloud
bucketUrl = 's3://alexander.williams.1967383.wcpe.sagebase.org';
[jobId, status, errMsg] = postCloudSimulation(...
    'simName', simName, ...
    'bucketUrl', bucketUrl, ...
    'parameterVals', parameters ...
    );

save('job-13-07-26b');