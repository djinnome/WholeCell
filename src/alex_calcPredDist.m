function predDist = alex_calcPredDist(tru,est)

% calculate means for singleCell params
time = min(length(tru.singleCell.growth), length(est.singleCell.growth));
tru.growth.mean = mean(tru.singleCell.growth(:,1:time),1)';
tru.growth.std = std(tru.singleCell.growth(:,1:time),0,1)';
tru.mass.mean = mean(tru.singleCell.mass(:,1:time),1)';
tru.mass.std = std(tru.singleCell.mass(:,1:time),0,1)';
tru.volume.mean = mean(tru.singleCell.volume(:,1:time),1)';
tru.volume.std = std(tru.singleCell.volume(:,1:time),0,1)';

est.growth.mean = mean(est.singleCell.growth(:,1:time),1)';
est.mass.mean = mean(est.singleCell.mass(:,1:time),1)';
est.volume.mean = mean(est.singleCell.volume(:,1:time),1)';

dSeq = calcSubDist(tru.dnaSeq.mean, est.dnaSeq.mean, tru.dnaSeq.std);
mCs = calcSubDist(tru.metConcs.mean, est.metConcs.mean, tru.metConcs.std);
pA = calcSubDist(tru.protArray.mean, est.protArray.mean, tru.protArray.std);
rA = calcSubDist(tru.rnaArray.mean, est.rnaArray.mean, tru.rnaArray.std);
rSeq = calcSubDist(tru.rnaSeq.mean, est.rnaSeq.mean, tru.rnaSeq.std);
flx = calcSubDist(tru.rxnFluxes.mean, est.rxnFluxes.mean, tru.rxnFluxes.std);
gr = calcSubDist(tru.growth.mean, est.growth.mean, tru.growth.std);
ma = calcSubDist(tru.mass.mean, est.mass.mean, tru.mass.std);
vol = calcSubDist(tru.volume.mean, est.volume.mean, tru.volume.std);

% Combine all variables
allvars = [dSeq; mCs; pA; rA; rSeq; flx; gr; ma; vol]; 
predDist = (1/length(allvars)) * sum(allvars);

function dists = calcSubDist (mean_tru, mean_est, std_tru)
std_tru(std_tru == 0) = NaN; %min(std_tru(std_tru~=0)); % avoid divide by zero
dists = ((mean_tru - mean_est)./std_tru).^2;
dists = dists(~isnan(dists));
