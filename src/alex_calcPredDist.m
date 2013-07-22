function predDist = alex_calcPredDist(tru,est)
% alex_calcPredDist: Estimates the "prediction distance" - the distance
%                    between all of the state variables (e.g. RNA array
%                    concentrations) for two whole cell simulations. See
%                    section 2.7.2 of the competition description.
%
% INPUT:  "tru" - a struct containing all the "high-throughput data" for a
%                 simulation. This is the "true" or "target" parameter set 
%                 that we are trying to replicate.
%
%         "est" - a struct containing all the "high-throughput data" for an
%                 estimated parameter set
%
% OUTPUT: "predDist" - a double that represents the prediction distance
%
% Author: Alex Williams

% calculate means for singleCell params (they aren't given to us)
tru.growth.mean = nanmean(tru.singleCell.growth,1)';
tru.growth.std = nanstd(tru.singleCell.growth,0,1)';
tru.mass.mean = nanmean(tru.singleCell.mass,1)';
tru.mass.std = nanstd(tru.singleCell.mass,0,1)';
tru.volume.mean = nanmean(tru.singleCell.volume,1)';
tru.volume.std = nanstd(tru.singleCell.volume,0,1)';

est.growth.mean = nanmean(est.singleCell.growth,1)';
est.mass.mean = nanmean(est.singleCell.mass,1)';
est.volume.mean = nanmean(est.singleCell.volume,1)';

% Calculate the sum of the squared z-scores for each of the data types
dSeq = calcSubDist(tru.dnaSeq.mean, est.dnaSeq.mean, tru.dnaSeq.std);
mCs = calcSubDist(tru.metConcs.mean, est.metConcs.mean, tru.metConcs.std);
pA = calcSubDist(tru.protArray.mean, est.protArray.mean, tru.protArray.std);
rA = calcSubDist(tru.rnaArray.mean, est.rnaArray.mean, tru.rnaArray.std);
rSeq = calcSubDist(tru.rnaSeq.mean, est.rnaSeq.mean, tru.rnaSeq.std);
flx = calcSubDist(tru.rxnFluxes.mean, est.rxnFluxes.mean, tru.rxnFluxes.std);
gr = calcSubDist(tru.growth.mean, est.growth.mean, tru.growth.std);
ma = calcSubDist(tru.mass.mean, est.mass.mean, tru.mass.std);
vol = calcSubDist(tru.volume.mean, est.volume.mean, tru.volume.std);
chipS = calcSubDist(tru.chipSeq.mean(:), est.chipSeq.mean(:), tru.chipSeq.std(:));

% Sum all of the partial sums (computed directly above), and divide by the
% total number of variables that we are considering.
allvars = [dSeq; mCs; pA; rA; rSeq; flx; gr; ma; vol; chipS]; 
predDist = (1/length(allvars)) * sum(allvars);

% This sub-function computes the distance for any given data type. For
% example, this function will calculate the distance between tru and est
% for the rnaSeq or protArray sub-struct.
function dists = calcSubDist (mean_tru, mean_est, std_tru)
assert(length(mean_tru) == length(mean_est));
ind = (std_tru ~= 0); % avoid divide by zero
dists = ((mean_tru(ind) - mean_est(ind))./std_tru(ind)).^2;
dists(isnan(dists)) = 0; % might be an unnecessary safeguard
