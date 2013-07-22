% This script attempts to build a neural network that aims to predict
% elements of perturbVec based on rnaArray data.

close all; clear; clc
set(0,'DefaultFigureWindowStyle','docked')

cd ..
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
cd analysis

% Load rnaArray.mean for all simulations.
basedir = 'C:\Users\Alex\Documents\Projects\DREAM8\WholeCell-parameter-estimation-DREAM-challenge-2013\';
outputdir = [basedir 'output\'];
unmod = load([basedir 'individual_perturbations\averaged_output\averaged_sim-0.mat']);
W = what(outputdir);
rnaArr_stacked = nan(length(W.mat),525); % matrix holding all rnaArray data
pV_stacked = nan(30,length(W.mat)); % matrix holding all perturbations

for a = 1:length(W.mat)
  S = load([outputdir W.mat{a}]);
  rnaArr_stacked(a,:) = S.rnaArray.mean;
  if ~isfield(S,'perturbVec')
    disp(['perturbVec missing: ' W.mat{a}])
  else
    pV_stacked(:,a) = S.perturbVec;
  end
end

% Normalize rnaArray data as z-scores. Use PCA to reduce dimensionality
rnaArr_norm = zscore(rnaArr_stacked,0,2);
[Coeff,Score,~,~,Explained] = pca(rnaArr_norm);

figure;
subplot(3,1,1)
plot(Coeff(:,1),'-b')
title('First Principal Component of RNA array')
subplot(3,1,2)
plot(Coeff(:,2),'-r')
title('Second Principal Component of RNA array')
subplot(3,1,3)
plot(Coeff(:,3),'-k')
title('Third Principal Component of RNA array')
xlabel('Different Genes')

figure;
plot(Explained(1:40))
xlabel('Principal Component Number')
ylabel('Variance Explained')

% Maximum r-squared value is still very low after PCA
figure;
max_rsq = 0;
for a = 1:50
  for p = 1:3:27
    y = prod(pV_stacked(p:(p+3),:))'; % Try to predict the product of perturbations
    x = Score(:,a);
    r_sqr = corr(x,y)^2;
    if r_sqr > max_rsq
      plot(x,y,'ob')
      max_rsq = r_sqr;
    end
  end
  disp(a)
end
h = lsline;
set(h,'Color','red');
xlabel('Principal Component Score')
ylabel('Perturbation')
title('Greatest predictive power')
