% This script attempts to predict paramters of simluations from the rnaSeq
% highthroughput data. This code was written to give me some insight into
% the data by playing around with different plots. It is mostly here to
% keep a record of what I did, and is not particularly useful anymore.

close all; clear; clc
set(0,'DefaultFigureWindowStyle','docked')
datatype = 'rnaArray';

cd ..
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
cd analysis

% Load rnaArray.mean for all simulations.
basedir = 'C:\Users\Alex\Documents\Projects\DREAM8\WholeCell-parameter-estimation-DREAM-challenge-2013\';
outputdir = [basedir 'output\'];
unmod = load([basedir 'individual_perturbations\averaged_output\averaged_sim-0.mat']);
W = what(outputdir);
data_stacked = nan(length(W.mat),length(unmod.(datatype).mean)); % matrix holding all rnaArray data
pV_stacked = nan(30,length(W.mat)); % matrix holding all perturbations

for a = 1:length(W.mat)
  S = load([outputdir W.mat{a}]);
  data_stacked(a,:) = S.(datatype).mean;
  if ~isfield(S,'perturbVec')
    disp(['perturbVec missing: ' W.mat{a}])
  else
    pV_stacked(:,a) = S.perturbVec;
  end
end

% Normalize rnaArray data as z-scores. Use PCA to reduce dimensionality
data_normlzd = zscore(data_stacked,0,2);
[Coeff,Score,~,~,Explained] = pca(data_normlzd);

figure;
subplot(3,1,1)
plot(Coeff(:,1),'-b')
title(['First Principal Component of ' datatype])
subplot(3,1,2)
plot(Coeff(:,2),'-r')
title(['Second Principal Component of ' datatype])
subplot(3,1,3)
plot(Coeff(:,3),'-k')
title(['Third Principal Component of ' datatype])
xlabel('Different Genes')

figure;
plot(Explained(1:40),'ob')
xlabel('Principal Component Number')
ylabel('Variance Explained')

% Maximum r-squared value is still very low after PCA
figure;
max_rsq = 0;
for a = 1:25
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
title(['Greatest predictive power: R^2 = ' num2str(max_rsq)])
