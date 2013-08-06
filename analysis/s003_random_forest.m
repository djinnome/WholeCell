% This script tries to use the "random forest" method from machine learning
% to predict the value of each parameter in the competition.

% Need all the ./compiled_data/f000_*.mat files to be made. 
% This can be done by running s002 script.

clear; close all; clc

datatypes = {'metConcs'
             'dnaSeq'
             'rnaSeq'
             'chipSeq'
             'rnaArray'
             'protArray'
             'rxnFluxes'};

all_data_stacked = [];
for a = 1:length(datatypes)
  S = load(['./compiled_data/f000_' datatypes{a} '.mat']);
  all_data_stacked = cat(2,all_data_stacked,S.data_stacked);
  %sPCA.(datatypes{a}).Score = S.Score(:,1:25);
  %sPCA.(datatypes{a}).Score = S.Explained(1:25);
  pV_stacked = S.pV_stacked;
end

minCols = min(all_data_stacked);
rngCols = range(all_data_stacked);
nr = size(all_data_stacked,1);
data_normlzd = zscore(all_data_stacked);
[Coeff,Score,~,~,Explained] = pca(data_normlzd);

% Figure suggesting that ~25 principal components and ~400 trees is enough
figure; hold all
for nComp = [5 25 100 length(Explained)]
  X = Score(1:(end-1),1:nComp);
  y = pV_stacked(3,:)';
  
  opts = statset('UseParallel','always');
  B = TreeBagger(5e2,X,y,'Method','regression','oobvarimp','on','Options',opts);
  plot(oobError(B))
  xlabel('number of grown trees')
  ylabel('out-of-bag classification error')
end
legend('5','10','15','25','50','100')

% Figure showing the random forest predictions for each variable:
figure; hold on
nComp = [5 50 100];
colors = ['b' 'r' 'k'];
offsets = linspace(-0.2,0.2,length(nComp));
nTrees = 500;
opts = statset('UseParallel','always');
for b = 1:length(nComp)
  perturbVec_est = nan(30,1);
  for param = 1:30
    y = pV_stacked(param,:)';
    X = Score(1:(end-1),1:nComp(b));
    x_pred = Score(end,1:nComp(b)); % The target 'mutant' phenotype
    
    B = TreeBagger(nTrees,X,y,'Method','regression','oobvarimp','on','Options',opts);
    y_est = B.predict(x_pred);
    y_err = oobError(B); y_err = y_err(end);
    
    plot(param+offsets(b),100*y_est,'o','Color',colors(b))
    plot((param+offsets(b))*ones(1,2),[(y_est-y_err); (y_est+y_err)],'-','Color',colors(b))
    perturbVec_est(param) = mean(y_est);
  end
  save(['./compiled_data/s003_perturbVec_est_nComp' num2str(nComp(b))],'perturbVec_est')
end
xlabel('Different Parameters')
ylabel('Estimate of Perturbation (% change)')
title('Whole Cell Parameters (1-30)')