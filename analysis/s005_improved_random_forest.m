clear; close all; clc

datatypes = {'metConcs'
             'dnaSeq'
             'rnaSeq'
             'chipSeq'
             'rnaArray'
             'protArray'}; % 'rxnFluxes' left out, zscore is bad..

all_data_stacked = [];
for a = 1:length(datatypes)
  S = load(['./compiled_data/f000_' datatypes{a} '.mat']);
  all_data_stacked = cat(2,all_data_stacked,S.data_stacked);
  %sPCA.(datatypes{a}).Score = S.Score(:,1:25);
  %sPCA.(datatypes{a}).Score = S.Explained(1:25);
  pV_stacked = S.pV_stacked;
end
S = load('./compiled_data/f000_rxnFluxes.mat');

data_normlzd = zscore(all_data_stacked,0,2);
[Coeff,Score,~,~,Explained] = pca(data_normlzd);

x1 = Score(1:(end-1),1:5); % first five principal components seem to work.
x2 = S.data_stacked(1:(end-1),:);
x3 = [pV_stacked(1:2,:); pV_stacked(4:end,:)]';
x1_pred = Score(end,1:5); % The target 'mutant' phenotype
x2_pred = S.data_stacked(end,:);
x3_pred = ones(29,1);

% Figure comparing random forests at 3 levels of complexity
figure; hold all
y = pV_stacked(3,:)';

subplot(2,1,1)
B1 = TreeBagger(5e2,x1,y,'Method','regression','oobvarimp','on');
B2 = TreeBagger(5e2,cat(2,x1,x2),y,'Method','regression','oobvarimp','on');
B3 = TreeBagger(5e2,cat(2,x1,x2,x3),y,'Method','regression','oobvarimp','on');
plot(oobError(B1))
plot(oobError(B2))
plot(oobError(B3))
xlabel('number of grown trees')
ylabel('out-of-bag classification error')
legend('without rxn flux','with rxn flux','with rxn flux & other perts')

subplot(2,1,2)
y_est = B1.predict(x1_pred);

save('./compiled_output/s005_treebags','B1','B2','B3')
