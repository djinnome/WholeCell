function f000_predictPerturbs_pairwiseCorr(datatype)

% Load rnaArray.mean for all simulations.
basedir = 'C:\Users\Alex\Documents\GitHub\WholeCell\';
outputdir = [basedir 'output\'];
unmod = load([basedir 'individual_perturbations\averaged_output\averaged_sim-0.mat']);
W = what(outputdir);
good_ind = unmod.(datatype).std(:) ~= 0;
data_stacked = nan(length(W.mat)+1,sum(good_ind)); % matrix holding all rnaArray data
pV_stacked = nan(30,length(W.mat)); % matrix holding all perturbations

w_this = what('.');
if any(strcmp(w_this.mat,['f000_' datatype '.mat']))
  load(['f000_' datatype '.mat']); % load pV_stacked, data_stacked
else
  for a = 1:length(W.mat)
    S = load([outputdir W.mat{a}]);
    mean_vec = S.(datatype).mean(:);
    data_stacked(a,:) = mean_vec(good_ind); % ignore variables with std=0
    if ~isfield(S,'perturbVec')
      disp(['perturbVec missing: ' W.mat{a}])
    else
      pV_stacked(:,a) = S.perturbVec;
    end
  end
end
% Add the mutant data to the data_stacked matrix
S = load([basedir 'mutant_WT_avgData']);
mean_vec = S.(datatype).mean(:);
data_stacked(end,:) = mean_vec(good_ind);


% Normalize rnaArray data as z-scores. Use PCA to reduce dimensionality
data_normlzd = zscore(data_stacked);
[Coeff,Score,~,~,Explained] = pca(data_normlzd);

save(['./compiled_data/f000_' datatype '.mat'],'pV_stacked',...
                      'data_stacked','Coeff','Score','Explained');

figure;
subplot(2,1,1)
plot(Explained(1:40),'ob')
xlabel('Principal Component Number')
ylabel('Variance Explained')
title(datatype)

% Maximum r-squared value is still very low after PCA
subplot(2,1,2); hold on
max_rsq = 0;
for a = 1:25
  for p = 1:3:28
    y = prod(pV_stacked(p:(p+2),:))'; % Try to predict the product of perturbations
    x = Score(1:(end-1),a);
    r_sqr = corr(x,y)^2;
    if r_sqr > max_rsq
      cla
      plot(x,y,'ob')
      max_rsq = r_sqr;
      xlabel(['Principal Component #' num2str(a)])
      ylabel(['Perturbation to Gene #' num2str(ceil(p/3))])
      ylimits = get(gca,'ylim');
      plot([Score(end,a) Score(end,a)], ylimits,'--g','LineWidth',1.5)
    end
  end
end
h = lsline;
set(h,'Color','red');
title(['Greatest predictive power: R^2 = ' num2str(max_rsq)])
