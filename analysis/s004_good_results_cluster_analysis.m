close all; clear; clc;

load('../whole-cell-table.mat') % Load results from cluster.

N_best = round(length(whole_cell_table.number)*0.2); % take top 10 percent

bestParams = whole_cell_table.perturbVecs(:,1:N_best)';
bestDists = whole_cell_table.predDist(1:N_best);

% Vizualize clusters using pca
[COEFF, SCORE] = pca(bestParams);

% PCs
figure; hold on
plot(COEFF(:,1)*100,'-ob')
plot(COEFF(:,2)*100,'-or')
plot(COEFF(:,3)*100,'-og')
xlim([0 31])
xlabel('Different Perturbations')
ylabel('PCA coefficients (centered)')
legend('PC #1','PC #2','PC #3')

% Scores - Not seeing any easily identifiable clusters..
figure;
plot3(SCORE(:,1),SCORE(:,2),SCORE(:,3),'o','MarkerFaceColor','b')
xlabel('Score PC #1')
ylabel('Score PC #2')
zlabel('Score PC #3')
grid on

% Scores - Not seeing any easily identifiable clusters..
figure;
for a = 1:4
  subplot(2,2,a)
  plot(SCORE(:,a),bestDists,'x')
  h = lsline;
  set(h,'Color','r')
  xlabel(['Score PC #' num2str(a)])
  ylabel('Prediction Distance')
end

% 3D color scatter
figure;
color_scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),bestDists);
grid on