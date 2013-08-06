function alex_compareToTarget(filenames)
% Loads a simulation in the "output" folder and compares it to the mutant
% wildtype data that was provided to us for free.

N = length(filenames);
basedir = 'C:\Users\Alex\Documents\GitHub\WholeCell\';
for a = 1:N
  est(a) = load([basedir 'output\' filenames{a}]);
end
targ = load([basedir 'mutant_WT_avgData.mat']);

% singleCell
figure(1); hold on
y_est = nan(65001,N);
for a = 1:N
  y_est(:,a) = est(a).singleCell.volume';
end
y_est(y_est==0) = nan;
y_targ = targ.singleCell.volume'; y_targ(y_targ==0) = nan;
plot(y_targ,'-b')
plot(y_est,'-r','LineWidth',1)
xlabel('time (s)')
ylabel('Volume')
title('30 different trials (end of lines = cell cycle length)')
%pos = get(gcf,'Position');
%set(gcf,'PaperPosition',[1 3 6.5 5],'PaperUnits','Inches','PaperSize',[8.5 11])
%saveas(gcf,'singleCell.pdf')

