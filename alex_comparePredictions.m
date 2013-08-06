function alex_comparePredictions()

clear; close all; clc;


est1 = load('./cloudOutput/perturbVec_7_25_10.predictions.mat');
est2 = load('whole-cell-sim-000185-058982772.mat');
targ = load('gold.predictions.mat');


est1_vec = getPredVec(est1)
est2_vec = getPredVec(est2)
targ_vec = getPredVec(targ)

predDist = alex_calcPredDist(targ,est1)
predDist = alex_calcPredDist(targ,est2)

% figure; hold on
% plot(targ_vec.mean,'-b','LineWidth',2)
% plot(targ_vec.mean+targ_vec.std,'-','Color',[0.7 0.7 1])
% plot(targ_vec.mean-targ_vec.std,'-','Color',[0.7 0.7 1])
% plot(est1_vec.mean,'-','Color',[1 0.7 0.7])
% plot(est2_vec.mean,'-','Color',[0.2 0.7 0.2])
% ylim([-1 1])

function predVec = getPredVec(predStruct)
cellCycleLength = 65000;
predVec.mean = [
                mean(nanmean(predStruct.singleCell.growth, 2), 1)
                mean(nanmean(predStruct.singleCell.mass, 2), 1)
                mean(nanmean(predStruct.singleCell.volume, 2), 1)
                
                mean(min(cellCycleLength, predStruct.singleCell.repInitTime))
                mean(min(cellCycleLength, predStruct.singleCell.repTermTime))
                mean(min(cellCycleLength, predStruct.singleCell.cellCycleLen))
                
                predStruct.metConcs.mean
                predStruct.dnaSeq.mean
                predStruct.rnaSeq.mean
                predStruct.chipSeq.mean(:)
                predStruct.rnaArray.mean
                predStruct.protArray.mean
                predStruct.rxnFluxes.mean
                ];

predVec.std = [
                sqrt(mean(nanvar(predStruct.singleCell.growth, 1, 2)))
                sqrt(mean(nanvar(predStruct.singleCell.mass, 1, 2)))
                sqrt(mean(nanvar(predStruct.singleCell.volume, 1, 2)))
                
                std(min(cellCycleLength, predStruct.singleCell.repInitTime))
                std(min(cellCycleLength, predStruct.singleCell.repTermTime))
                std(min(cellCycleLength, predStruct.singleCell.cellCycleLen))
                
                predStruct.metConcs.std
                predStruct.dnaSeq.std
                predStruct.rnaSeq.std
                predStruct.chipSeq.std(:)
                predStruct.rnaArray.std
                predStruct.protArray.std
                predStruct.rxnFluxes.std
                ];

% % SingleCell
% figure;
% 
% subplot(2,2,1); hold on
% y1 = est1.singleCell.volume'; y1(y1 == 0) = nan;
% y2 = targ.singleCell.volume'; y2(y2 == 0) = nan;
% plot(y2(:,1),'-b'); plot(y1(:,1),'-r','LineWidth',2)
% plot(y2(:,2:end),'-b'); plot(y1(:,2:end),'-r','LineWidth',2)
% legend('target','estimate')
% ylabel('volume')
% 
% subplot(2,2,2); hold on
% y1 = est1.singleCell.mass'; y1(y1 == 0) = nan;
% y2 = targ.singleCell.mass'; y2(y2 == 0) = nan;
% plot(y2(:,1),'-b'); plot(y1(:,1),'-r','LineWidth',2)
% plot(y2(:,2:end),'-b'); plot(y1(:,2:end),'-r','LineWidth',2)
% legend('target','estimate')
% ylabel('mass')
% 
% subplot(2,2,3); hold on
% y1 = nanmean(est1.singleCell.volume)'; y1(y1 == 0) = nan;
% y2 = nanmean(targ.singleCell.volume)'; y2(y2 == 0) = nan;
% plot(y2,'-b','LineWidth',2); plot(y1,'-r','LineWidth',2)
% plot(y2+nanstd(targ.singleCell.volume)','--b');
% plot(y2-nanstd(targ.singleCell.volume)','--b');
% legend('target','estimate')
% ylabel('volume')
% title('Average')
% 
% subplot(2,2,4); hold on
% y1 = nanmean(est1.singleCell.mass)'; y1(y1 == 0) = nan;
% y2 = nanmean(targ.singleCell.mass)'; y2(y2 == 0) = nan;
% plot(y2,'-b','LineWidth',2); plot(y1,'-r','LineWidth',2)
% plot(y2+nanstd(targ.singleCell.mass)','--b');
% plot(y2-nanstd(targ.singleCell.mass)','--b');
% legend('target','estimate')
% ylabel('mass')
% title('Average')
% 
% 
% % protArray
% figure;
% hold on
% y_est = est1.protArray.mean; est_std = est1.protArray.std;
% y_targ = targ.protArray.mean; targ_std = targ.protArray.std;
% plot(y_targ,'-b','LineWidth',2)
% plot(y_targ+targ_std,'--b'); plot(y_targ-targ_std,'--b')
% plot(y_est,':r','LineWidth',2)
% plot(est2.protArray.mean,'-g')
% %plot(y_est+est_std,'--r'); plot(y_est-est_std,'--r')
% 
% predDist = alex_calcPredDist(targ,est1)
% predDist = alex_calcPredDist(targ,est2)
