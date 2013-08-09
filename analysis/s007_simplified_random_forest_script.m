% This script was an attempt to find any pairwise correlations between the
% principal components of the high-throughput data, and the perturbations
% that produced that data. This doesn't work very well.
%
% This script is still useful because it compiles all the highthroughput
% data into the f000_*.mat files, which contain the relevant features of
% each dataset in a simple matrix format. This is used by later scripts.

close all; clear; clc
set(0,'DefaultFigureWindowStyle','docked')

datatypes = {'metConcs'
             'dnaSeq'
             'rnaSeq'
             'chipSeq'
             'rnaArray'
             'protArray'
             'rxnFluxes'};

% Calculate principal components for each datatype and save them into the
% "compiled_data" directory
for a = 1:length(datatypes)
  f000_predictPerturbs_pairwiseCorr(datatypes{a})
end

% Run the script to build estimates for perturbVec using random forest. 
s006_improved_random_forest_est