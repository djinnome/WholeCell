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

% delete the output pdf so we can make it again.
delete('s002_output.pdf') 
pdfnames = {};

% Calculate principal components for each datatype
for a = 1:length(datatypes)
  f000_predictPerturbs_pairwiseCorr(datatypes{a})
  export_fig(['s002_' num2str(gcf)],'-pdf')
  pdfnames{end+1} = ['s002_' num2str(gcf) '.pdf'];
end

% Make the output pdf
append_pdfs('s002_output.pdf',pdfnames{:});
for a = 1:length(pdfnames)
  delete(pdfnames{a})
end