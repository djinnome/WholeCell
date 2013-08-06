% This script summarizes an analysis of how promoter affinity and RNA half
% life affect the data of the rnaSeq and rnaArray data types

% Description of output:
%
%  Figures 1-10:
%    Each figure shows the effect of manipulating parameters of a single
%    gene. Thus, there are ten figures, one for each gene. Each panel on a
%    figure shows the effect of manipulating each parameter for that gene
%    (top panel = promotor affinity, middle panel = RNA half life, and
%    bottom panel = reaction kcat). The blue lines show the ABSOLUTE change
%    in rnaArray.mean when the specified parameter is increased. The red
%    lines show the ABSOLUTE change in rnaArray.mean when the specified
%    parameter is decreased. The green dots show the location of the gene
%    of interest.
%
%  Figures 11-20:
%    Same as the first ten figures, but each panel shows the PERCENT change
%    in rnaArray.mean when the specified parameter is increased/decreased.

% Interpretation:
%
%   Looking at the percent change in mRNA expression seems to be slightly
%   more informative. However, the rna array data in general looks pretty
%   complicated and non-linear. There doesn't appear to be an obvious way
%   of gleaning parameter information from these plots. We may need to try
%   something more complicated (e.g. a PCA analysis to reduce the
%   dimensionality of this dataset, see s001)

close all; clear; clc
set(0,'DefaultFigureWindowStyle','docked')

cd ..
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
cd analysis

% Labels
genesTusRxns = {
    'MG_006' 'TU_003' 'Tmk'
    'MG_023' 'TU_011' 'Fba'
    'MG_047' 'TU_027' 'MetK'
    'MG_111' 'TU_069' 'Pgi'
    'MG_272' 'TU_180' 'AceE'
    'MG_299' 'TU_203' 'Pta'
    'MG_330' 'TU_233' 'CmkA2'
    'MG_357' 'TU_260' 'AckA'
    'MG_407' 'TU_294' 'Eno'
    'MG_431' 'TU_307' 'TpiA'
    };
parameterTypes = {
    'PromAffinity'
    'HalfLife'
    'RxnKcat'
    };
gene_ids = sim.gene.wholeCellModelIDs;
  
% Specify path to averaged-output
basedir = '/Users/zucker/src/WholeCell/';
outputdir = [basedir '/output/'];
filename = {};

% Unmodified model (original params) and mutant data
unmod = load([outputdir 'averaged_sim-0.mat']);
mutant = load([basedir 'mutant_WT_avgData.mat']);

% First 10 figures:
a = 1; % index from 1-60 for each perturbation
for g = 1:10
  figure;
  set(gcf,'Name','Change in RNA expression')
  for p = 1:3
    subplot(3,1,p); hold on
    down = load([outputdir 'averaged_sim-' num2str(a)]);
    up = load([outputdir 'averaged_sim-' num2str(a+1)]);
    gene_index = strcmp(genesTusRxns{g , 1},gene_ids);
    up_profile = up.rnaArray.mean - unmod.rnaArray.mean;
    down_profile = down.rnaArray.mean - unmod.rnaArray.mean;
    x = 1:length(up_profile);
    plot(x,down_profile,'-r')
    plot(x,up_profile,'-b')
    plot(x(gene_index),up_profile(gene_index),'.g','MarkerSize',15)
    plot(x(gene_index),down_profile(gene_index),'.g','MarkerSize',15)
    
    xlim([x(1) x(end)])
    a = a+2;
    if p==1
      title(['Promoter Affinity : ' genesTusRxns{g,1} '  //  ' genesTusRxns{g,2} '  //  ' genesTusRxns{g,3}],'interpreter','none')
    end
    if p==2
      title('mRNA Half Life')
      ylabel({'Change in Concentration','from original/unmodified model'})
    end
    if p==3
      title('Reaction kcat')
      legend('0.5x','2.0x')
      xlabel('Different mRNAs')
    end
  end
  export_fig(['s000_' num2str(gcf)],'-pdf')
  filename{end+1} = ['s000_' num2str(gcf) '.pdf'];
end

% Second 10 figures:
a = 1; % index from 1-60 for each perturbation
for g = 1:10
  figure;
  set(gcf,'Name','Percent change in RNA expression')
  for p = 1:3
    subplot(3,1,p); hold on
    down = load([outputdir 'averaged_sim-' num2str(a)]);
    up = load([outputdir 'averaged_sim-' num2str(a+1)]);
    gene_index = strcmp(genesTusRxns{g , 1},gene_ids);
    up_profile = 100.*(up.rnaArray.mean - unmod.rnaArray.mean)./unmod.rnaArray.mean;
    down_profile = 100.*(down.rnaArray.mean - unmod.rnaArray.mean)./unmod.rnaArray.mean;
    x = 1:length(up_profile);
    plot(x,down_profile,'-r')
    plot(x,up_profile,'-b')
    plot(x(gene_index),up_profile(gene_index),'.g','MarkerSize',15)
    plot(x(gene_index),down_profile(gene_index),'.g','MarkerSize',15)
    
    xlim([x(1) x(end)])
    a = a+2;
    if p==1
      title(['Promoter Affinity : ' genesTusRxns{g,1} '  //  ' genesTusRxns{g,2} '  //  ' genesTusRxns{g,3}],'interpreter','none')
    end
    if p==2
      title('mRNA Half Life')
      ylabel({'Percent change in Concentration','from original/unmodified model'})
    end
    if p==3
      title('Reaction kcat')
      legend('0.5x','2.0x')
      xlabel('Different mRNAs')
    end
  end
  export_fig(['s000_' num2str(gcf)],'-pdf')
  filename{end+1} = ['s000_' num2str(gcf) '.pdf'];
end

append_pdfs('s000_output.pdf',filename{:});

for a = 1:length(filename)
  delete(filename{a})
end

