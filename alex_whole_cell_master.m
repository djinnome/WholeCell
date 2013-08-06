% alex_whole_cell_master: This script loads the "master table" that holds
%                         all of the previous simulation results, and adds
%                         new simulation data, that haven't been
%                         incorporated into the table yet.
%
% REMARKS:
%  -- This script calculates the "prediction distance" for each simulation.
%  -- The whole-cell-table is sorted by the "prediction distance"
%  -- This script is called by alex_runWholeCell_inf after a simulation
%       finishes.

load whole-cell-table.mat
w = what('./output');

for a = 1:length(w.mat)
  if ~any(strcmp(whole_cell_table.filename,w.mat{a}))
    % Calculate prediction distance
    est = load(['./output/' w.mat{a}]);
    if isfield(est,'perturbVec')
      tru = load('./mutant_WT_avgData');
      predDist = alex_calcPredDist(tru,est);
      est.predDist = predDist;
      save(['./output/' w.mat{a}],'predDist','-append');
      % Incorporate new information into the table
      whole_cell_table.perturbVecs = [whole_cell_table.perturbVecs est.perturbVec];
      whole_cell_table.predDist = [whole_cell_table.predDist; est.predDist];
      whole_cell_table.number = [whole_cell_table.number; length(whole_cell_table.number)+1];
      whole_cell_table.filename{end+1,1} = w.mat{a};
      % Sort everything based on whole_cell_table.predDist
      [~,I] = sort(whole_cell_table.predDist);
      whole_cell_table.predDist = whole_cell_table.predDist(I);
      whole_cell_table.perturbVecs = whole_cell_table.perturbVecs(:,I);
      whole_cell_table.number = whole_cell_table.number(I);
      whole_cell_table.filename = whole_cell_table.filename(I);
    end
  end
end
save('whole-cell-table.mat','whole_cell_table');
