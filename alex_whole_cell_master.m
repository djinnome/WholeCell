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

load averaged-whole-cell-table.mat
w = what('./averaged_output');

for a = 1:length(w.mat)
  if ~any(strcmp(avg_cell_table.filename,w.mat{a}))
    % Calculate prediction distance
    est = load(['./averaged_output/' w.mat{a}]);
    if isfield(est,'perturbVec')
      tru = load('./mutant_WT_avgData');
      predDist = alex_calcPredDist(tru,est);
      est.predDist = predDist;
      save(['./averaged_output/' w.mat{a}],'predDist','-append');
      % Incorporate new information into the table
      avg_cell_table.perturbVecs = [avg_cell_table.perturbVecs est.perturbVec];
      avg_cell_table.predDist = [avg_cell_table.predDist; est.predDist];
      avg_cell_table.number = [avg_cell_table.number; length(avg_cell_table.number)+1];
      avg_cell_table.filename{end+1,1} = w.mat{a};
      % Sort everything based on whole_cell_table.predDist
      [~,I] = sort(avg_cell_table.predDist);
      avg_cell_table.predDist = avg_cell_table.predDist(I);
      avg_cell_table.perturbVecs = avg_cell_table.perturbVecs(:,I);
      avg_cell_table.number = avg_cell_table.number(I);
      avg_cell_table.filename = avg_cell_table.filename(I);
    end
  end
end
save('averaged-whole-cell-table.mat','avg_cell_table');
