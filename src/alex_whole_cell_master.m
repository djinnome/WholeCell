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
