% alex_whole_cell_master: This script continuously runs whole cell
%                         simulations. Parameters can be fed in through the
%                         "good_guesses" mat file. If that is empty, this
%                         script uses differential crossover to generate a
%                         new set of parameters.
%
% REMARKS:
%  -- The parameters of the differential crossover can be modified by
%       changing the mat files "fitness_quantile_cutoff" and
%       "chance_to_constrain_perturbs"
%  -- fitness_quantile_cutoff specifies which models can be picked by the
%       differential crossover algorithm. I started with the top 50% of
%       models. This could be decreased to increase the "greediness" of the
%       algorthm, or increased to search a broader space.
%  -- chance_to_constrain_perturbs will be important to modify towards the
%       end of the competition when we are trying to decrease parameter
%       distance.

while true
  rng('shuffle')
  % Create perturb vector
  load('good_guesses.mat')
  if ~isempty(guesses)
    % if we have a guess, use that as a starting point. Hope for no race
    % condition. (Not the end of the world if two identical simulations).
    perturbVec = guesses(:,1);
    guesses = guesses(:,2:end);
    save('good_guesses.mat','guesses')
  else
    load('whole-cell-table.mat')
    popSize = size(whole_cell_table.perturbVecs,2);
    if popSize < 100
      perturbVec = alex_randPerturb;
    else
      load('fitness_quantile_cutoff.mat')
      % sample from top 50% of population (can be changed by altering the
      % mat file 'fitness_quantile_cutoff'
      N = round(popSize*fitness_quantile_cutoff);
      passed = whole_cell_table.perturbVecs(:,1:N); % perturbVecs is sorted
      p_ind = randsample(N,3); % Three parents
      p(:,1) = whole_cell_table.perturbVecs(:,min(p_ind)); % best parent
      p(:,2) = whole_cell_table.perturbVecs(:,median(p_ind)); % mid parent
      p(:,3) = whole_cell_table.perturbVecs(:,max(p_ind)); % worst parent
      couple = sort(randsample(3,2));
      mutation_vec = p(:,couple(1))-p(:,couple(2)); % better minus worse
      child = find((1:3 ~= couple(1)) & (1:3 ~= couple(2)));
      perturbVec = p(:,child) + (rand*0.8+0.2).*mutation_vec;
      
      % Do we wish to constrain the perturbation to meet the specifications
      % of the contest?
      load('chance_to_constrain_perturbs.mat') % on unit interval
      perturbVec(perturbVec > 1.906) = 1.906;
      perturbVec(perturbVec < 0.0660) = 0.0660;
      r = rand;
      if r < chance_to_hard_constrain
        perturbVec = alex_constrainPerturb(perturbVec);
      elseif r < (chance_to_loose_constrain + chance_to_hard_constrain)
        % loose constraint - only change 17 params, delete positive
        % perturbations first
        while sum(perturbVec ~= 1) < 17
          if sum(perturbVec > 3)
            minPositive = min(perturbVec(perturbVec>1));
            perturbVec(perturbVec == minPositive) = 1;
          else
            minNegative = max(perturbVec(perturbVec<1));
            perturbVec(perturbVec == minNegative) = 1;
          end
        end
      end
    end
  end
  
  % Run the simulation (results are saved to './output' directory)
  alex_perturb_wholeCell(perturbVec);
  
  % run the "master" script, to update whole_cell_table
  clear;
  alex_whole_cell_master
  
  % Clear everything and do it again.
  clear;
end