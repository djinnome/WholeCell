function perturbVec = alex_randPerturb()
% alex_randPerturb: Generates a random 30x1 perturbation vector that obeys
%                   the constraints outlined in the competition
%                   description.
%
% CONSTRAINTS:
%    -- 5 genes have 1 associated parameter whose value was modified
%    -- 5 genes have 2 associated parameters whose value was modified
%    -- 13 out of the 15 parameters were decreases
%    -- Maximum positive perturbation is 1.906
%    -- Minimum positive perturbation is 1.117
%    -- Maximum negative perturbation is 0.066
%    -- Minimum negative perturbation is 0.972
%
% Author: Alex Williams



perturbVec = ones(30,1);

% Generate the "stash" of possible perturbations
stash(1) = (1-0.028); % 2.8 percent decrease
stash(2) = (1-0.934); % 93.4 percent decrease
stash(3) = (1+0.117); % 11.7 percent increase
stash(4) = (1+0.906); % 90.6 percent increase

% The rest of the perturbations are decreases between 2.8 and 93.4
stash(5:15) = 1 - (0.028+rand(11,1).*(0.934-0.028));

% Sprinkle the perturbations in a random order
stash = stash(randperm(15));
double_pert_genes = logical([ones(5,1); zeros(5,1)]);
double_pert_genes = double_pert_genes(randperm(10));

% For the first five genes, perturb 2 of the 3 parameters
i_pert = 1;
i_stash = 1;
for i_gene = 1:10
  if double_pert_genes(i_gene)
    i_param = randsample((0:2)',2);
    perturbVec(i_pert+i_param(1)) = stash(i_stash);
    perturbVec(i_pert+i_param(2)) = stash(i_stash+1);
    i_stash = i_stash+2;
  else
    perturbVec(i_pert+(randi(3)-1)) = stash(i_stash);
    i_stash = i_stash+1;
  end
  i_pert = i_pert+3;
end