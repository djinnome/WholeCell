function pV = alex_constrainPerturb(pV)
% alex_constrainPerturb: Constrains a given perturbation vector to meet the
%                        criteria listed in the DREAM challenge.
%
% INPUT:  "pV" - a 30x1 perturbation vector
%
% OUTPUT: "pV" - a 30x1 perturbation vector, modified such that 13 of the
%                entries are less than one, 2 of the entries are greater
%                than one, and 15 of the entries are equal to one. In
%                addition, none of the entries are greater than 1.906 and
%                none of the entries are less than 0.0660, and 5 genes have
%                1 non-zero perturbation and 5 other genes have 2 non-zero
%                perturbations. I tried to meet these specifications by
%                making the fewest number of manipulations.
%
% Author: Alex Williams

% pV is "perturbVec" - abbreviated for clarity in this function.
assert(length(pV) == 30);

% Establish bounds of [0.0660, 1.906]
pV(pV > 1.906) = 1.906;
pV(pV < 0.0660) = 0.0660;

NUM_NEG_PAR = 13;
NUM_POS_PAR = 2;
NUM_ONEPERT = 5; % Five genes have one perturbed param
NUM_TWOPERT = 5; % Five genes have two perturbed param

% set diffs to non-zero to start the while loop
negDiff = 1; posDiff = 1; oneDiff = 1; twoDiff = 1;
while negDiff ~= 0 || posDiff ~=0 || oneDiff ~= 0 || twoDiff ~= 0 
  
  % Ensure each gene has either one or two perturbations
  a = 1;
  for i_gene = 1:10
    if sum(pV(a:a+2) == 1) == 0
      % the perturbations to the first two parameters are somewhat
      % redundant, keep the higher one, erase the lower
      if abs(pV(a)-1) < abs(pV(a+1)-1)
        pV(a) = 1;
      else
        pV(a+1) = 1;
      end
    elseif sum(pV(a:a+2) == 1) == 3
      % add a random positive or negative perturbation
      if rand>0.5
        pV(a:a+(randi(3)-1)) = 0.9720; 
      else
        pV(a:a+(randi(3)-1)) = 1.1170;
      end
    end
    a = a+3;
  end
  
  % check number of pos/neg params and perturbations per gene
  negDiff = sum(pV<1) - NUM_NEG_PAR; % actual minus desired neg params
  posDiff = sum(pV>1) - NUM_POS_PAR; % actual minus desired pos params
  a = 1; s1 = 0; s2 = 0;
  for i_gene = 1:10
    if sum(pV(a:a+2) == 1) == 1
      s2 = s2 + 1; % this gene has 2 perturbations
    elseif sum(pV(a:a+2) == 1) == 2
      s1 = s1 + 1; % this gene has 1 perturbations
    end
    a = a+3;
  end
  oneDiff = s1 - NUM_ONEPERT; % actual minus desired genes with 1 perturb
  twoDiff = s2 - NUM_TWOPERT; % actual minus desired genes with 2 perturb
  
  % construct candidates to erase:
  eraseC = false(size(pV)); a = 1;
  for i_gene = 1:10
    if sum(pV(a:a+2) == 1) <= 1 % check if there are two perturbations
      eraseC(a:a+2) = (pV(a:a+2) ~= 1);
    end
    a = a+3;
  end
  if any(eraseC)
    for x = 1:negDiff
      min_pert = max(pV((pV < 1) & eraseC)-1);
      min_ind = find(pV == (min_pert+1));
      min_ind = min_ind(randi(length(min_ind)));
      eraseC(min_ind) = false;
      pV(min_ind) = 1; % erase the smallest candidate
    end
    for x = 1:posDiff
      min_pert = min(pV((pV > 1) & eraseC)-1);
      min_ind = find(pV == (min_pert+1));
      min_ind = min_ind(randi(length(min_ind)));
      eraseC(min_ind) = false;
      pV(min_ind) = 1; % erase the smallest candidate
    end
  end
  
  % construct candidates to add:
  addCands = false(size(pV)); % find genes with 0 or 1 perturbations
  a = 1;
  for i_gene = 1:10
    if sum(pV(a:a+2) == 1) >= 2 % check if there is 1 or more perturbations
      addCands(a:a+2) = (pV(a:a+2) == 1);
    end
    a = a+3;
  end
  addPerts = [];
  if negDiff < 0
    addPerts = [addPerts; ones(abs(negDiff),1)*0.9720]; %#ok<AGROW>
  end
  if posDiff < 0
    addPerts = [addPerts; ones(abs(posDiff),1)*1.1170]; %#ok<AGROW>
  end
  if ~isempty(addPerts) && sum(addCands)~=0
    len = min([length(addPerts) sum(addCands)]);
    pV(randsample(find(addCands),len)) = addPerts(1:len);
  end
end

% minimum postive and negative perturbations
minPositive = min(pV(pV>1));
pV(pV == minPositive) = 1.1170;
minNegative = max(pV(pV<1));
pV(pV == minNegative) = 0.9720;
