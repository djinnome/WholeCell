% This is a script that I ran to generate some initial guesses for the
% population of random models. I based these guesses on the simple
% experiments where each parameter was perturbed by 0.5x or 2.0x. This
% isn't really that important anymore.
guesses(:,1) = alex_randPerturb;

while size(guesses,2) < 5
  rp = alex_randPerturb;
  if sum(rp(1:3)<1)==2 && sum(rp([4 14 19])>1)==2 && rp(23)>0
    guesses(:,end+1) = rp;
  end
end
while size(guesses,2) < 15
  rp = alex_randPerturb;
  if sum(rp(1:3)<1)==2 && sum(rp([4 14 19])>1)==2 
    guesses(:,end+1) = rp;
  end
end
while size(guesses,2) < 45
  rp = alex_randPerturb;
  if ~any(rp(1:3)<1) && sum(rp([4 14 19])>1) >= 1 
    guesses(:,end+1) = rp;
  end
end