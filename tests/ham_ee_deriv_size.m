function ok = ham_ee_deriv_size()

% simple spin system
Sys(1).S = [1/2 1/2];
Sys(1).ee = [1 2 3];
Sys(1).ee2 = 100;

% more complex one
Sys(2).S = [1/2 1 1];
Sys(2).ee = [1 2 3; 4 5 6 ; 7 8 9];
Sys(2).ee2 = [100 200 300];

for k = 1:numel(Sys)
  [~,dHee] = ham_ee(Sys(k));
  elSpins = numel(Sys(k).S);
  ok(k,1) = iscell(dHee);  %are these cell?
  ok(k,2) = all(size(dHee)==[1 elSpins*(elSpins-1)/2]);  % what about their size
  nElements = cellfun(@(x)numel(x),dHee);
  ok(k,3) = all(sum(nElements(:)) == elSpins*(elSpins-1)/2*4); %ensure the size are correct
end
