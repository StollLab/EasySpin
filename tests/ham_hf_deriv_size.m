function ok = test()

Sys(1).S = 1/2;
Sys(1).Nucs = '1H';
Sys(1).A = [1 2 3];

Sys(2).S = 1/2;
Sys(2).Nucs = '1H,14N';
Sys(2).A = [1 2 3; 4 5 6];

Sys(3).S = [1/2 1];
Sys(3).Nucs = '1H';
Sys(3).A = [1 2 3, 4 5 6];
Sys(3).J = 100;

Sys(4).S = [1/2 1];
Sys(4).Nucs = '1H,13C';
Sys(4).A = [1 2 3, 4 5 6; 7 8 9, 11 15 19];
Sys(4).J = 100;

for k = 1:numel(Sys)
  [~,dHhf] = ham_hf(Sys(k));
  elSpins = numel(Sys(k).S);
  I = nucspin(Sys(k).Nucs);
  nucSpins = numel(I);
  ok(k,1) = iscell(dHhf);
  ok(k,2) = all(size(dHhf)==[elSpins nucSpins]);
  nElements = cellfun(@(x)numel(x),dHhf);
  ok(k,3) = all(nElements(:)==3);
end
