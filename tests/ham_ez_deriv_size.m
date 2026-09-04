function ok = ham_ez_deriv_size()

% check the output with a simple case
Sys(1).S = 1/2;
Sys(1).g = 2*[1 1 1];

% check the output with a two electrons and a relative orientation between g tensors
Sys(2).S = [1/2 1/2];
Sys(2).g = [2*[1 1.9 1.1]; 3*[1 1 2]];
Sys(2).gFrame = [23 96 30 ; 10 6 81]*pi/180;
Sys(2).J=0;

for k = 1:numel(Sys)
  [muMx,muMy,muMz,dmuMx,dmuMy,dmuMz] = ham_ez(Sys(k));
  elSpins = numel(Sys(k).S);
  ok(k,1) = iscell(dmuMx);
  ok(k,2) = iscell(dmuMy);
  ok(k,3) = iscell(dmuMz);
  ok(k,4) = all(size(dmuMx)==[1 elSpins]);
  ok(k,5) = all(size(dmuMy)==[1 elSpins]);
  ok(k,6) = all(size(dmuMz)==[1 elSpins]);
  nElements = cellfun(@(x)numel(x),dmuMx);
  ok(k,7) = all(sum(nElements(:)) == elSpins*3);
  nElements = cellfun(@(x)numel(x),dmuMy);
  ok(k,8) = all(sum(nElements(:)) == elSpins*3);
  nElements = cellfun(@(x)numel(x),dmuMz);
  ok(k,9) = all(sum(nElements(:)) == elSpins*3);
end
