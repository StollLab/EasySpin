function ok = ham_zf_deriv_size()

Sys(1).S = [1];
Sys(1).D = 100*[-1 -1 2];

Sys(2).S = [3/2 5/2];
D = [1 2]';
E = [0.5 0.3]';
Sys(2).D = D*[-1,-1,2]/3 + E*[1,-1,0];
Sys(2).DFrame = [23 96 30 ; 10 6 81]*pi/180;
Sys(2).J=0;

for k = 1:numel(Sys)
  [~,dHzf] = ham_zf(Sys(k));
  elSpins = numel(Sys(k).S);
  ok(k,1) = iscell(dHzf);
  ok(k,2) = all(size(dHzf)==[1 elSpins]);
  nElements = cellfun(@(x)numel(x),dHzf);
  ok(k,3) = all(sum(nElements(:)) == elSpins*3);
end
