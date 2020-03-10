function ok = test()

% Test cases where wigner3j implements explicit dedicated expressions

J = 17;
M = -3;

% All cases with J3 = 2, M3>=0, J1>=J2
jm{1} = [J+2, J, 2; M, -M-2, 2];
jm{2} = [J+2, J, 2; M, -M-1, 1];
jm{3} = [J+2, J, 2; M, -M  , 0];
jm{4} = [J+1, J, 2; M, -M-2, 2];
jm{5} = [J+1, J, 2; M, -M-1, 1];
jm{6} = [J+1, J, 2; M, -M  , 0];
jm{7} = [J  , J, 2; M, -M-2, 2,];
jm{8} = [J  , J, 2; M, -M-1, 1];
jm{9} = [J  , J, 2; M, -M  , 0];

% All cases with J3 = 1, M3>=0, J1>=J2
jm{10} = [J+1, J, 1; M, -M-1, 1];
jm{11} = [J+1, J, 1; M, -M  , 0];
jm{12} = [J  , J, 1; M, -M-1, 1];
jm{13} = [J  , J, 1; M, -M  , 0];

% All cases with J3 = 0
jm{14} = [J  , J, 0; M, -M  , 0];

% All cases with J3 = 1/2, M3>=0, J1>=J2
M = M+1/2;
jm{15} = [J+1/2, J, 1/2; M, -M-1/2, 1/2];

% All cases with J3 = 1/2, M3>=0, J1>=J2
jm{16} = [J+3/2, J, 3/2; M, -M-3/2, 3/2];
jm{17} = [J+3/2, J, 3/2; M, -M-1/2, 1/2];
jm{18} = [J+1/2, J, 3/2; M, -M-3/2, 3/2];
jm{19} = [J+1/2, J, 3/2; M, -M-1/2, 1/2];

nBaseCases = numel(jm);

% All cases, reordered all possible ways
for k = 1:nBaseCases
  jm{end+1} = jm{k}(:,[1 3 2]);
  jm{end+1} = jm{k}(:,[2 1 3]);
  jm{end+1} = jm{k}(:,[2 3 1]);
  jm{end+1} = jm{k}(:,[3 1 2]);
end

% All cases, inverted m
for k = 1:nBaseCases
  jm{end+1} = jm{k}.*[1;-1];
end

for k = 1:numel(jm)
  a0 = wigner3j(jm{k},'f');
  a1 = wigner3j(jm{k},'f+');
  ok(k) = abs(a0-a1)<1e-10;
end
