function ok = test()

% Test rank 0, 1, and 2 ISTOs for spin-1 against explicit matrices

S = 1;

thr = 1e-10;

Op{1} = isto(S,[0 0]);
Op{2} = isto(S,[1 1]);
Op{3} = isto(S,[1 0]);
Op{4} = isto(S,[1 -1]);
Op{5} = isto(S,[2 2]);
Op{6} = isto(S,[2 1]);
Op{7} = isto(S,[2 0]);
Op{8} = isto(S,[2 -1]);
Op{9} = isto(S,[2 -2]);

% Reference
Op0{1} = eye(3);
Op0{2} = [0 -1 0; 0 0 -1; 0 0 0];
Op0{3} = diag([1 0 -1]);
Op0{4} = [0 0 0; 1 0 0; 0 1 0];
Op0{5} = [0 0 1; 0 0 0; 0 0 0];
Op0{6} = [0 -1 0; 0 0 1; 0 0 0]/sqrt(2);
Op0{7} = diag([sqrt(1/6) -sqrt(2/3) sqrt(1/6)]);
Op0{8} = [0 0 0; 1 0 0; 0 -1 0]/sqrt(2);
Op0{9} = [0 0 0; 0 0 0; 1 0 0];

for k = 1:numel(Op0)
  ok(k) = areequal(Op{k},Op0{k},thr,'abs');
end
