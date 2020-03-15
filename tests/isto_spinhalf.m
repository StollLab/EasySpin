function ok = test()

% Test rank 0 and 1 ISTOs for spin-1/2 against explicit matrices

S = 1/2;

thr = 1e-10;

Op{1} = isto(S,[0 0]);
Op{2} = isto(S,[1 1]);
Op{3} = isto(S,[1 0]);
Op{4} = isto(S,[1 -1]);

% Reference
Op0{1} = eye(2);
Op0{2} = [0 -1; 0 0]/sqrt(2);
Op0{3} = diag([1/2 -1/2]);
Op0{4} = [0 0; 1 0]/sqrt(2);

for k = 1:numel(Op0)
  ok(k) = areequal(Op{k},Op0{k},thr,'abs');
end
