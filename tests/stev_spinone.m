function ok = test()

% All rank-1 operators for S=1

a = sqrt(2)/2;
Op{1} = 1i*[0 -a 0; a 0 -a; 0 a 0];
Op{2} = [1 0 0; 0 0 0; 0 0 -1];
Op{3} = [0 a 0; a 0 a; 0 a 0];

S = 1;
k = 1;
q = -k:k;
for iq = 1:numel(q)
  Opp = stev(S,[k,q(iq)]);
  ok(iq) = areequal(Opp,Op{iq},1e-10,'abs');
end
