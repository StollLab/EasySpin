function ok = test()

% All operators for S=1

a = sqrt(2)/2;
Op{1} = 1i*[0 -a 0; a 0 -a; 0 a 0];
Op{2} = [1 0 0; 0 0 0; 0 0 -1];
Op{3} = [0 a 0; a 0 a; 0 a 0];

q = -1:1;
for iq = 1:numel(q)
  Opp = stev(1,1,q(iq));
  ok(iq) = areequal(Opp,Op{iq},1e-10,'abs');
end
