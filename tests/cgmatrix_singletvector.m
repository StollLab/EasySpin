function ok = test()

% Check whether cgmatrix() correctly selects the coupled state with the
% specific total spin given in the third input argument.

S_list = [1/2 1 3/2];

for k = 1:numel(S_list)
  S = S_list(k);
  Ufull = cgmatrix(S,S);
  Usinglet = cgmatrix(S,S,0);
  ok(k) = areequal(Ufull(end,:),Usinglet,1e-10,'abs');
end
