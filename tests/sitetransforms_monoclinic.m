function ok = test()

% Determine whether correct unique axes are assigned with various
% point and space group symbols for monoclinic group 3

ID = {3, 'C2', 'P2', 'P112', 'P121', 'P211'};

for k = 1:numel(ID)
  a = runprivate('sitetransforms',ID{k});
  result(k) = find(diag(a{2})==1);
end

ok = result == [2 2 2 3 2 1];
