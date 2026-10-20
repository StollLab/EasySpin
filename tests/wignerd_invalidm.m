function ok = test()

% Invalid m1 or m2 must throw an error that names the offending index

beta = 0.7;
invalid = {[1 0.5 0], [1 2 0], [5 0.5 0], [11/2 1 1/2], ...
           [1 0 0.5], [1 0 -2], [5 0 1.5], [11/2 1/2 1]};
name = {'m1','m1','m1','m1','m2','m2','m2','m2'};
for k = numel(invalid):-1:1
  try
    wignerd(invalid{k},beta);
    ok(k) = false;
  catch err
    ok(k) = startsWith(err.message,name{k});
  end
end

% Valid half-integer indices must not throw
wignerd([11/2 -3/2 1/2],beta);
ok(end+1) = true;
