function ok = test()

% Check that invalid detection harmonics raise errors

Sys.g = 2;
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.Range = [330 350];

harmonics = {-1, 3, 0.5, [1 2]};

for k = numel(harmonics):-1:1
  Exp.Harmonic = harmonics{k};
  try
    pepper(Sys,Exp);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,'Exp.Harmonic must be 0, 1 or 2');
  end
end

end
