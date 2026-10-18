function ok = test()

% Check that a nonzero harmonic without line width raises an error
% instead of silently returning the absorption spectrum

Sys.Nucs = '1H';
Sys.A = [2 3 4];
Exp.Field = 350;
Exp.Range = [10 20];
Exp.Harmonic = 1;

try
  salt(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'No broadening given');
end

end
