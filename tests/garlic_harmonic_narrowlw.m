function ok = test()

% Check that a harmonic with a line width narrower than two increments raises
% an error, instead of differentiating the stick spectrum numerically

Sys.g = 2;
Exp.mwFreq = 9.5;
Exp.Range = [330 350];
Exp.nPoints = 1024;
Exp.Harmonic = 1;

dx = diff(Exp.Range)/(Exp.nPoints-1);
Sys.lw = dx/2;

try
  garlic(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'smaller than 2 increments');
end

end
