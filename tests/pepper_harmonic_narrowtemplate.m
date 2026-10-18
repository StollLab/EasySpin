function ok = test()

% Check that line widths narrower than two increments give the same error for
% isotropic spectra (line templates) as for powder spectra (convolution)

Exp.mwFreq = 9.5;
Exp.Range = [330 350];
Exp.nPoints = 1001;
Exp.Harmonic = 1;

dx = diff(Exp.Range)/(Exp.nPoints-1);

lw = {[dx/2 0], [0 dx/2], [dx/2 dx/2]};
for k = 1:numel(lw)
  Sys = struct('g',2,'lw',lw{k});
  try
    pepper(Sys,Exp);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,'smaller than 2 increments');
  end
end

% A wide line is fine
Sys = struct('g',2,'lw',[0 5*dx]);
try
  pepper(Sys,Exp);
  ok(end+1) = true;
catch
  ok(end+1) = false;
end

end
