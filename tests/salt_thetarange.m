function ok = test(opt)

% ThetaRange option

Sy = struct('S',1/2,'Nucs','1H','g',2,...
  'A',2*[-1,-1,2],'lwEndor',0.05);
Ex.Range = [10 18];
Ex.Field = 326.5;
Op.nKnots = 40;
Op.Verbosity = opt.Verbosity;
Op.Method = 'perturb1';
th = linspace(0,pi/2,5);
for k = 1:numel(th)-1
  Op.ThetaRange = th(k:k+1);
  [x,y(k,:)] = salt(Sy,Ex,Op);
end

if (opt.Display)
  plot(x,y);
  title('Orientation selection via Opt.ThetaRange');
end

ok = true;
