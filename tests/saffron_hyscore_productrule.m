function ok = test()

% HYSCORE with and without product rule give the same signal

Sys = struct('Nucs','1H,1H','A',[1 1 3; 4 4 2]);

Exp = struct('Field', 1213.2,...
  'Sequence' , 'HYSCORE',...
  'dt' , 0.004,...
  'nPoints',10,...
  'tau' ,0.1);

Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 pi/4 0];

Opt.TimeDomain = true;  % exact time-domain signal instead of histogram
Opt.ProductRule = 0;
[~,y0] = saffron(Sys,Exp,Opt);
Opt.ProductRule = 1;
[~,y1] = saffron(Sys,Exp,Opt);

ok = areequal(y1,y0,1e-10,'rel');
