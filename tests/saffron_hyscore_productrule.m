function ok = test() %#ok<FNDEF>

% HYSCORE of two 1H nuclei in a single orientation gives the same signal
% with the product rule (nuclei simulated separately, signals multiplied)
% and without it (both nuclei simulated together). Exp.nPoints is small
% enough (8x8 points, at most 64) that both calculations evaluate each point
% exactly in the time domain, without frequency binning errors.

Sys = struct('Nucs','1H,1H','A',[1 1 3; 4 4 2]);  % A in MHz

% Field in mT, dt and tau in µs
Exp = struct('Field', 1213.2,...
  'Sequence' , 'HYSCORE',...
  'dt' , 0.004,...
  'nPoints',8,...
  'tau' ,0.1);

% Single orientation, with the field tilted away from the tensor z axes
Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 pi/4 0];

Opt.ProductRule = 0;
[~,y0] = saffron(Sys,Exp,Opt);
Opt.ProductRule = 1;
[~,y1] = saffron(Sys,Exp,Opt);

ok = areequal(y1,y0,1e-10,'rel');
