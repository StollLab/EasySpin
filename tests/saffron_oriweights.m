function ok = test()

% The orientation weights add up to 1: signals are averaged over the powder,
% and over the sample orientations and sites of a crystal

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz
Sys.AFrame = [0 0.6 0];

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',2);
Exp.ExciteWidth = Inf;
Exp.Sequence = '2pESEEM';
Exp.tau = 0.1;  % µs
Exp.dt = 0.02;  % µs
Exp.nPoints = 4;

Opt.TimeDomain = true;  % exact time-domain signal instead of histogram

% Powder without modulation: echo 1/2
Sys0 = Sys;
Sys0.A = [0 0];
[~,y] = saffron(Sys0,Exp,Opt);
ok(1) = areequal(y,1/2*ones(1,4),1e-10,'rel');

% Crystal: a repeated sample orientation gives the same signal
Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 pi/4 0];
[~,y1] = saffron(Sys,Exp,Opt);
Exp.SampleFrame = [0 pi/4 0; 0 pi/4 0];
[~,y2] = saffron(Sys,Exp,Opt);
ok(2) = areequal(y2,y1,1e-10,'rel');

% Crystal: average over symmetry-related sites
Exp.SampleFrame = [0 pi/4 0];
Exp.CrystalSymmetry = 'C2h';
[~,yAll] = saffron(Sys,Exp,Opt);
Opt.Sites = 1;
[~,ySite1] = saffron(Sys,Exp,Opt);
Opt.Sites = 2;
[~,ySite2] = saffron(Sys,Exp,Opt);
ok(3) = areequal(yAll,(ySite1+ySite2)/2,1e-10,'rel');

% Single crystal orientation with molecular frame = lab frame: same as
% spidyan (thyme mode, echo pathway selected)
p90 = struct('Flip',pi/2,'Phase',pi/2,'tp',0.016);
p180 = struct('Flip',pi,'Phase',pi/2,'tp',0.032);
ExpT = rmfield(Exp,{'ExciteWidth','Sequence','tau','dt','CrystalSymmetry'});
ExpT.SampleFrame = [0 0 0];
ExpT.Sequence = {p90 0.1 p180 0.1};
ExpT.Dim1 = {'d1,d2' 0.02};
ExpT.nPoints = 3;
ExpT.DetWindow = 0;
ExpT.PhaseCycle = {[0 1; pi/2 1i; pi -1; 3*pi/2 -1i]};
OptT.SimulationMode = 'thyme';
[~,ySaffron] = saffron(Sys,ExpT,OptT);
[~,ySpidyan] = spidyan(Sys,rmfield(ExpT,{'MolFrame','SampleFrame'}));
ok(4) = areequal(ySaffron(:),ySpidyan(:),1e-6,'rel');
