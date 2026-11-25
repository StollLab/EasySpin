function ok = test() %#ok<FNDEF>

% The manifold engine simulates ENDOR as population transfer by the RF
% pulse, so it needs the electron spin in a population state during the RF
% pulse. Without a pathway selection, pathways with electron coherence
% during the RF pulse are included, and saffron stops with an error. The
% echo pathways exclude them, since the spins evolve during the RF pulse.

Sys.Nucs = '1H';
Sys.A = [1 1 3];  % MHz

p90.Flip = pi/2;
rf.Channel = 'rf';
rf.tp = 5;  % µs
rf.Flip = pi;

Exp.Field = 350;  % mT
Exp.MolFrame = [0 0.3 0];
Exp.CrystalSymmetry = 1;
Exp.rfFreq = 12;  % MHz
Exp.Sequence = {p90 0.2 p90 0 rf 0 p90 0.2};
Exp.nPoints = 41;
Exp.Dim1 = {'rfFreq',0.1};

% all pathways (default): error
try
  saffron(Sys,Exp);
  ok(1) = false;
catch ME
  ok(1) = contains(ME.message,'coherence order 0 during the RF pulse');
end

% population pathways during the RF pulse: manifold engine
Exp.Pathways = '+00-';
[~,y,info] = saffron(Sys,Exp);
ok(2) = strcmp(info.Method,'manifold') && any(abs(y-y(1))>0);

% echo pathways: the same as '+00-' (no delays around the RF pulse, so
% only its length prevents refocusing of coherence during it)
Exp.Pathways = 'echo';
[~,yEcho,info] = saffron(Sys,Exp);
ok(3) = strcmp(info.Method,'manifold') && areequal(yEcho,y,1e-12,'rel');
