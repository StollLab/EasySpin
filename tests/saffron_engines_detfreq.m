function ok = test()  %#ok<FNDEF>

% Detection frame (Exp.DetFreq) in the manifold engine: the same signals as
% the full engine, with the phase from the time between the start of
% the sequence and the detection point
% - delays before the first pulse count
% - indirect dimension over delays (the frequency of the signal shifts)
% - indirect dimension over a pulse length (computed point by point)
% - two dimensions (HYSCORE)
% - DetFreq = mwFreq and empty DetFreq give the signal without detection frame

Sys.g = 2.005;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz, spins 2.5 MHz off resonance
Exp.MolFrame = [0.3 0.6 0.2];
Exp.CrystalSymmetry = 1;
delta = 0.003;  % GHz
Exp.DetFreq = Exp.mwFreq + delta;

man = struct('Method','manifold');
full = struct('Method','full','SimFreq',9.7);

% Two-pulse ESEEM with ideal pulses and a delay before the first pulse
p90.Flip = pi/2;
p180.Flip = pi;
E = Exp;
E.Sequence = {0.05 p90 0.2 p180 0.2};
E.nPoints = 16;
E.Dim1 = {'d2,d3' 0.01};  % µs
[t,yMan,info] = saffron(Sys,E,man);
yFull = saffron(Sys,E,full);
ok(1) = strcmp(info.Method,'manifold') && areequal(yMan,yFull,1e-4,'rel');

% The detection frame only multiplies by the phase factor
E0 = rmfield(E,'DetFreq');
y0 = saffron(Sys,E0,man);
tSeq = 0.05 + 2*t;  % µs, from the start of the sequence to detection (t: tau)
ok(2) = areequal(yMan,y0.*exp(-2i*pi*delta*1e3*tSeq),1e-12,'abs') && ...
  max(abs(yMan-y0))>0.1*max(abs(y0));

% Finite pulses, indirect dimension over the length of the first pulse
p90f = struct('Flip',pi/2,'tp',0.01);  % µs
p180f = struct('Flip',pi,'tp',0.02);  % µs
E = Exp;
E.Sequence = {p90f 0.2 p180f 0.2};
E.nPoints = 3;
E.Dim1 = {'p1.tp' 0.002};  % µs
[~,yMan,info] = saffron(Sys,E,man);
yFull = saffron(Sys,E,full);
ok(3) = strcmp(info.Method,'manifold') && areequal(yMan,yFull,1e-3,'rel');

% Two dimensions (HYSCORE): the phase changes along both
E = Exp;
E.Sequence = {p90 0.12 p90 0 p180 0 p90 0.12};
E.nPoints = [6 5];
E.Dim1 = {'d2' 0.02};  % µs
E.Dim2 = {'d3' 0.03};  % µs
E.Pathways = '+00-';
yMan = saffron(Sys,E,man);
yFull = saffron(Sys,E,full);
y0 = saffron(Sys,rmfield(E,'DetFreq'),man);
[t1,t2] = ndgrid(0.02*(0:5),0.03*(0:4));
tSeq = 0.24 + t1 + t2;  % µs
ok(end+1) = isequal(size(yMan),[6 5]) && areequal(yMan,yFull,1e-4,'rel') && ...
  areequal(yMan,y0.*exp(-2i*pi*delta*1e3*tSeq),1e-12,'abs');

% DetFreq = mwFreq: no detection frame
E = rmfield(E,'DetFreq');
y0 = saffron(Sys,E,man);
E.DetFreq = Exp.mwFreq;
y1 = saffron(Sys,E,man);
ok(end+1) = isequal(y1,y0);

% Empty DetFreq: no detection frame, in both engines
E.DetFreq = [];
ok(end+1) = isequal(saffron(Sys,E,man),y0) && ...
  areequal(saffron(Sys,E,full),saffron(Sys,rmfield(E,'DetFreq'),full),1e-12,'abs');

end
