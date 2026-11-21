function ok = test()  %#ok<FNDEF>

% Single-point detection shifted from the end of the sequence (scalar
% Exp.DetWindow) in the manifold engine: the same signals as the full
% engine
% - positive and negative shifts, ideal and rectangular pulses
% - a shift after a final pulse is the same as a final delay
% - with Exp.DetFreq, the phase includes the shift
% - the shift takes the relaxation setting of the last delay
%   (Opt.Relaxation)
% - echo times stay relative to the end of the sequence
% - transients and 'sequence' need the full engine
% - errors: invalid windows, detection points in a pulse

Sys.g = 2.005;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz, spins 2.5 MHz off resonance
Exp.MolFrame = [0.3 0.6 0.2];
Exp.CrystalSymmetry = 1;

man = struct('Method','manifold');
full = struct('Method','full','SimFreq',9.7);
ok = false(1,0);

% Two-pulse ESEEM with ideal pulses, shifted detection point: the
% off-resonance spins dephase, so the shift changes the signal
p90.Flip = pi/2;
p180.Flip = pi;
E = Exp;
E.Sequence = {p90 0.2 p180 0.2};
E.nPoints = 16;
E.Dim1 = {'d1,d2' 0.01};  % µs
y0 = saffron(Sys,E,man);
for shift = [0.03 -0.03]
  E.DetWindow = shift;
  [~,yMan,info] = saffron(Sys,E);
  yFull = saffron(Sys,E,full);
  ok(end+1) = strcmp(info.Method,'manifold') && areequal(yMan,yFull,1e-4,'rel') && ...
    max(abs(yMan-y0))>0.1*max(abs(y0)); %#ok<AGROW>
end

% The same point from a longer last delay, or from a delay after a final
% pulse; two equal values are a single point
E.DetWindow = -0.05;
E.Sequence = {p90 0.2 p180 0.25};
y1 = saffron(Sys,E,man);
E.DetWindow = 0.2;
E.Sequence = {p90 0.2 p180};
E.Dim1 = {'d1' 0.01};
y2 = saffron(Sys,E,man);
E.DetWindow = [0.2 0.2];
y3 = saffron(Sys,E,man);
ok(end+1) = areequal(y1,y0,1e-12,'abs') && areequal(y2(1),y0(1),1e-12,'abs') && isequal(y3,y2);

% Invalid windows: more than two values, NaN, reversed start and end
E.DetWindow = 0.2 + 0.01*(0:15);
ok(end+1) = throwsError(@()saffron(Sys,E,man),'Exp.DetWindow must be');
E.DetWindow = NaN;
ok(end+1) = throwsError(@()saffron(Sys,E,man),'NaN');
E.DetWindow = [0.5 0];
ok(end+1) = throwsError(@()saffron(Sys,E,full),'start <= end');

% Echo times relative to the end of the sequence, also for the echo
% selection: the echo forms 0.2 µs after the final pulse
E = Exp;
E.Sequence = {p90 0.2 p180};
E.DetWindow = 0.2;
E.Pathways = 'echo';
[~,yMan,info] = saffron(Sys,E,man);
yFull = saffron(Sys,E,full);
ok(end+1) = isequal(info.Pathways,'+-') && areequal(info.EchoTimes,0.2,1e-12,'abs') && ...
  areequal(yMan,yFull,1e-4,'rel');

% Rectangular pulses, three-pulse ESEEM, and the detection frame
% (DetFreq): its phase includes the shift
p90f = struct('Flip',pi/2,'tp',0.01);  % µs
E = Exp;
E.Sequence = {p90f 0.15 p90f 0.1 p90f 0.15};
E.nPoints = 6;
E.Dim1 = {'d2' 0.04};  % µs
E.Pathways = '+0-';
E.DetWindow = 0.02;
[t,yMan,info] = saffron(Sys,E);
yFull = saffron(Sys,E,full);
ok(end+1) = strcmp(info.Method,'manifold') && areequal(yMan,yFull,1e-3,'rel');
delta = 0.003;  % GHz
E.DetFreq = Exp.mwFreq + delta;
yDet = saffron(Sys,E,man);
yDetFull = saffron(Sys,E,full);
tSeq = 3*0.01 + 0.15 + t + 0.15 + 0.02;  % µs, from the start of the sequence to detection (t: d2)
ok(end+1) = areequal(yDet,yMan.*exp(-2i*pi*delta*1e3*tSeq),1e-12,'abs') && ...
  areequal(yDet,yDetFull,1e-3,'rel');

% Relaxation: the shift takes the setting of the last delay. With only T2,
% the detected -1 coherence of the last interval decays by
% exp(-(tau+shift)/T2).
SysR = Sys;
SysR.T2 = 0.5;  % µs
E = Exp;
E.Sequence = {p90 0.2 p180 0.2};
E.DetWindow = 0.05;
OptOff = struct('Relaxation',[1 1 1 0]);
yOn = saffron(SysR,E,man);
yOff = saffron(SysR,E,setfield(OptOff,'Method','manifold')); %#ok<SFLD>
yOffFull = saffron(SysR,E,setfield(full,'Relaxation',OptOff.Relaxation)); %#ok<SFLD>
ok(end+1) = areequal(yOn,yOff*exp(-0.25/SysR.T2),1e-12,'rel') && areequal(yOff,yOffFull,1e-4,'rel');

% Transients and detection during the whole sequence: full engine
E = Exp;
E.Sequence = {p90 0.2 p180 0.2};
E.DetWindow = 'sequence';
[~,~,info] = saffron(Sys,E,rmfield(full,'Method'));
ok(end+1) = strcmp(info.Method,'full') && contains(info.MethodReason,'whole sequence');
E.DetWindow = [-0.01 0.01];
ok(end+1) = throwsError(@()saffron(Sys,E,man),'transient detection');

% A detection point before the last delay, or before a final pulse: error
% in both engines
E.DetWindow = -0.3;
ok(end+1) = throwsError(@()saffron(Sys,E,man),'into a pulse') && ...
  throwsError(@()saffron(Sys,E,full),'into a pulse');
E.Sequence = {p90 0.2 p180};
E.DetWindow = -0.1;
ok(end+1) = throwsError(@()saffron(Sys,E,man),'during a pulse');

end

function tf = throwsError(fun,msg)
try
  fun();
  tf = false;
catch err
  tf = contains(err.message,msg);
end
end
