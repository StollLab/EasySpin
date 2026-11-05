function ok = test() %#ok<FNDEF>

% info.EventTimes gives start and end of each element of Exp.Sequence, on
% the time axis of the transient: from the start of the sequence for
% DetWindow = 'sequence', from the end of the sequence for detection
% windows. Cropping the whole-sequence signal to the last delay gives the
% same signal as a detection window over that delay.

Sys.S = 1/2;
Sys.g = gfree*mhz2mt(9505)/mhz2mt(9500);  % 5 MHz above the carrier

P1.Type = 'rectangular';
P1.tp = 0.02;  % µs
P1.Flip = pi/2;
P2 = P1;
P2.tp = 0.04;  % µs
P2.Flip = pi;

Exp.Field = mhz2mt(9500);  % mT
Exp.mwFreq = 9.5;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetOperator = {'+1'};

tol = 1e-9;

% Whole sequence
E = Exp;
E.Sequence = {P1 0.1 P2 0.2};
E.DetWindow = 'sequence';
[t,S,info] = saffron(Sys,E);
ok(1) = areequal(info.EventTimes,[0 0.02; 0.02 0.12; 0.12 0.16; 0.16 0.36],tol,'abs');

% Detection window over the last delay (sequence ends with the pulse)
E = Exp;
E.Sequence = {P1 0.1 P2};
E.DetWindow = [0 0.2];  % µs
[tW,SW,infoW] = saffron(Sys,E);
ok(2) = areequal(infoW.EventTimes,[-0.16 -0.14; -0.14 -0.04; -0.04 0],tol,'abs');

% Cropping the whole-sequence signal to the last delay
idx = t >= info.EventTimes(4,1)-tol & t <= info.EventTimes(4,2)+tol;
tCrop = t(idx) - info.EventTimes(4,1);
SCrop = S(idx);
ok(3) = numel(tCrop)==numel(tW) && areequal(tCrop,tW,tol,'abs');
ok(4) = numel(SCrop)==numel(SW) && areequal(SCrop,SW,1e-6,'abs');

% Single point at the end: times relative to the end of the sequence
E = Exp;
E.Sequence = {P1 0.1 P2 0.2};
E.DetWindow = 0;
[~,~,info0] = saffron(Sys,E);
ok(5) = areequal(info0.EventTimes(end,:),[-0.2 0],tol,'abs');

end
