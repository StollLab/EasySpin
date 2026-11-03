function ok = test() %#ok<FNDEF>

% Predefined HYSCORE: pulses along x and detection phase pi/2, giving a
% positive echo; a nonzero Exp.DetPhase and Exp.PhaseCycle are rejected

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz
Sys.AFrame = [0 0.6 0];

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',2);
Exp.tau = 0.12;  % µs
Exp.t1 = 0.07;  % µs
Exp.t2 = 0.23;  % µs
Exp.dt = 0.013;  % µs
Exp.nPoints = 3;

% Predefined HYSCORE gives a positive echo (generic timings, so that no
% other echo is detected)
Exp.Sequence = 'HYSCORE';
[~,yPredefined] = spidyan(Sys,Exp);
ok(1) = all(real(yPredefined(:))>0);

% Same as user-defined sequence with pulses along x and detection phase pi/2
P90 = struct('tp',0.016,'Flip',pi/2);
P180 = struct('tp',0.032,'Flip',pi);
ExpCustom = rmfield(Exp,{'tau','t1','t2','dt'});
ExpCustom.Sequence = {P90 Exp.tau P90 Exp.t1 P180 Exp.t2 P90 Exp.tau};
ExpCustom.Dim1 = {'d2' Exp.dt};
ExpCustom.Dim2 = {'d3' Exp.dt};
ExpCustom.nPoints = [3 3];
ExpCustom.DetWindow = 2*P90.tp/pi;  % echo of rectangular pulses
ExpCustom.DetPhase = pi/2;
sel = [0 1; pi/2 1i; pi -1; 3*pi/2 -1i];  % echo pathway after the first pulse
ExpCustom.PhaseCycle = {sel [] [0 1; pi 1] [0 1; pi -1]};
[~,yCustom] = spidyan(Sys,ExpCustom);
ok(2) = areequal(yCustom,yPredefined,1e-10,'rel');

% Exp.PhaseCycle cannot be used with predefined experiments
ExpPC = Exp;
ExpPC.PhaseCycle = {[0 1; pi -1]};
try
  spidyan(Sys,ExpPC);
  ok(4) = false;
catch err
  ok(4) = contains(err.message,'Exp.PhaseCycle cannot be used with predefined experiments');
end

% Exp.DetPhase cannot be changed for predefined experiments
Exp.DetPhase = pi/2;
try
  spidyan(Sys,Exp);
  ok(3) = false;
catch err
  ok(3) = contains(err.message,'Exp.DetPhase cannot be used with predefined experiments');
end
