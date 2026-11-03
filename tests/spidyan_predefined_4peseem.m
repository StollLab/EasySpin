function ok = test() %#ok<FNDEF>

% Predefined 4pESEEM: pulses along y and detection phase pi, giving a
% positive echo, same as the equivalent user-defined sequence

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz
Sys.AFrame = [0 0.6 0];

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',2);
Exp.tau = 0.12;  % µs
Exp.T = 0.07;  % µs
Exp.dt = 0.013;  % µs
Exp.nPoints = 3;

% Predefined 4pESEEM gives a positive echo
Exp.Sequence = '4pESEEM';
[~,yPredefined] = spidyan(Sys,Exp);
ok(1) = all(real(yPredefined(:))>0);

% Same as user-defined sequence with pulses along y and detection phase pi
P90 = struct('tp',0.016,'Flip',pi/2,'Phase',pi/2);
P180 = struct('tp',0.032,'Flip',pi,'Phase',pi/2);
ExpCustom = rmfield(Exp,{'tau','T','dt'});
ExpCustom.Sequence = {P90 Exp.tau P90 Exp.T P180 Exp.T P90 Exp.tau};
ExpCustom.Dim1 = {'d2,d3' Exp.dt};
ExpCustom.DetWindow = 2*P90.tp/pi;  % echo of rectangular pulses
ExpCustom.DetPhase = pi;
sel = [0 1; pi/2 1i; pi -1; 3*pi/2 -1i];  % echo pathway after the first pulse
ExpCustom.PhaseCycle = {sel [] [0 1; pi 1] [0 1; pi -1]};
[~,yCustom] = spidyan(Sys,ExpCustom);
ok(2) = areequal(yCustom,yPredefined,1e-10,'rel');
