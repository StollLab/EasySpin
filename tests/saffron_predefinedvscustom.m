function ok = test() %#ok<FNDEF>

% Fast mode: predefined experiments give the same signal as the same
% sequence entered as a custom sequence (same coherence pathways), including
% the sign of the pulse transfer prefactors. Predefined experiments use
% pulses along x with Exp.DetPhase = -pi/2 (2pESEEM, 3pESEEM) or +pi/2
% (4pESEEM, HYSCORE); equivalently, pulses along y with Exp.DetPhase = 0
% or pi.

Sys.S = 1/2;
Sys.g = [2.00 2.02 2.04];
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Base.mwFreq = 9.5;  % GHz
Base.Field = unitconvert(Base.mwFreq,'GHz->mT',2.02);
Base.ExciteWidth = Inf;

Opt.GridSize = 7;
Opt.SimulationMode = 'fast';

P90.Flip = pi/2;
P180.Flip = pi;
tau = 0.15;  % µs
T = 0.1;  % µs
dt = 0.02;  % µs

% 2pESEEM
E = Base; E.Sequence = '2pESEEM'; E.tau = tau; E.dt = dt; E.nPoints = 8;
yPredefined = saffron(Sys,E,Opt);
E = Base; E.Sequence = {P90 tau P180 tau}; E.Dim1 = {'d1,d2' dt}; E.nPoints = 8;
E.DetPhase = -pi/2;
yCustom = saffron(Sys,E,Opt);
ok(1) = areequal(yCustom,yPredefined,1e-10,'rel');

% 3pESEEM
E = Base; E.Sequence = '3pESEEM'; E.tau = tau; E.T = T; E.dt = dt; E.nPoints = 8;
yPredefined = saffron(Sys,E,Opt);
E = Base; E.Sequence = {P90 tau P90 T P90 tau}; E.Dim1 = {'d2' dt}; E.nPoints = 8;
E.DetPhase = -pi/2;
yCustom = saffron(Sys,E,Opt);
ok(2) = areequal(yCustom,yPredefined,1e-10,'rel');

% 4pESEEM
E = Base; E.Sequence = '4pESEEM'; E.tau = tau; E.T = T; E.dt = dt; E.nPoints = 8;
yPredefined = saffron(Sys,E,Opt);
E = Base; E.Sequence = {P90 tau P90 T P180 T P90 tau}; E.Dim1 = {'d2,d3' dt}; E.nPoints = 8;
E.Pathways = ['+ab-';'+ba-'];
E.DetPhase = pi/2;
yCustom = saffron(Sys,E,Opt);
ok(3) = areequal(yCustom,yPredefined,1e-10,'rel');

% HYSCORE
E = Base; E.Sequence = 'HYSCORE'; E.tau = tau; E.t1 = T; E.t2 = T; E.dt = dt; E.nPoints = [4 4];
yPredefined = saffron(Sys,E,Opt);
E = Base; E.Sequence = {P90 tau P90 T P180 T P90 tau}; E.Dim1 = {'d2' dt}; E.Dim2 = {'d3' dt};
E.nPoints = [4 4];
E.DetPhase = pi/2;
yCustom = saffron(Sys,E,Opt);
ok(4) = areequal(yCustom,yPredefined,1e-10,'rel');

% HYSCORE with pulses along y and Exp.DetPhase = pi
P90y = P90; P90y.Phase = pi/2;
P180y = P180; P180y.Phase = pi/2;
E.Sequence = {P90y tau P90y T P180y T P90y tau};
E.DetPhase = pi;
yCustom = saffron(Sys,E,Opt);
ok(5) = areequal(yCustom,yPredefined,1e-10,'rel');

% Predefined experiments give a real positive echo
ok(6) = all(abs(imag(yPredefined(:)))<1e-10*abs(yPredefined(:))) && all(real(yPredefined(:))>0);
