function ok = test()

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Gaussian pulse cascades
% G3
clear Params
Params.tp = 0.800; % us
Params.Type = 'G3';
Params.TimeStep = 0.001; % us
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Emsley, L., Bodenhausen, G. Gaussian pulse cascades: New analytical
% functions for rectangular selective inversion and in-phase excitation
% in NMR. Chem. Phys. Lett. 165, 469–476 (1990).
% DOI: 10.1016/0009-2614(90)87025-M
% Table 1 on p. 473, Cycle 3
A0 = [-1 1.37 0.49];
x0 = [0.287 0.508 0.795]*Params.tp;
FWHM = [0.189 0.183 0.243]*Params.tp;
A = zeros(1,numel(t0));
for j = 1:3
  A = A + A0(j)*exp(-(4*log(2)/FWHM(j)^2)*(t0-x0(j)).^2);
end
IQ0 = A/max(A);

[t,IQ] = pulse(Params);

ok(1) = areequal(IQ0,IQ,1e-12,'abs');

% Gaussian pulse cascades
% G4
clear Params
Params.tp = 0.800; % us
Params.Type = 'G4';
Params.TimeStep = 0.001; % us
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Emsley, L., Bodenhausen, G. Gaussian pulse cascades: New analytical
% functions for rectangular selective inversion and in-phase excitation
% in NMR. Chem. Phys. Lett. 165, 469–476 (1990).
% DOI: 10.1016/0009-2614(90)87025-M
% Table 1 on p. 474, Cycle 11
A0 = [0.62 0.72 -0.91 -0.33];
x0 = [0.177 0.492 0.653 0.892]*Params.tp;
FWHM = [0.172 0.129 0.119 0.139]*Params.tp;
A = zeros(1,numel(t0));
for j = 1:4
  A = A + A0(j)*exp(-(4*log(2)/FWHM(j)^2)*(t0-x0(j)).^2);
end
IQ0 = A/max(A);

[t,IQ] = pulse(Params);

ok(2) = areequal(IQ0,IQ,1e-12,'abs');

% Q3
clear Params
Params.tp = 1.200; % us
Params.Type = 'Q3';
Params.TimeStep = 0.01; % us
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Emsley, L., Bodenhausen, G. Optimization of shaped selective pulses
% for NMR using a quaternion description of their overall propagators.
% J. Magn. Reson. 97, 135–148 (1992).
% DOI: 10.1016/0022-2364(92)90242-Y
% Table 2 on p. 142
A0 = [-4.39 4.57 2.60];
x0 = [0.306 0.545 0.804]*Params.tp;
FWHM = [0.180 0.183 0.245]*Params.tp;
A = zeros(1,numel(t0));
for j = 1:3
  A = A + A0(j)*exp(-(4*log(2)/FWHM(j)^2)*(t0-x0(j)).^2);
end
IQ0 = A/max(A);

[t,IQ] = pulse(Params);

ok(3) = areequal(IQ0,IQ,1e-12,'abs');

% Q5
clear Params
Params.tp = 1.200; % us
Params.Type = 'Q5';
Params.TimeStep = 0.01; % us
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Emsley, L., Bodenhausen, G. Optimization of shaped selective pulses
% for NMR using a quaternion description of their overall propagators.
% J. Magn. Reson. 97, 135–148 (1992).
% DOI: 10.1016/0022-2364(92)90242-Y
% Table 2 on p. 142
A0 = [-1.48 -4.34 7.33 -2.30 5.66];
x0 = [0.162 0.307 0.497 0.525 0.803]*Params.tp;
FWHM = [0.186 0.139 0.143 0.290 0.137]*Params.tp;
A = zeros(1,numel(t0));
for j = 1:5
  A = A + A0(j)*exp(-(4*log(2)/FWHM(j)^2)*(t0-x0(j)).^2);
end
IQ0 = A/max(A);

[t,IQ] = pulse(Params);

ok(4) = areequal(IQ0,IQ,1e-12,'abs');

% User-defined
clear Params
Params.tp = 0.600; % us
Params.Type = 'GaussianCascade';
Params.TimeStep = 0.01; % us
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% random values
A0 = [-0.84 1.53 -2.12 5.21];
x0 = [0.154 0.358 0.521 0.766]*Params.tp;
FWHM = [0.180 0.183 0.195 0.235]*Params.tp;
A = zeros(1,numel(t0));
for j = 1:4
  A = A + A0(j)*exp(-(4*log(2)/FWHM(j)^2)*(t0-x0(j)).^2);
end
IQ0 = A/max(A);

Params.A0 = A0;
Params.x0 = x0/Params.tp;
Params.FWHM = FWHM/Params.tp;

[t,IQ] = pulse(Params);

ok(5) = areequal(IQ0,IQ,1e-12,'abs');
