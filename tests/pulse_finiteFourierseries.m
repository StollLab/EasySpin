function ok = test()

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Fourier series pulses
% I-BURP 1
Params.tp = 0.500; % µs
Params.Type = 'I-BURP 1';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Geen, H., Freeman, R. Band-selective radiofrequency pulses. 
% J. Magn. Reson. 93, 93-141 (1991). 
% DOI: 10.1016/0022-2364(91)90034-Q
% Table 5 on p. 117, Np = 256        
A0 = 0.5;
An = [0.70 -0.15 -0.94 0.11 -0.02 -0.04 0.01 -0.02 -0.01];
Bn = [-1.54 1.01 -0.24 -0.04 0.08 -0.04 -0.01 0.01 -0.01];
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = 1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% I-BURP 2
Params.tp = 0.700; % µs
Params.Type = 'I-BURP 2';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Geen, H., Freeman, R. Band-selective radiofrequency pulses. 
% J. Magn. Reson. 93, 93-141 (1991). 
% DOI: 10.1016/0022-2364(91)90034-Q
% Table 6 on p. 119, Np = 256
A0 = 0.5;
An = [0.81 0.07 -1.25 -0.24 0.07 0.11 0.05 -0.02 -0.03 -0.02 0.00];
Bn = [-0.68 -1.38 +0.20 0.45 0.23 0.05 -0.04 -0.04 0.00 0.01 0.01];
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% E-BURP 1
Params.tp = 0.500; % µs
Params.Type = 'E-BURP 1';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Geen, H., Freeman, R. Band-selective radiofrequency pulses. 
% J. Magn. Reson. 93, 93-141 (1991). 
% DOI: 10.1016/0022-2364(91)90034-Q
% Table 3 on p. 112, Np = 256
A0 = 0.23;
An = [0.88 -1.04 -0.24 0.14 0.03 0.04 -0.03 0.00];
Bn = [-0.40 -1.42 0.77 0.06 0.03 -0.04 -0.02 0.01];
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% E-BURP 2
Params.tp = 1.000; % µs
Params.Type = 'E-BURP 2';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Geen, H., Freeman, R. Band-selective radiofrequency pulses. 
% J. Magn. Reson. 93, 93-141 (1991). 
% DOI: 10.1016/0022-2364(91)90034-Q
% Table 4 on p. 115, Np = 256
A0 = 0.26;
An = [0.91 0.29 -1.28 -0.05 0.04 0.02 0.06 0.00 -0.02 0.00];
Bn = [-0.16 -1.82 0.18 0.42 0.07 0.07 -0.01 -0.04 0.00 0.00];
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% U-BURP
Params.tp = 1.000; % µs
Params.Type = 'U-BURP';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Geen, H., Freeman, R. Band-selective radiofrequency pulses. 
% J. Magn. Reson. 93, 93-141 (1991). 
% DOI: 10.1016/0022-2364(91)90034-Q
% Table 7 on p. 122, Np = 256
A0 = 0.27;
An = [-1.42 -0.37 -1.84 4.40 -1.19 0 -0.37 0.50 -0.31 0.18 -0.21 0.23 -0.12 0.07 -0.06 0.06 -0.04 0.03 -0.02 0.02];
Bn = zeros(size(An));
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% RE-BURP
Params.tp = 1.000; % µs
Params.Type = 'RE-BURP';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Geen, H., Freeman, R. Band-selective radiofrequency pulses. 
% J. Magn. Reson. 93, 93-141 (1991). 
% DOI: 10.1016/0022-2364(91)90034-Q
% Table 8 on p. 124, Np = 256
A0 = 0.49;
An = [-1.02 1.11 -1.57 0.83 -0.42 0.26 -0.16 0.10 -0.07 0.04 -0.03 0.01 -0.02 0.00 -0.01];
Bn = zeros(size(An));
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% SNOB i2
Params.tp = 0.200; % µs
Params.Type = 'SNOB i2';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Kupce, E., Boyd, J., Campbell, I. D. Short Selective Pulses for
% Biochemical Applications. J. Magn. Reson. B 106, 300-303 (1995).
% DOI: 10.1006/jmrb.1995.1049
% Table 1 on p. 300
A0 = 0.5;
An = [-0.2687 -0.2972 0.0989 -0.0010 -0.0168 0.0009 -0.0017 -0.0013 -0.0014];
Bn = [-1.1461 0.4016 0.0736 -0.0307 0.0079 0.0062 0.0003 -0.0002 0.0009];
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% SNOB i3
Params.tp = 0.200; % µs
Params.Type = 'SNOB i3';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Kupce, E., Boyd, J., Campbell, I. D. Short Selective Pulses for
% Biochemical Applications. J. Magn. Reson. B 106, 300-303 (1995).
% DOI: 10.1006/jmrb.1995.1049
% Table 1 on p. 300
A0 = 0.5;
An = [0.2801 -0.9995 0.1928 0.0967 -0.0480 -0.0148 0.0088 -0.0002 -0.0030];
Bn = [-1.1990 0.4893 0.2439 -0.0816 -0.0409 0.0234 0.0036 -0.0042 0.0001];
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');

% User-defined (SLURP-1, T2/T = 2)
Params.tp = 2.000; % µs
Params.Type = 'FourierSeries';
Params.TimeStep = 0.001; % µs
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
% Nuzillard, J. M., Freeman, R. Band-Selective Pulses Designed to 
% Accommodate Relaxation. J. Magn. Reson. A. 107, 113-118 (1994).
% DOI: 10.1006/jmra.1994.1056
% Table 1
A0 = 0.308;
An = [ 1.017 -0.480 -1.033 0.078 0.103 0.109  0.027 -0.043 -0.018 0.000 0.005  0.004];
Bn = [-0.384 -1.894  0.574 0.409 0.098 0.009 -0.079 -0.024  0.014 0.010 0.003 -0.001];
A = zeros(1,numel(t0)) + A0;
for j = 1:numel(An)
  A = A + An(j)*cos(j*2*pi*t0/Params.tp) + Bn(j)*sin(j*2*pi*t0/Params.tp);
end
A = A/max([-min(A) max(A)]);
IQ0 = A/max(A);

Params.A0 = A0;
Params.An = An;
Params.Bn = Bn;

[~,IQ] = pulse(Params);

ipulse = ipulse+1;
ok(ipulse) = areequal(IQ0,IQ,1e-12,'abs');
