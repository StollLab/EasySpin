function ok = test()

% Compare time-domain and frequency-domain methods
%-------------------------------------------------------

g = 2;  % g value
T1 = 1; % longitudinal relaxation time, µs
T2 = 1; % transverse relaxation time, µs
deltaB0 = -0.1; % in mT
B1 = 0.002; % microwave field, in mT
Bm = 0.3; % peak-to-peak field modulation amplitude, in mT
fm = 50; % field modulation frequency, in kHz

[t1,My] = blochsteady(g,T1,T2,deltaB0,B1,Bm,fm);
[t1,Mx,My,Mz] = blochsteady(g,T1,T2,deltaB0,B1,Bm,fm);

Options.Method = 'fft';
[t1,My] = blochsteady(g,T1,T2,deltaB0,B1,Bm,fm,Options);
[t1,Mx,My,Mz] = blochsteady(g,T1,T2,deltaB0,B1,Bm,fm,Options);

ok = true;
