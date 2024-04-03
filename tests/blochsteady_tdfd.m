function ok = test(opt)

% Compare time-domain and frequency-domain methods
%-------------------------------------------------------

g = 2;  % g value
T1 = 1; % longitudinal relaxation time, µs
T2 = 1; % transverse relaxation time, µs
deltaB0 = 0.1; % in mT
B1 = 0.002; % microwave field, in mT
Bm = 0.3; % peak-to-peak field modulation amplitude, in mT
fm = 50; % field modulation frequency, in kHz

Options.Method = 'td';
[t1,My1] = blochsteady(g,T1,T2,deltaB0,B1,Bm,fm,Options);

Options.Method = 'fft';
[t2,My2] = blochsteady(g,T1,T2,deltaB0,B1,Bm,fm,Options);

ok = areequal(My1,My2,1e-10,'abs');

if (opt.Display)
  plot(t1,My1,'.',t2,My2);
  xlabel('time (\mus)');
  axis tight
  legend('td','fft');
  legend boxoff
  grid on
end
