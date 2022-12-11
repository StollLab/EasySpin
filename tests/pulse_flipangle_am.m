function ok = test()

% Check flip angle to amplitude conversion for the pulse() function

% Amplitude modulated pulses
Params(1).Type = 'rectangular';
Params(1).tp = 0.060;  % µs

Params(2).Type = 'gaussian';
Params(2).tFWHM = 0.060;  % µs
Params(2).tp = 0.200;  % µs

Params(3).Type = 'sinc';
Params(3).zerocross = 0.050;  % µs
Params(3).tp = 0.200;  % µs

Params(4).Type = 'quartersin';
Params(4).trise = 0.020;  % µs
Params(4).tp = 0.100;  % µs

Params(5).Type = 'sech';
Params(5).beta = 12;  % 1/µs
Params(5).tp = 0.500;  % µs

Params(6).Type = 'WURST';
Params(6).nwurst = 20;
Params(6).tp = 0.300;  % µs

offsets = 0;

tol = 1e-12;
ok = false(1,numel(Params));
for i = 1:numel(Params)
  Params(i).Flip = pi/2;
  [t,IQ] = pulse(Params(i));
  [offsets,M90] = exciteprofile(t,IQ,offsets);
  Params(i).Flip = pi;
  [t,IQ] = pulse(Params(i));
  [offsets,M180] = exciteprofile(t,IQ,offsets);
  ok(i) = M90(3)<tol && M180(3)>(-1-tol) && M180(3)<(-1+tol);  
end
