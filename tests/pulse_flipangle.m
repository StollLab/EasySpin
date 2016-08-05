function [err,data] = test(opt,olddata)

% Check flip angle to amplitude conversion for the pulse() function
%--------------------------------------------------------------------------

% Amplitude modulated pulses
Params(1).Type = 'rectangular';
Params(1).tp = 0.060; % us

Params(2).Type = 'gaussian';
Params(2).tFWHM = 0.060; % us
Params(2).tp = 0.200; % us

Params(3).Type = 'sinc';
Params(3).zerocross = 0.050; % us
Params(3).tp = 0.200; % us

Params(4).Type = 'quartersin';
Params(4).trise = 0.020; % us
Params(4).tp = 0.100; % us

Params(5).Type = 'sech';
Params(5).beta = 12;
Params(5).tp = 0.500; % us

Params(6).Type = 'WURST';
Params(6).nwurst = 20;
Params(6).tp = 0.300; % us

offsets = 0;

tol = 1e-12;
suberr = zeros(1,numel(Params));
for i = 1:numel(Params)
  Params(i).Flip = pi/2;
  [t,IQ] = pulse(Params(i));
  [offsets,M90] = exciteprofile(t,IQ,offsets);
  Params(i).Flip = pi;
  [t,IQ] = pulse(Params(i));
  [offsets,M180] = exciteprofile(t,IQ,offsets);
  if M90(3)<tol && (M180(3)>(-1-tol) && M180(3)<(-1+tol))  
    suberr(i) = 0;
  else
    suberr(i) = 1;
  end
end

err(1) = any(suberr);

% Amplitude and frequency modulated pulses
clear Params
Params(1).Type = 'quartersin/linear';
Params(1).tp = 0.200; % us
Params(1).trise = 0.050; % us
Params(1).Frequency = [-250 250]; % MHz
Params(2).Type = 'WURST/linear';
Params(2).tp = 0.200; % us
Params(2).nwurst = 30; % us
Params(2).Frequency = [150 -150]; % MHz
Params(3).Type = 'sech/tanh';
Params(3).tp = 0.400; % us
Params(3).beta = 10;
Params(3).Frequency = [-35 35]; % MHz
Params(4).Type = 'sech*WURST/tanh';
Params(4).tp = 0.500; % us
Params(4).beta = 4;
Params(4).nwurst = 8;
Params(4).Frequency = [60 -60]; % MHz
Params(5).Type = 'sech/uniformQ';
Params(5).tp = 0.300; % us
Params(5).beta = 10;
Params(5).n = 4;
Params(5).Frequency = [-100 100]; % MHz
Params(6).Type = 'sech/uniformQ';
Params(6).tp = 0.400; % us
Params(6).beta = 7;
Params(6).n = 12;
Params(6).Frequency = [-125 125]; % MHz

offsets = 0;

tol = 1e-2;
suberr = zeros(1,numel(Params));
for i = 1:numel(Params)
  Params(i).Flip = pi/2;
  [t,IQ] = pulse(Params(i));  
  [offsets,M90] = exciteprofile(t,IQ,offsets);
  Params(i).Flip = pi;
  [t,IQ] = pulse(Params(i));
  [offsets,M90] = exciteprofile(t,IQ,offsets);
  if M90(3)<tol && (M180(3)>(-1-tol) && M180(3)<(-1+tol))  
    suberr(i) = 0;
  else
    suberr(i) = 1;
  end
end

err(2) = any(suberr);

err = any(err);

data = [];