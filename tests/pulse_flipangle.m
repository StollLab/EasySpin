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

Opt.Offsets = 0;


tol = 1e-12;
suberr = zeros(1,numel(Params));
for i = 1:numel(Params)
  Params(i).Flip = pi/2;
  [~,~,p90] = pulse(Params(i),Opt);
  Params(i).Flip = pi;
  [~,~,p180] = pulse(Params(i),Opt);
  if p90.Mz<tol && (p180.Mz>(-1-tol) && p180.Mz<(-1+tol))  
    suberr(i) = 0;
  else
    suberr(i) = 1;
  end
end

err(1) = any(suberr);

% Amplitude and frequency modulated pulses
clearvars -except err
Params(1).Type = 'quartersin/linear';
Params(1).tp = 0.200; % us
Params(1).trise = 0.050; % us
Params(1).BW = 500; % MHz
Params(2).Type = 'WURST/linear';
Params(2).tp = 0.200; % us
Params(2).nwurst = 30; % us
Params(2).BW = 300; % MHz
Params(3).Type = 'sech/tanh';
Params(3).tp = 0.400; % us
Params(3).beta = 10;
Params(3).BW = 70; % MHz
Params(4).Type = 'sech*WURST/tanh';
Params(4).tp = 0.500; % us
Params(4).beta = 4;
Params(4).nwurst = 8;
Params(4).BW = 120; % MHz
Params(5).Type = 'sech/uniformQ';
Params(5).tp = 0.300; % us
Params(5).beta = 10;
Params(5).n = 4;
Params(5).BW = 200; % MHz
Params(6).Type = 'sech/uniformQ';
Params(6).tp = 0.400; % us
Params(6).beta = 7;
Params(6).n = 12;
Params(6).BW = 250; % MHz

Opt.Offsets = 0;

tol = 1e-2;
suberr = zeros(1,numel(Params));
for i = 1:numel(Params)
  Params(i).Flip = pi/2;
  [~,~,p90] = pulse(Params(i),Opt);  
  Params(i).Flip = pi;
  [~,~,p180] = pulse(Params(i),Opt);
  if p90.Mz<tol && (p180.Mz>(-1-tol) && p180.Mz<(-1+tol))  
    suberr(i) = 0;
  else
    suberr(i) = 1;
  end
end

err(2) = any(suberr);

err = any(err);

data = [];