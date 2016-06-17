function [err,data] = test(opt,olddata)

% Check flip angle to amplitude conversion for the pulse() function
%--------------------------------------------------------------------------

% Amplitude modulated pulses
Exp.tp = [0.060 0.200 0.200 0.100 0.500 0.300]; % us
Exp.PulseShape(1).Type = 'rectangular';
Exp.PulseShape(2).Type = 'gaussian';
Exp.PulseShape(2).tFWHM = 0.060; % us
Exp.PulseShape(3).Type = 'sinc';
Exp.PulseShape(3).zerocross = 0.050; % us
Exp.PulseShape(4).Type = 'quartersin';
Exp.PulseShape(4).trise = 0.020; % us
Exp.PulseShape(5).Type = 'sech';
Exp.PulseShape(5).beta = 12;
Exp.PulseShape(6).Type = 'WURST';
Exp.PulseShape(6).nwurst = 20;

Opt.Offsets = 0;

Exp.Flip = pi/2;
[~,~,p90] = pulse(Exp,Opt);

Exp.Flip = pi;
[~,~,p180] = pulse(Exp,Opt);

tol = 1e-12;
suberr = zeros(1,numel(Exp.tp));
for i = numel(Exp.tp)
  if p90(i).Mz<tol && (p180(i).Mz>(-1-tol) && p180(i).Mz<(-1+tol))  
    suberr(i) = 0;
  else
    suberr(i) = 1;
  end
end

err(1) = any(suberr);

% Amplitude and frequency modulated pulses
clearvars -except err
Exp.tp = [0.200 0.200 0.400 0.500 0.300 0.400]; % us
Exp.PulseShape(1).Type = 'quartersin/linear';
Exp.PulseShape(1).trise = 0.050; % us
Exp.PulseShape(1).BW = 500; % MHz
Exp.PulseShape(2).Type = 'WURST/linear';
Exp.PulseShape(2).nwurst = 30; % us
Exp.PulseShape(2).BW = 300; % MHz
Exp.PulseShape(3).Type = 'sech/tanh';
Exp.PulseShape(3).beta = 10;
Exp.PulseShape(3).BW = 70; % MHz
Exp.PulseShape(4).Type = 'sech*WURST/tanh';
Exp.PulseShape(4).beta = 4;
Exp.PulseShape(4).nwurst = 8;
Exp.PulseShape(4).BW = 120; % MHz
Exp.PulseShape(5).Type = 'sech/uniformQ';
Exp.PulseShape(5).beta = 10;
Exp.PulseShape(5).n = 4;
Exp.PulseShape(5).BW = 200; % MHz
Exp.PulseShape(6).Type = 'sech/uniformQ';
Exp.PulseShape(6).beta = 7;
Exp.PulseShape(6).n = 12;
Exp.PulseShape(6).BW = 250; % MHz

Opt.Offsets = 0;

Exp.Flip = pi/2;
[~,~,p90] = pulse(Exp,Opt);

Exp.Flip = pi;
[~,~,p180] = pulse(Exp,Opt);

tol = 1e-2;
suberr = zeros(1,numel(Exp.tp));
for i = numel(Exp.tp)
  if p90(i).Mz<tol && (p180(i).Mz>(-1-tol) && p180(i).Mz<(-1+tol))  
    suberr(i) = 0;
  else
    suberr(i) = 1;
  end
end

err(2) = any(suberr);

err = any(err);

data = [];