function [err,data] = test(opt,olddata)

Sys.S = 1/2;
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95];
Sys.Nucs = '14N';
Sys.lwpp = 10;
Sys.T1 = 2;
Sys.T2 = 0.5;

Exp.Field = 1240; % run the experiment at Q band
Exp.mwFreq = 34.78;

Chirp90.Type = 'quartersin/linear';
Chirp90.trise = 0.030;
Chirp90.tp = 0.200;
Chirp90.Flip = pi/2;
Chirp90.Frequency = [-300 300]; % excitation band, GHz

Chirp180.Type = 'quartersin/linear';
Chirp180.trise = 0.030;
Chirp180.tp = 0.100;
Chirp180.Flip = pi;
Chirp180.Frequency = [-300 300];

Exp.Sequence = {Chirp90 0.25 Chirp180 0.25}; 
Exp.DetWindow = [-0.1 0.1] + Chirp180.tp;
Exp.DetPhase = 0;

Opt.nKnots = 7;
Opt.SimulationMode = 'thyme';

[x, y] = saffron(Sys,Exp,Opt);

data.x = x;
data.y = y;

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,real(y),'r',x,real(olddata.y),'b');
    axis tight
    legend('new','old');
    title(mfilename);
    subplot(3,1,3);
    plot(x,real(olddata.y'-y));
    axis tight
    xlabel('time [us]');
  end
end

if ~isempty(olddata)
  err = ~areequal(y,olddata.y',1e-4);
else
  err = [];
end

