function [ok,data] = test(opt,olddata)

% Spin System
Sys.S = 1/2;
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95];
Sys.Nucs = '14N';
Sys.lwpp = 10;

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

Exp.nPoints = 2;
Exp.Dim1 = {'d1,d2' 0.5};

Opt.GridSize = 7;
Opt.SimulationMode = 'thyme';

[x1, y1] = saffron(Sys,Exp,Opt);

data.x1 = x1;
data.y1 = y1;

Exp.nPoints = [3 2];
Exp.Dim1 = {'d1,d2' 0.5};
Exp.Dim2 = {'p1.Flip' pi/3};

[x2, y2] = saffron(Sys,Exp,Opt);

data.x2 = x2;
data.y2 = y2;

if opt.Display
  if ~isempty(olddata)
    p1 = subplot(3,1,[1 2]);
    plot(x1{2},real(y1(1,:)),x1{2},real(olddata.y1(1,:)));
    axis tight
    legend('new','old');
    title(mfilename);
    p2 = subplot(3,1,3);
    plot(x1{2},real(olddata.y1(1,:)-y1(1,:)));
    axis tight
    xlabel('time (µs)');
    title('old - new')
    linkaxes([p1,p2],'x')
  end
end

if ~isempty(olddata)
  ok = areequal(y1,olddata.y1,1e-3,'abs') && areequal(y2,olddata.y2,1e-3,'abs');
else
  ok = [];
end
