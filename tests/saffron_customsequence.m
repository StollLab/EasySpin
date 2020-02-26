function [ok,data] = test(opt,olddata)

Sys.Nucs = '14N';
Sys.A = [-1 -1 2]*0.5+0.8;
Sys.Q = 1;

% User-defined HYSCORE experiment
tau = 0.080;  % us
% Exp.Flip = [1 1 2 1];
% Exp.t = [tau 0 0 tau];
% Exp.Inc = [0 1 2 0];
% Exp.dt = 0.025;
% 
% Exp.Field = 330;
% Exp.tau = 0.001;
% Exp.dt = 0.1;
% Exp.nPoints = 200;

p90.Flip = pi/2;
p180.Flip = pi;

Exp.Field = 330;
Exp.Sequence = {p90 tau p90 0 p180 0 p90 tau};
Exp.nPoints = [200 200];
Exp.Dim1 = {'d2' 0.1};
Exp.Dim2 = {'d3' 0.1};

Opt.Verbosity = 0;
y = saffron(Sys,Exp,Opt);

data.y = y;

if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-4,'abs');
else
  ok = [];
end
