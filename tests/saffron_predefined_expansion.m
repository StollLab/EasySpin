function ok = test() %#ok<FNDEF>

% A predefined experiment is a shorthand: the expanded experiment returned in
% info.Exp, run as a user-defined experiment, gives the same result, and the
% axis is the value of the first parameter in Exp.Dim1.

Sys.Nucs = '14N';
Sys.A = [1.2 0.3 0.4];  % MHz
Sys.Q = 0.8;  % MHz

Exp.Field = 350;  % mT
Exp.Sequence = '3pESEEM';
Exp.tau = 0.136;  % µs
Exp.T = 0.05;  % µs
Exp.dt = 0.016;  % µs
Exp.nPoints = 64;

Opt.GridSize = 7;

[x1,y1,info] = saffron(Sys,Exp,Opt);
ExpCustom = rmfield(info.Exp,'predefinedExp');
[x2,y2] = saffron(Sys,ExpCustom,Opt);

ok(1) = areequal(y1,y2,1e-12,'rel');
ok(2) = areequal(x1,x2,1e-12,'rel');
ok(3) = areequal(x1,Exp.T + Exp.dt*(0:Exp.nPoints-1),1e-12,'rel');
ok(4) = iscell(info.Exp.Sequence) && strcmp(info.Exp.Pathways,'+0-');
