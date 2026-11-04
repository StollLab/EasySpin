function ok = test() %#ok<FNDEF>

% Options that saffron does not support give an error instead of being
% ignored, with any Opt.Method: a pulse field Operator, Opt.StateTrajectories,
% and Sys.lwAnisoG when the full-space engine is needed

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

p90.Flip = pi/2;

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.Sequence = {p90 0.1 p90 0.1};

auto = struct;
full = struct('Method','fullspace');

% Pulse field Operator
E = Exp;
E.Sequence{1}.Operator = 'x1';
ok(1) = throwsWith(@() saffron(Sys,E,auto),'Operator');
ok(2) = throwsWith(@() saffron(Sys,E,full),'Operator');

% Opt.StateTrajectories
ok(3) = throwsWith(@() saffron(Sys,Exp,struct('StateTrajectories',true)),'StateTrajectories');
ok(4) = throwsWith(@() saffron(Sys,Exp,setfield(full,'StateTrajectories',true)),'StateTrajectories'); %#ok<SFLD>

% Sys.lwAnisoG: the subspace engine with 'auto', an error when the
% full-space engine is needed or requested
SysL = Sys;
SysL.lwAnisoG = 5;  % MHz
[~,~,info] = saffron(SysL,Exp,auto);
ok(5) = strcmp(info.Method,'subspace');
ok(6) = throwsWith(@() saffron(SysL,Exp,full),'lwAnisoG');
E = Exp;
E.DetOperator = {'z1'};  % needs the full-space engine
ok(7) = throwsWith(@() saffron(SysL,E,auto),'lwAnisoG');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
