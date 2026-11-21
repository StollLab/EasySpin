function ok = test() %#ok<FNDEF>

% Options that saffron does not support give an error instead of being
% ignored, with any Opt.Method: a pulse field Operator, Opt.StateTrajectories,
% Sys.lwAnisoG and Opt.IdealPulses when the full engine is needed,
% and Opt.DetOperator

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

p90.Flip = pi/2;

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.Sequence = {p90 0.1 p90 0.1};

auto = struct;
full = struct('Method','full');

% Pulse field Operator
E = Exp;
E.Sequence{1}.Operator = 'x1';
ok(1) = throwsWith(@() saffron(Sys,E,auto),'Operator');
ok(2) = throwsWith(@() saffron(Sys,E,full),'Operator');

% Opt.StateTrajectories
ok(3) = throwsWith(@() saffron(Sys,Exp,struct('StateTrajectories',true)),'StateTrajectories');
ok(4) = throwsWith(@() saffron(Sys,Exp,setfield(full,'StateTrajectories',true)),'StateTrajectories'); %#ok<SFLD>

% Sys.lwAnisoG: the manifold engine with 'auto', an error when the
% full engine is needed or requested
SysL = Sys;
SysL.lwAnisoG = 5;  % MHz
[~,~,info] = saffron(SysL,Exp,auto);
ok(5) = strcmp(info.Method,'manifold');
ok(6) = throwsWith(@() saffron(SysL,Exp,full),'lwAnisoG');
E = Exp;
E.DetOperator = {'z1'};  % needs the full engine
ok(7) = throwsWith(@() saffron(SysL,E,auto),'lwAnisoG');

% Opt.DetOperator (detection operators belong in Exp.DetOperator)
ok(8) = throwsWith(@() saffron(Sys,Exp,struct('DetOperator','z1')),'Exp.DetOperator');
ok(9) = throwsWith(@() saffron(Sys,Exp,setfield(full,'DetOperator',{'z1'})),'Exp.DetOperator'); %#ok<SFLD>

% Opt.IdealPulses: the manifold engine with 'auto', an error when the
% full engine is needed or requested
EF = Exp;
EF.Sequence = {struct('Flip',pi/2,'tp',0.016) 0.1 struct('Flip',pi/2,'tp',0.016) 0.1};
[~,~,info] = saffron(Sys,EF,struct('IdealPulses',true));
ok(10) = strcmp(info.Method,'manifold');
ok(11) = throwsWith(@() saffron(Sys,EF,setfield(full,'IdealPulses',true)),'IdealPulses'); %#ok<SFLD>
E = EF;
E.DetOperator = {'z1'};  % needs the full engine
ok(12) = throwsWith(@() saffron(Sys,E,struct('IdealPulses',true)),'IdealPulses');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
