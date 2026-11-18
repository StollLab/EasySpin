function ok = test() %#ok<FNDEF>

% Detection fields in the full-space engine: Exp.DetFreq defines the
% detection frame (also in the subspace engine), Exp.DetOperator gives other
% detection operators, and Exp.DetSequence is replaced by
% Exp.DetWindow = 'sequence'

Sys.S = 1/2;
Sys.g = 2;

Exp.Field = 350;
Exp.mwFreq = 9.8;

Pulse.Type = 'rectangular';
Pulse.tp = 0.016;  % µs
Pulse.Flip = pi/2;  % rad
Exp.Sequence = {Pulse 0.1};
Exp.DetWindow = [0 0.05];  % µs

Opt.GridSize = 7;
Opt.Method = 'fullspace';

% Exp.mwFreq is the default detection frequency
[~,y0] = saffron(Sys,Exp,Opt);
Exp.DetFreq = Exp.mwFreq;
[~,y1] = saffron(Sys,Exp,Opt);
ok(1) = areequal(y1,y0,1e-10,'rel');

% A different detection frequency shifts the signal by the frequency difference
Exp.DetFreq = Exp.mwFreq + 0.01;  % GHz
[t,y2] = saffron(Sys,Exp,Opt);
ratio = y2./y0;
ok(2) = areequal(ratio/ratio(1),exp(-2i*pi*10*(t-t(1))),1e-8,'abs');

% Exp.DetOperator: <Sz> after a pi pulse (initial state -Sz)
ExpD = Exp;
ExpD.Sequence{1}.Flip = pi;
ExpD.DetOperator = {'z1'};
ExpD.DetWindow = 0;
[~,yz] = saffron(Sys,ExpD,Opt);
ok(3) = isscalar(yz) && real(yz)>0.9;
ExpD = rmfield(Exp,'DetWindow');
ExpD.DetSequence = 'last';
ok(4) = throwsWith(@() saffron(Sys,ExpD,Opt),'DetSequence');

% The subspace engine gives a single point in the same detection frame
ExpF = rmfield(Exp,'DetWindow');
OptF = Opt;
OptF.Method = 'subspace';
yF = saffron(Sys,ExpF,Opt);
yS = saffron(Sys,ExpF,OptF);
ok(5) = areequal(yS,yF,1e-3,'rel');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
