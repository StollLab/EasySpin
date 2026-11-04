function ok = test()

% Negative magnetic field is rejected, and Sys.ZeemanFreq is no longer
% supported

Pulse.Type = 'rectangular';
Pulse.tp = 0.016;  % µs
Pulse.Flip = pi/2;  % rad

Exp.Sequence = {Pulse 0.1};
Exp.mwFreq = 9.5;  % GHz

Sys.S = 1/2;
Sys.g = 2;
Exp.Field = -340;  % mT
ok(1) = throwsWith(@() spidyan(Sys,Exp),'Exp.Field cannot be negative');

Sys = struct('S',1/2,'ZeemanFreq',9.5);  % GHz
Exp.Field = 340;
ok(2) = throwsWith(@() spidyan(Sys,Exp),'ZeemanFreq');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
