function ok = test()

% Sweep type detection; NaN or empty range fields of the other sweep type
% are ignored; zero field only where allowed

Exp = struct('nPoints',1024,'mwFreq',9.5,'mwRange',NaN,'mwCenterSweep',[]);
ok(1) = runprivate('p_sweeptype',Exp,false);

Exp = struct('nPoints',1024,'Field',340,'mwFreq',[],'Range',NaN);
ok(2) = ~runprivate('p_sweeptype',Exp,false);

Exp = struct('nPoints',1024,'Field',0);
ok(3) = ~runprivate('p_sweeptype',Exp,true);

end
