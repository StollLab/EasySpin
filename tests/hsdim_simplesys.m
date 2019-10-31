function err = test(opt,olddata)

% Spin system with electrons only
%================================================================
Sys = struct('S',[1 3/2],'g',[3 4 5; 2 2 2],'J',100);

N = hsdim(Sys);

Nref = prod(2*Sys.S+1);

err = ~(N==Nref);
