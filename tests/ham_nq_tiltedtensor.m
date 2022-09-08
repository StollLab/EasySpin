function ok = test()

% tilted quadrupole tensor

Sys = struct('S',1/2,'g',[2 2 2],'Nucs','2H','A',[1 2 3]);
Sys.Q = [-1.5 -0.5 2];
Sys.QFrame = [30 50 80]*pi/180;
HQ = ham_nq(Sys);

I{1} = sop(Sys,'ex');
I{2} = sop(Sys,'ey');
I{3} = sop(Sys,'ez');

R = erot(Sys.QFrame);
QTens = R.'*diag(Sys.Q)*R;

HQ1 = 0;
for k = 1:3
  for q = 1:3
    HQ1 = HQ1 + QTens(k,q)*I{k}*I{q};
  end
end

ok = areequal(HQ,HQ1,1e-10,'abs');
