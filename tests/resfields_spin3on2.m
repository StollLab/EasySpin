function [err,data] = test(opt,olddata)

% S=3/2 regression test

Sys = struct('S',3/2,'g',[2 2 2],'D',[-1 -1 2]*3e3);
Exp = struct('mwFreq',9.5,'Range',[0 2000]);

[phi,theta] = sphgrid('Dinfh',31);
chi = zeros(size(phi));
Exp.CrystalOrientation = [phi(:) theta(:) chi(:)];

[p,i] = resfields(Sys,Exp);

if opt.Display
  plot(theta,p,'.');
end

data.p = p;
data.i = i;

if isempty(olddata)
  err = [];
else
  idxp = ~isnan(p);
  idxi = ~isnan(i);
  thr = 1e-4;
  ok = areequal(p(idxp),olddata.p(idxp),thr,'abs') && ...
       areequal(i(idxi),olddata.i(idxi),thr,'abs');
  err = ~ok;
end
