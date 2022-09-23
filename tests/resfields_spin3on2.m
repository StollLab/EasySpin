function [ok,data] = test(opt,olddata)

% S=3/2 regression test

Sys = struct('S',3/2,'g',[2 2 2],'D',[-1 -1 2]*3e3);
Exp = struct('mwFreq',9.5,'Range',[0 2000]);

grid = sphgrid('Dinfh',31);
chi = zeros(size(grid.phi));
Exp.CrystalOrientation = [grid.phi(:) grid.theta(:) chi(:)];

[p,i] = resfields(Sys,Exp);

if opt.Display
  theta = grid.theta*180/pi;
  subplot(2,1,1)
  plot(theta,p,'r.',theta,olddata.p,'ko');
  xlabel('theta (deg)')
  ylabel('resonance field (mT)')
  subplot(2,1,2)
  plot(theta,i,'r.',theta,olddata.i,'ko');
  xlabel('theta (deg)')
  ylabel('intensity (arb.u.)')
end

data.p = p;
data.i = i;

if isempty(olddata)
  ok = [];
else
  idxp = ~isnan(p);
  idxi = ~isnan(i);
  thr = 1e-4;
  ok = areequal(p(idxp),olddata.p(idxp),thr,'abs') && ...
       areequal(i(idxi),olddata.i(idxi),thr,'abs');
end
