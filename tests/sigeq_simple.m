function [err,data] = test(opt,olddata)

%======================================================================
% Test 2: Regression test
%======================================================================
Sys = struct('S',3/2,'g',[2 2 2],'D',[-1 -1 2]*3000);
H = sham(Sys,[0 0 350]);
Temps = [1 2 5 10];
clear s
for k=1:numel(Temps)
  s{k} = sigeq(H,Temps(k));
end

data.s = s;

if ~isempty(olddata)
  ok = 1;
  for k=1:numel(s)
    ok = ok & areequal(s{k},olddata.s{k},1e-6);
  end
  err = ~ok;
else
  err = [];
end

