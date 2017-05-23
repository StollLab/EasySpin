function [err,data] = test(opt,olddata)

Sys = struct('S',1/2,'g',[2.01 2.2]);
for k = 1:15
  Sys = nucspinadd(Sys,'H',rand(3));
end

idx{1} = 1;
idx{2} = [5 6 7 15];
err = false;
for k = 1:numel(idx)
  Sys_  = nucspinkeep(Sys,idx{k});
  n = numel(idx{k});
  err = err || size(Sys_.A,1)~=3*n || size(Sys_.A,2)~=3;
  err = err || size(Sys_.Q,1)~=3*n || size(Sys_.Q,2)~=3;
end

data = [];
