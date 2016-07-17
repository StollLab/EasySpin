function [err,data] = test(opt,olddata)

% Test 3: bin a lot of peaks
%=======================================================
err = 1;
N = 20e3;
pos = rand(1,N); amp = rand(1,N);
[x,s] = makespec([-1,2],2000,pos,amp);
if abs(sum(amp)-sum(s))<N*1e-10
  err = 0;
end
data = [];
