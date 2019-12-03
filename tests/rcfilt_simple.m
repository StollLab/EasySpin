function [err,data] = test(opt,olddata)

%=======================================================
% Simple signals
%=======================================================
x = linspace(0,23,201);
y{1} = gaussian(x,10,3);
y{2} = zeros(1,1000);

for k = 1:numel(y)
  yy{k} = rcfilt(y{k},1,10);
end

data.yy = yy;

if ~isempty(olddata)
  ok = 1;
  for k = 1:numel(yy)
    ok = ok && areequal(yy{k},olddata.yy{k},1e-10,'rel');
  end
  err = ~ok;
else
  err = [];
end
