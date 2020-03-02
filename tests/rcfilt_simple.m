function [ok,data] = test(opt,olddata)

% Simple signals

x = linspace(0,23,201);
y{1} = gaussian(x,10,3);
y{2} = zeros(1,1000);

for k = 1:numel(y)
  yy{k} = rcfilt(y{k},1,10);
end

data.yy = yy;

if ~isempty(olddata)
  for k = 1:numel(yy)
    ok(k) = areequal(yy{k},olddata.yy{k},1e-10,'rel');
  end
else
  ok = [];
end
