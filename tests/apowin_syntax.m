function [err,data] = test(opt,olddata)

Type = {'bla','bar','con','cos','ham','han','wel'};
N = 200;
for k = 1:numel(Type)
  w = apowin(Type{k},N);
  w = apowin([Type{k} '+'],N);
  w = apowin([Type{k} '-'],N);
end

err = 0;
data = [];
