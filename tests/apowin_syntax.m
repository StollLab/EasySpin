function ok = test()

% Test all possible window types

Type = {'bla','bar','con','cos','ham','han','wel'};
N = 200;
for k = 1:numel(Type)
  w = apowin(Type{k},N);
  w = apowin([Type{k} '+'],N);
  w = apowin([Type{k} '-'],N);
end

ok = true;

