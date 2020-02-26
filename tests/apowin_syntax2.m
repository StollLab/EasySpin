function ok = test()

% Test inputs for windows that need parameters

Type = {'exp','gau','kai'};
alpha = [3 1 5];
N = 340;
for k = 1:numel(Type)
  w = apowin(Type{k},N,alpha(k));
  w = apowin([Type{k} '+'],N,alpha(k));
  w = apowin([Type{k} '-'],N,alpha(k));
end

ok = true;
