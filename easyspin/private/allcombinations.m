function Result = allcombinations(varargin)

if nargin>0 && ischar(varargin{end})
  Elements = varargin(1:end-1);
  mode = varargin{end};
else
  Elements = varargin;
  mode = 'c';
end

switch mode
  case '+'
    fcn = @(a,b)a+b;
    Result = 0;
  case '*'
    fcn = @(a,b)a.*b;
    Result = 1;
  case 'c'
    fcn = @(a,b)[a b];
    Result = [];
  case 'i'
    fcn = @(a,b)[a b];
    Result = [];
    Elements = cellfun(@(x)1:numel(x),Elements,'UniformOutput',false);
end

nGroups = numel(Elements);
nElements = cellfun(@(g)numel(g),Elements);

n1 = prod(nElements);
n2 = 1;
for g = 1:nGroups
  n1 = n1/nElements(g);
  R_ = repmat(Elements{g},n1,n2);
  Result = fcn(Result,R_(:));
  n2 = n2*nElements(g);
end

end
