% Detect sample type (disordered, partially ordered, crystal) and process
% sample orientation.
%
%  Required fields in input structures:
%    Exp.SampleFrame
%    Exp.CrystalSymmetry
%    Exp.MolFrame
%    Exp.Ordering
%
%  Output:
%    Exp.SampleFrame
%    Exp.CrystalSymmetry
%    Exp.MolFrame
%    Exp.Ordering
%
%    Opt.partiallyOrderedSample
%    Opt.disorderedSample
%    Opt.crystalSample
%    Opt.R_L2S
%    Opt.rotatedSample

function [Exp,Opt] = p_sampletype(Exp,Opt)

% Determine sample type
partiallyOrderedSample = ~isempty(Exp.Ordering);
crystalSample = ~partiallyOrderedSample && (~isempty(Exp.MolFrame) || ~isempty(Exp.CrystalSymmetry));
disorderedSample = ~partiallyOrderedSample && ~crystalSample;

% Check Exp.MolFrame and Exp.CrystalSymmetry, supplement if needed
if partiallyOrderedSample
  if ~isempty(Exp.MolFrame)
    error('Exp.MolFrame cannot be used for partially ordered samples (Exp.Ordering given).');
  elseif ~isempty(Exp.CrystalSymmetry)
    error('Exp.CrystalSymmetry cannot be used for partially ordered samples (Exp.Ordering given).');
  end
end
if crystalSample
  if isempty(Exp.CrystalSymmetry)
    Exp.CrystalSymmetry = 'P1';
  end
  if isempty(Exp.MolFrame)
    Exp.MolFrame = [0 0 0];
  end
end

% Supply Exp.SampleFrame if absent
if isempty(Exp.SampleFrame)
  Exp.SampleFrame = [0 0 0];
end
if ~isnumeric(Exp.SampleFrame)
  error('Exp.SampleFrame must be an Nx3 array.');
end

% Check whether sample is rotated
Opt.rotatedSample = any(Exp.SampleFrame(:)) || (isfield(Exp,'SampleRotation') && ~isempty(Exp.SampleRotation));

% Process Exp.Ordering
if partiallyOrderedSample
  if isnumeric(Exp.Ordering) && numel(Exp.Ordering)==1 && isreal(Exp.Ordering)
    lambda = Exp.Ordering;
    Exp.Ordering = @(alpha,beta,gamma) exp(lambda*plegendre(2,0,cos(beta)));
    logmsg(1,'  partial ordering (built-in function, coefficient = %g)',lambda);
  elseif isa(Exp.Ordering,'function_handle')
    logmsg(1,'  partial ordering (user-supplied function)');
  else
    error('Exp.Ordering must be either a single number or a function handle.');
  end
  if nargin(Exp.Ordering)==1
    Exp.Ordering = @(alpha,beta,gamma) Exp.Ordering(beta).*ones(size(gamma));
  elseif nargin(Exp.Ordering)~=3
    logmsg(1,'  Ordering function in Exp.Ordering must take 1 input argument (beta) or 3 input arguments (alpha,beta,gamma).');
  end
end

% Calculate sample orientation
[Opt.R_L2S,Opt.rotatedSample] = p_samplerotation(Exp);

if ~crystalSample && numel(Opt.R_L2S)>1
  error('For disordered and partially ordered samples, only a single sample orientation (Exp.SampleFrame) can be used.');
end


% Collect output
Opt.partiallyOrderedSample = partiallyOrderedSample;
Opt.crystalSample = crystalSample;
Opt.disorderedSample = disorderedSample;

end
