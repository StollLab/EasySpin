% Detect sample type (disordered, partially ordered, crystal)
%
%  Required fields in input structures:
%    Exp.Ordering
%    Exp.MolFrame
%    Exp.CrystalSymmetry
%    Exp.SampleFrame
%
%  Output:
%    Exp.MolFrame
%    Exp.CrystalSymmetry
%    Exp.SampleFrame
%    Exp.Ordering
%
%    Opt.partiallyOrderedSample
%    Opt.disorderedSample
%    Opt.crystalSample
%    Opt.GridIntegration

function [Exp,Opt] = p_sampletype(Exp,Opt)

partiallyOrderedSample = ~isempty(Exp.Ordering);
crystalSample = ~partiallyOrderedSample && (~isempty(Exp.MolFrame) || ~isempty(Exp.CrystalSymmetry));
disorderedSample = ~partiallyOrderedSample && ~crystalSample;

% Check fore presence/absence of Exp.MolFrame and Exp.CrystalSymmetry
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

% Supply ExpSampleFrame if absent
if ~disorderedSample && isempty(Exp.SampleFrame)
  Exp.SampleFrame = [0 0 0];
end

Opt.partiallyOrderedSample = partiallyOrderedSample;
Opt.crystalSample = crystalSample;
Opt.disorderedSample = disorderedSample;
Opt.GridIntegration = disorderedSample || partiallyOrderedSample;  % for communication with p_*

% Process Exp.Ordering
if partiallyOrderedSample
  if isnumeric(Exp.Ordering) && numel(Exp.Ordering)==1 && isreal(Exp.Ordering)
    lambda = Exp.Ordering;
    Exp.Ordering = @(beta) exp(lambda*plegendre(2,0,cos(beta)));
    logmsg(1,'  partial ordering (built-in function, coefficient = %g)',lambda);
  elseif isa(Exp.Ordering,'function_handle')
    logmsg(1,'  partial ordering (user-supplied function)');
  else
    error('Exp.Ordering must be either a single number or a function handle.');
  end
  if nargin(Exp.Ordering)==1
    Exp.Ordering = @(beta,gamma) Exp.Ordering(beta).*ones(size(gamma));
  elseif nargin(Exp.Ordering)>2
    logmsg(1,'  Ordering function in Exp.Ordering must take 1 argument (beta) or 2 arguments (beta,gamma).');
  end
end

% Prepare sample rotation matrix (needed to rotate grid in pepper/salt)
if partiallyOrderedSample
  [Exp.R_sample,Opt.rotatedSample] = p_samplerotmatrix(Exp.SampleRotation);
end

end
