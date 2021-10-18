% larmorfrq   Larmor frequency of nuclear spins 
%
%   Freqs = larmorfrq(Nucs,Fields)
%
%   Computes Larmor frequencies (in MHz) for the
%   nuclear spins given in Nucs at magnetic fields
%   (in mT) given in Fields.
%
%   Nucs is a string containing a comma-separated
%   list of nuclear isotope symbols, e.g. '14N'
%   or '63Cu,65Cu'.
%
%   Fields is a vector of magnetic field values in mT,
%   e.g. 350 or 300:10:380.
%
%   Each column in Freqs contains the Larmor frequencies
%   in MHz for one of the nuclei listed in Nucs.
%   
%   Example:
%     nu = larmorfrq('1H',3500)
%     nu = larmorfrq('1H,14N',300:10:350)

function LarmFreq = larmorfrq(Nucs,Fields)

if nargin==0
  help(mfilename);
  return
end

if nargin<2 || nargin>2, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>1, error('Too many output arguments.'); end

if isempty(Nucs) || isempty(Fields)
  LarmFreq = [];
  return
end

if iscell(Nucs)
  Nucs = nuclist2string(Nucs);
end

gn = abs(nucgval(Nucs));

prefac = gn*nmagn/1e3/planck/1e6;
for iF = 1:numel(Fields)
  LarmFreq(iF,:) = prefac*Fields(iF);
end

return
