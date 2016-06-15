% nucspin  Nuclear spin quantum number 
%
%   I = nucspin(Isotopes)
%
%   Returns the spin quantum number of one or
%   several nuclei.
% 
%   Isotopes is a string specifying the
%   nucleus, e.g. '1H', '13C', '63Cu', '191Ir'.
%   If Isotopes is a comma-separated list of nuclei
%   like '14N,14N,14N,1H,63Cu', a vector is returned.
%   Attention, case-sensitive!

function I = nucspin(varargin)

if (nargin==0), help(mfilename); return; end

I = nucdata(varargin);

return
