% nucgval  Nuclear g value 
%
%   gn = nucgval(Isotopes)
%
%   Returns the g value of one or several nuclei.
%
%   Isotopes is a string specifying the
%   nucleus, e.g. '1H', '13C', '63Cu', '191Ir'.
%   If Isotopes is a comma-separated list of nuclei
%   like '14N,14N,14N,1H,63Cu', a vector is returned.
%   Attention, case-sensitive!

function gn = nucgval(varargin)

if (nargin==0), help(mfilename); return; end

[ignore,gn] = nucdata(varargin);

return
