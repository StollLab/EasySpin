% nucabund  Natural abundances of nuclear isotopes 
%
%   w = nucabund(Isotopes)
%
%   Returns the natural abundance w of one or
%   several nuclei, as a number between 0
%   (not occuring naturally) and 1 (100%,
%   the only occuring isotope of that element).
%
%   Isotopes is a string specifying the
%   nucleus, e.g. '1H', '13C', '63Cu', '191Ir'.
%   If Isotopes is a comma-separated list of nuclei
%   like '14N,14N,14N,1H,63Cu', a vector is returned.
%   Attention, case-sensitive!


function w = nucabund(varargin)

if (nargin==0), help(mfilename); return; end

[ignore,ignore,ignore,w] = nucdata(varargin);

return
