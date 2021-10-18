% nucabund  Nuclear electric quadrupole moments 
%
%   w = nucqmom(Isotopes)
%
%   Returns the nuclear electric quadrupole
%   moment of one or several nuclei, in barn
%   (1 barn = 10^-28 m^2 = 10^-24 cm^2).
%
%   Isotopes is a string specifying the
%   nucleus, e.g. '1H', '13C', '63Cu', '191Ir'.
%   If Isotopes is a comma-separated list of nuclei
%   like '14N,14N,14N,1H,63Cu', a vector is returned.
%   Attention, case-sensitive!

function Q = nucqmom(varargin)

if (nargin==0), help(mfilename); return; end

[ignore,ignore,Q] = nucdata(varargin);

return
