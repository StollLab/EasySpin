% nucspinkeep  Remove nuclear spins from spin system
%
%   NewSys = nucspinrmv(Sys,keepidx)
%
%   Removes one or more nuclei from the spin system
%   Sys and returns the result in NewSys. keepidx a vector
%   of nuclear indices as they appear in Sys. E.g. keepidx=1
%   removes all except the first nucleus, and keepidx=[2 4]
%   removes all except the second and the fourth nucleus.
%
%   Example:
%    Sys.Nucs = '13C,1H,14N,14N';
%    Sys.A = [14 5 7 8.2];
%    Sys = nucspinkeep(Sys,2);  % keeps only 1H

function NewSys = nucspinkeep(Sys,keepidx)

if (nargin==0), help(mfilename); return; end

NewSys = Sys;

if (nargin<2), return; end
if isempty(keepidx), return; end

if isfield(Sys,'nn') && ~isempty(Sys.nn) && any(Sys.nn(:)~=0)
  error('nucspinkeep does not work if Sys.nn is given.');
end

Nucs = nucstring2list(NewSys.Nucs);
nNuclei = numel(Nucs);

if any(keepidx>nNuclei) || any(keepidx<=0)
  error('There are only %d nuclei in the spin system. Index out of range.',nNuclei);
end

rmvidx = 1:nNuclei;
rmvidx(keepidx) = [];

NewSys = nucspinrmv(Sys,rmvidx);

return
