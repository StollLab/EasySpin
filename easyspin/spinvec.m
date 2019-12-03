% spinvec  Spin quantum numbers 
%
%   vec = spinvec(System)
%
%   Returns a vector of the quantum
%   numbers of the spins constituting
%   the spin system 'System'.

function vec = spinvec(System)

if nargin==0, help(mfilename); return; end

[Sys,err] = validatespinsys(System);
error(err);

vec = Sys.Spins;

return
