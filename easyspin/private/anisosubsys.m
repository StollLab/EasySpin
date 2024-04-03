% anisosubsys    extract sub system responsible for anisotropy
%
%   SubSys = anisosubsys(Sys,lw);
%   SubSys = anisosubsys(Sys);

function subSys = anisosubsys(Sys,lw)

subSys = validatespinsys(Sys);

if (nargin==1)
  if isfield(Sys,'HStrain')
    lw = Sys.HStrain;
  else
    lw = [0 0 0];
  end
end

nNuclei = subSys.nNuclei;
if (nNuclei==0), return; end

splitting = 2*subSys.I.*max(abs(subSys.A),[],2).';
rmv = (splitting<min(lw));

% Increase HStrain for each nucleus deleted
% to account for the additional broadening
% caused by this nucleus.
splitting(~rmv) = 0;
subSys.HStrain = sqrt(lw.^2 + sum(splitting.^2));

subSys = nucspinrmv(subSys,find(rmv));

subSys.processed = 0;
subSys = validatespinsys(subSys);

end
