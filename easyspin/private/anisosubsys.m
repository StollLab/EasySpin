% anisosubsys    extract sub system responsible for anisotropy
%
%   SubSys = anisosubsys(Sys,lw);
%   SubSys = anisosubsys(Sys);

function subSys = anisosubsys(Sys,lw);

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

return
%====================================================================

%-------------------------------
% test area

clear

Field = 350;
Nuc = '35Cl';
nuI = larmorfrq(Nuc,Field)
alpha = linspace(0,4,101);
mw = bmagn*Field*1e-3*2/planck/1e6;
Sys = struct('S',1/2,'g',[2 2 2],'Nucs',Nuc);
for k=1:numel(alpha)
  Sys.A = [1 1 1]*nuI*alpha(k);
  H = sham(Sys,[0 0 Field]); E = eig(H);
  dE = alldiff(E); dE = dE(dE>mw/2);
  freqs(:,k) =  dE-mw;
  spread(k) = max(dE)-min(dE);
end

plot(alpha,spread/nuI);

%-------------------------------

Sys = struct('S',1/2,'g',[2 2 2],'Nucs','1H,1H','A',[[1 1 1]*30;[1 1 1]*5]);
Sys.HStrain = [1 1 1]*30;

anisosubsys(Sys)
