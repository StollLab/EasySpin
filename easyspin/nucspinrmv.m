% nucspinrmv   Remove nuclear spin from spin system
%
%   NewSys = nucspinrmv(Sys,rmvidx)
%
%   Removes one or more nuclei from the spin system
%   Sys and returns the result in NewSys. rmvidx a vector
%   of nuclear numbers as they appear in Sys. E.g. rmvidx=1
%   removes the first nucleus, rmvidx=[2 4] removes the
%   second and the fourth.
%
%   Example:
%    Sys = struct('S',1/2,'g',[2 2 3],'Nucs','14N,1H','A',[4 6 10; 1 1 2]);
%    Sys = nucspinrmv(Sys,1);  % removes the 1H

function NewSys = nucspinrmv(Sys,rmvidx)

if (nargin==0), help(mfilename); return; end

NewSys = Sys;

if (nargin<2), return; end
if isempty(rmvidx), return; end

Nucs = nucstring2list(NewSys.Nucs);
nNuclei = numel(Nucs);

if any(rmvidx>nNuclei) || any(rmvidx<=0)
  error('There are only %d nuclei in the spin system. Index out of range.',nNuclei);
end

% nuclei -----------------------------------------------
Nucs(rmvidx) = [];
Nucs = nuclist2string(Nucs);
if isempty(Nucs)
  NewSys = rmfield(NewSys,'Nucs');
else
  NewSys.Nucs = Nucs;
end

% multiplicities ---------------------------------------
if isfield(NewSys,'n')
  n = NewSys.n;
  n(rmvidx) = [];
  if isempty(Nucs)
    NewSys = rmfield(NewSys,'n');
  else
    NewSys.n = n;
  end
end


% hyperfine ---------------------------------------
if ~isfield(NewSys,'fullA')
  fullA = size(NewSys.A,1)==3*nNuclei;
else
  fullA = NewSys.fullA;
end

if ~fullA
  F = {'A','AFrame'};
  for k = 1:2
    Field = F{k};
    if isfield(NewSys,Field)
      v = NewSys.(Field);
      v(rmvidx,:) = [];
      if isempty(v)
        NewSys = rmfield(NewSys,Field);
      else
        NewSys.(Field) = v;
      end
    end
  end
else
  v = NewSys.A;
  for iNuc = numel(rmvidx):-1:1
    v((rmvidx(iNuc)-1)*3+(1:3),:) = [];
  end
  if isempty(v)
    NewSys = rmfield(NewSys,'A');
  else
    NewSys.A = v;
  end
  if isfield(NewSys,'AFrame')
    NewSys = rmfield(NewSys,'AFrame');
  end
end

% nuclear quadrupole ----------------------------------
if ~isfield(NewSys,'fullQ')
  if isfield(NewSys,'Q')
    fullQ = size(NewSys.Q,1)==3*nNuclei;
  else
    fullQ = false;
  end
else
  fullQ = NewSys.fullQ;
end

if ~fullQ
  F = {'Q','QFrame'};
  for k = 1:2
    Field = F{k};
    if isfield(NewSys,Field)
      v = NewSys.(Field);
      v(rmvidx,:) = [];
      if isempty(v)
        NewSys = rmfield(NewSys,Field);
      else
        NewSys.(Field) = v;
      end
    end
  end
else
  v = NewSys.Q;
  for iNuc = numel(rmvidx):-1:1
    v((rmvidx(iNuc)-1)*3+(1:3),:) = [];
  end
  if isempty(v)
    NewSys = rmfield(NewSys,'Q');
  else
    NewSys.Q = v;
  end
  if isfield(NewSys,'AFrame')
    NewSys = rmfield(NewSys,'AFrame');
  end
end

return
