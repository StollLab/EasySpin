% nucspinrmv   Remove nuclear spin from spin system
%
%   NewSys = nucspinrmv(Sys,idx)
%
%   Removes one or more nuclei from the spin system
%   Sys and returns the result in NewSys. idx a vector
%   of nuclear numbers as they appear in Sys. E.g. idx=1
%   removes the first nucleus, idx=[2 4] removes the
%   second and the fourth.
%
%   Example:
%    Sys = struct('S',1/2,'g',[2 2 3],'Nucs','14N','A',[4 6 10]);
%    Sys = nucspinrmv(Sys,1);

function NewSys = nucspinrmv(Sys,idx)

if (nargin==0), help(mfilename); return; end

NewSys = Sys;

if (nargin<2), return; end
if isempty(idx), return; end

Nucs = nucstring2list(NewSys.Nucs);
nNuclei = numel(Nucs);

if any(idx>nNuclei) || any(idx<=0)
  error('There are only %d nuclei in the spin system. Index out of range.',nNuclei);
end

Nucs(idx) = [];

Nucs = nuclist2string(Nucs);

if isempty(Nucs)
  NewSys = rmfield(NewSys,'Nucs');
else
  NewSys.Nucs = Nucs;
end

% multiplicities ---------------------------------------
if isfield(NewSys,'n')
  n = NewSys.n;
  n(idx) = [];
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
  for k=1:2
    Field = F{k};
    if isfield(NewSys,Field)
      v = NewSys.(Field);
      v(idx,:) = [];
      if isempty(v)
        NewSys = rmfield(NewSys,Field);
      else
        NewSys.(Field) = v;
      end
    end
  end
else
  v = NewSys.A;
  for iNuc=1:numel(idx)
    v((idx(iNuc)-1)*3+(1:3),:) = NaN;
  end
  v(isnan(v))=[];
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
F = {'Q','QFrame'};
for k=1:2
  Field = F{k};
  if isfield(NewSys,Field)
    v = NewSys.(Field);
    v(idx,:) = [];
    if isempty(v)
      NewSys = rmfield(NewSys,Field);
    else
      NewSys.(Field) = v;
    end
  end
end

return
