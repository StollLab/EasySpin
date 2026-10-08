function ok = test()

% Analytical tensor derivatives of strain parameters against finite
% differences of the molecular-frame tensors

ang = [0.3 0.7 1.1];
c = 0;

% g: isotropic, axial, principal values, symmetric, frame angles
c=c+1; S{c} = struct('g',2.1,'gFrame',ang); P{c} = {'g'};
c=c+1; S{c} = struct('g',[2.0 2.2],'gFrame',ang); P{c} = {'g(1)','g(2)','gFrame(1)','gFrame(2)'}; % gamma has no effect on an axial tensor
c=c+1; S{c} = struct('g',[2.0 2.1 2.2],'gFrame',ang); P{c} = {'g(2)','gFrame(1)','gFrame(2)','gFrame(3)'};
c=c+1; S{c} = struct('g',[2.0 2.1 2.2 0.01 0.02 0.03]); P{c} = {'g(1)','g(4)','g(5)','g(6)'};
% D: D, [D E], principal values, frame angles
c=c+1; S{c} = struct('S',1,'D',300,'DFrame',ang); P{c} = {'D'};
c=c+1; S{c} = struct('S',1,'D',[300 50],'DFrame',ang); P{c} = {'D(1)','D(2)','DFrame(1)','DFrame(2)','DFrame(3)'};
c=c+1; S{c} = struct('S',1,'D',[-100 -50 150],'DFrame',ang); P{c} = {'D(1)','D(3)'};
% A: iso (1 x nNuc), axial, principal values, symmetric, frames
c=c+1; S{c} = struct('Nucs','1H,1H','A',[5 8]); P{c} = {'A(2)'};
c=c+1; S{c} = struct('Nucs','1H','A',[5 8],'AFrame',ang); P{c} = {'A(1)','A(2)','AFrame(2)'};
c=c+1; S{c} = struct('Nucs','1H','A',[5 8 12],'AFrame',ang); P{c} = {'A(3)','AFrame(1)','AFrame(3)'};
c=c+1; S{c} = struct('Nucs','1H','A',[5 8 12 1 2 3]); P{c} = {'A(5)'};
% Q: eeqQ, [eeqQ eta], principal values, frames
c=c+1; S{c} = struct('Nucs','14N','A',3,'Q',2,'QFrame',ang); P{c} = {'Q','QFrame(2)'};
c=c+1; S{c} = struct('Nucs','14N','A',3,'Q',[2 0.3],'QFrame',ang); P{c} = {'Q(1)','Q(2)','QFrame(1)'};
c=c+1; S{c} = struct('Nucs','14N','A',3,'Q',[-0.5 -0.3 0.8],'QFrame',ang); P{c} = {'Q(2)'};
% sigma
c=c+1; S{c} = struct('Nucs','1H','A',3,'sigma',[1 1.001 1.002],'sigmaFrame',ang); P{c} = {'sigma(2)','sigmaFrame(2)'};
% ee: iso, principal values, frames; J, dip (1, 2, 3 columns), dvec
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2],'ee',100,'eeFrame',ang); P{c} = {'ee'};
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2],'ee',[90 100 120],'eeFrame',ang); P{c} = {'ee(3)','eeFrame(2)'};
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2],'J',100,'dip',20,'eeFrame',ang); P{c} = {'J','dip','eeFrame(1)'};
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2],'J',100,'dip',[20 5]); P{c} = {'dip(1)','dip(2)'};
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2],'J',100,'dip',[-10 -5 15]); P{c} = {'dip(2)'};
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2],'J',100,'dvec',[1 2 3],'eeFrame',ang); P{c} = {'dvec(1)','dvec(2)','dvec(3)','eeFrame(2)'};
c=c+1; S{c} = struct('S',[1/2 1/2],'g',[2 2],'J',100,'dvec',[0 0 0]); P{c} = {'dvec(2)'};
% nn
c=c+1; S{c} = struct('Nucs','1H,1H','A',[5 8],'nn',[1 2 3],'nnFrame',ang); P{c} = {'nn(2)','nnFrame(3)'};

h = 1e-6;
ok = true(1,c);
for k = 1:c
  Sys = S{k};
  Sys.StrainPars = P{k};
  Sys.StrainFWHM = ones(1,numel(P{k}));
  [Sys_,err] = runprivate('validatespinsys',Sys);
  if ~isempty(err), ok(k) = false; continue; end
  Sys = rmfield(Sys,{'StrainPars','StrainFWHM'});
  for p = 1:numel(P{k})
    D = Sys_.StrainData.Deriv(p);
    [f,idx] = parseref(P{k}{p});
    Sp = Sys; Sp.(f)(idx) = Sp.(f)(idx) + h;
    Sm = Sys; Sm.(f)(idx) = Sm.(f)(idx) - h;
    dT = (moltensor(Sp,D.type,D.idx) - moltensor(Sm,D.type,D.idx))/(2*h);
    ok(k) = ok(k) && areequal(D.dT,dT,1e-6*max(1,norm(dT)),'abs');
  end
end

end

%-------------------------------------------------------------------------------
function [f,idx] = parseref(str)
tok = regexp(str,'^(\w+)(\((\d+)\))?','tokens','once');
f = tok{1};
idx = 1;
if numel(tok)>1 && ~isempty(tok{2}), idx = str2double(tok{2}(2:end-1)); end
end

%-------------------------------------------------------------------------------
% Interaction tensor in the molecular frame, as used by the ham_* functions
function T = moltensor(Sys,type,idx)
Sys = runprivate('validatespinsys',Sys);
switch type
  case {'g','D'}
    e = idx; full_ = Sys.(['full' type]); rows = 3*(e-1)+(1:3);
    if full_, T = Sys.(type)(rows,:); else, T = diag(Sys.(type)(e,:)); end
    ang = Sys.([type 'Frame'])(e,:);
  case {'Q','sigma'}
    n = idx; full_ = Sys.(['full' type]); rows = 3*(n-1)+(1:3);
    if full_, T = Sys.(type)(rows,:); else, T = diag(Sys.(type)(n,:)); end
    ang = Sys.([type 'Frame'])(n,:);
  case 'A'
    e = idx(1); n = idx(2); cols = 3*(e-1)+(1:3);
    if Sys.fullA, T = Sys.A(3*(n-1)+(1:3),cols); else, T = diag(Sys.A(n,cols)); end
    ang = Sys.AFrame(n,cols);
  case {'ee','nn'}
    if strcmp(type,'ee'), N = Sys.nElectrons; else, N = Sys.nNuclei; end
    pairs = nchoosek(1:N,2);
    p = find(pairs(:,1)==idx(1) & pairs(:,2)==idx(2));
    if Sys.(['full' type]), T = Sys.(type)(3*(p-1)+(1:3),:); else, T = diag(Sys.(type)(p,:)); end
    ang = Sys.([type 'Frame'])(p,:);
end
R = erot(ang).';
T = R*T*R.';
end
