% nucspinadd  Adds a nuclear spin to a spin system
%
%    NewSys = nucspinadd(Sys,Nuc,A)
%    NewSys = nucspinadd(Sys,Nuc,A,Apa)
%    NewSys = nucspinadd(Sys,Nuc,A,Apa,Q)
%    NewSys = nucspinadd(Sys,Nuc,A,Apa,Q,Qpa)
%
%    NewSys = nucspinadd(Sys,Nuc,Afull)
%    NewSys = nucspinadd(Sys,Nuc,Afull,Qfull)
%
%    Add the nuclear isotope Nuc (e.g. '14N') to
%    the spin system structure, with the hyperfine
%    values A, the hyperfine tilt angles Apa, the
%    quadrupole values Q and the quadrupole tilt angles
%    Qpa. Any missing parameter is assumed to be [0 0 0].
%
%    Alternatively, full 3x3 hyperfine and quadrupole
%    matrices can be specified in Afull and Qfull.
%
%    Examples:
%     Sys = struct('S',1/2,'g',[2 2 2.2]);
%     Sys = nucspinadd(Sys,'Cu',[50 50 520]);
%     Sys = nucspinadd(Sys,'14N',[20 0 0; 0 30 0; 0 0 50]);

function NewSys = nucspinadd(Sys,Nuc,A,Apa,Q,Qpa)

if (nargin==0), help(mfilename); return; end

if ~isstruct(Sys)
  error('First input argument must be a spin system structure!');
end

if isfield(Sys,'S')
  if numel(Sys.S)>1
    error('nucspinadd works only for 1 electron spin.');
  end
end

if (nargin<3)
  error('At least the hyperfine values are needed!');
end

if ~any(numel(A)==[1 2 3 9])
  error('Wrong size of hyperfine tensor (3rd input argument).');
end

fullA = numel(A)==9;

if fullA
  if (nargin<4), Q = zeros(3,3); else, Q = Apa; end
  if (nargin>4)
    error('If you give a full hyperfine matrix (3rd argument), give Q as fourth argument, and no more.');
  end
else
  if (nargin<4), Apa = []; end
  if (nargin<5), Q =   []; end
  if (nargin<6), Qpa = []; end
  if isempty(Apa), Apa = [0 0 0]; end
  if isempty(Qpa), Qpa = [0 0 0]; end
  if isempty(Q), Q = [0 0 0]; end
end

NewSys = Sys;

if isfield(NewSys,'Nucs')
  Nucs = nucstringparse(NewSys.Nucs);
  nNuclei = numel(Nucs);
else
  nNuclei = 0;
end
newNuc = nNuclei + 1;

if (nNuclei==0)
  NewSys.Nucs = Nuc;
else
  NewSys.Nucs = [NewSys.Nucs ',' Nuc];
end

if isfield(NewSys,'n')
  NewSys.n(newNuc) = 1;
end

if fullA
  if nNuclei==0
    NewSys.A = A;
    NewSys.Q = Q;
  else
    NewSys = fullifyA(NewSys,nNuclei);
    NewSys = fullifyQ(NewSys,nNuclei);
    NewSys.A = [NewSys.A; A];
    NewSys.Q = [NewSys.Q; Q];
  end
else
  NewSys.A(newNuc,:) = tensorexpand(A);
  NewSys.Apa(newNuc,:) = Apa;
  NewSys.Q(newNuc,:) = Q;
  NewSys.Qpa(newNuc,:) = Qpa;
end

return

function T_ = tensorexpand(T)
T_ = T;
switch numel(T)
  case 1, T_ = T([1 1 1]);
  case 2, T_ = T([1 1 2]);
end

function SysFull = fullifyA(Sys,nNuclei)
% fullify hyperfine tensors
SysFull = Sys;
if isfield(SysFull,'Apa')
  SysFull = rmfield(SysFull,'Apa');
end
if ~isfield(SysFull,'A')
  SysFull.A = zeros(nNuclei*3,3);
  return;
end
if numel(Sys.A)==9*nNuclei, return; end

fullA_ = [];
for iNuc=1:nNuclei
  if numel(Sys.A)==nNuclei
    A_ = [1 1 1]*Sys.A(iNuc);
  elseif numel(Sys.A)==2*nNuclei
    A_ = Sys.A(iNuc,[1 1 2]);
  else
    A_ = Sys.A(iNuc,:);
  end
  A_ = diag(A_);
  if isfield(Sys,'Apa')
    if ~isempty(Sys.Apa)
      R_ = erot(Sys.Apa(iNuc,:));
      A_ = R_*A_*R_';
    end
  end
  fullA_ = [fullA_; A_];
end
SysFull.A = fullA_;

return

function SysFull = fullifyQ(Sys,nNuclei)
% fullify quadrupole tensors
SysFull = Sys;

if isfield(SysFull,'Qpa')
  SysFull = rmfield(SysFull,'Qpa');
end

if ~isfield(SysFull,'Q')
  SysFull.Q = zeros(nNuclei*3,3);
  return;
end

if numel(Sys.Q)==9*nNuclei, return; end

if numel(Sys.Q)~=nNuclei*3
  error('nucspinadd works only if Sys.Q has 3 or 9 values per nucleus.');
end

fullQ_ = [];
for iNuc = 1:nNuclei
  Q_ = Sys.Q(iNuc,:);
  Q_ = diag(Q_);
  if isfield(Sys,'Qpa')
    R_ = erot(Sys.Qpa(iNuc,:));
    Q_ = R_*Q_*R_';
  end
  fullQ_ = [fullQ_; Q_];
end
SysFull.Q = fullQ_;

return
