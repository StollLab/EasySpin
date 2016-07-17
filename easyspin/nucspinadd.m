% nucspinadd  Adds a nuclear spin to a spin system
%
%    NewSys = nucspinadd(Sys,Nuc,A)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q,QFrame)
%
%    NewSys = nucspinadd(Sys,Nuc,Afull)
%    NewSys = nucspinadd(Sys,Nuc,Afull,Qfull)
%
%    Add the nuclear isotope Nuc (e.g. '14N') to
%    the spin system structure, with the hyperfine
%    values A, the hyperfine tilt angles AFrame, the
%    quadrupole values Q and the quadrupole tilt angles
%    QFrame. Any missing parameter is assumed to be [0 0 0].
%
%    Alternatively, full 3x3 hyperfine and quadrupole
%    matrices can be specified in Afull and Qfull.
%
%    Examples:
%     Sys = struct('S',1/2,'g',[2 2 2.2]);
%     Sys = nucspinadd(Sys,'Cu',[50 50 520]);
%     Sys = nucspinadd(Sys,'14N',[20 0 0; 0 30 0; 0 0 50]);

function NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q,QFrame)

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
  if (nargin<4), Q = zeros(3,3); else, Q = AFrame; end
  if (nargin>4)
    error('If you give a full hyperfine matrix (3rd argument), give Q as fourth argument, and no more.');
  end
else
  if (nargin<4), AFrame = []; end
  if (nargin<5), Q =   []; end
  if (nargin<6), QFrame = []; end
  if isempty(AFrame), AFrame = [0 0 0]; end
  if isempty(QFrame), QFrame = [0 0 0]; end
  if isempty(Q), Q = [0 0 0]; end
end

NewSys = Sys;

if isfield(NewSys,'Nucs')
  Nucs = nucstring2list(NewSys.Nucs);
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
  NewSys.AFrame(newNuc,:) = AFrame;
  NewSys.Q(newNuc,:) = Q;
  NewSys.QFrame(newNuc,:) = QFrame;
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
if isfield(SysFull,'AFrame')
  SysFull = rmfield(SysFull,'AFrame');
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
  if isfield(Sys,'AFrame')
    if ~isempty(Sys.AFrame)
      R_A2M = erot(Sys.AFrame(iNuc,:)).'; % A frame -> molecular frame
      A_ = R_A2M*A_*R_A2M.';
    end
  end
  fullA_ = [fullA_; A_];
end
SysFull.A = fullA_;

return

function SysFull = fullifyQ(Sys,nNuclei)
% fullify quadrupole tensors
SysFull = Sys;

if isfield(SysFull,'QFrame')
  SysFull = rmfield(SysFull,'QFrame');
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
  if isfield(Sys,'QFrame')
    R_Q2M = erot(Sys.QFrame(iNuc,:)).'; % Q frame -> molecular frame
    Q_ = R_Q2M*Q_*R_Q2M.';
  end
  fullQ_ = [fullQ_; Q_];
end
SysFull.Q = fullQ_;

return
