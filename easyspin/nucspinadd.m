% nucspinadd  Adds a nuclear spin to a spin system
%
%    NewSys = nucspinadd(Sys,Nuc,A)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q)
%    NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q,QFrame)
%
%    Add the nuclear isotope Nuc (e.g. '14N') to
%    the spin system structure, with the hyperfine
%    values A, the hyperfine tilt angles AFrame, the
%    quadrupole values Q and the quadrupole tilt angles
%    QFrame. Any missing parameter is assumed to be [0 0 0].
%
%    Alternatively, full 3x3 hyperfine and quadrupole
%    matrices can be specified in A and Q, or symmetric matrices
%    as 6 elements [xx yy zz xy xz yz]. AFrame and QFrame must
%    then be empty or zero.
%
%    If the existing and the added tensors are given in different
%    forms, all tensors are converted to the most general form
%    present: principal values < symmetric matrices (6 elements)
%    < full matrices. AFrame and QFrame are then absorbed into the
%    matrices and removed.
%
%    Examples:
%     Sys = struct('S',1/2,'g',[2 2 2.2]);
%     Sys = nucspinadd(Sys,'Cu',[50 50 520]);
%     Sys = nucspinadd(Sys,'14N',[20 0 0; 0 30 0; 0 0 50]);
%     Sys = nucspinadd(Sys,'1H',[3 4 5 0.5 0 0]);

function NewSys = nucspinadd(Sys,Nuc,A,AFrame,Q,QFrame)

if nargin==0, help(mfilename); return; end

if nargin<2
  error('Second input (isotope) is required.');
end
if nargin<3
  error('Third input (hyperfine tensor) is required.');
end
if nargin<4, AFrame = []; end
if nargin<5, Q = []; end
if nargin<6, QFrame = []; end

% Check Sys
if ~isstruct(Sys)
  error('First input argument must be a spin system structure!');
end
if isfield(Sys,'S')
  if numel(Sys.S)>1
    error('nucspinadd does not work if the system contains more than one electron spin.');
  end
end
if isfield(Sys,'nn') && ~isempty(Sys.nn) && any(Sys.nn(:)~=0)
  error('nucspinadd does not work if Sys.nn is present.');
end

% Check Nuc
if ~ischar(Nuc)
  error('Second input (Nuc) must be a character array, such as ''14N''.');
end

% Supplement AFrame and QFrame
if isempty(AFrame), AFrame = [0 0 0]; end
if isempty(QFrame), QFrame = [0 0 0]; end

% Check A, AFrame, Q, QFrame
if ~any(numel(A)==[1 2 3 6 9])
  error('Wrong size of hyperfine tensor (3rd input argument).');
end

if numel(AFrame)~=3
  error('Wrong size of AFrame (4th input argument).');
end

if ~isempty(Q) && ~any(numel(Q)==[1 2 3 6 9])
  error('Wrong size of quadrupole tensor (5th input argument).');
end

if numel(QFrame)~=3
  error('Wrong size of QFrame (6th input argument).');
end

if numel(A)==9 && any(AFrame(:))
  error('A full hyperfine matrix cannot be combined with nonzero AFrame.');
end
if numel(Q)==9 && any(QFrame(:))
  error('A full quadrupole matrix cannot be combined with nonzero QFrame.');
end
if numel(A)==6 && any(AFrame(:))
  error('A symmetric hyperfine matrix cannot be combined with nonzero AFrame.');
end
if numel(Q)==6 && any(QFrame(:))
  error('A symmetric quadrupole matrix cannot be combined with nonzero QFrame.');
end

% Symmetric matrices [xx yy zz xy xz yz] are stored as rows
if numel(A)==6
  if ~isvector(A)
    error('A symmetric hyperfine matrix must be given as a 6-element vector [xx yy zz xy xz yz].');
  end
  A = A(:).';
end
if numel(Q)==6
  if ~isvector(Q)
    error('A symmetric quadrupole matrix must be given as a 6-element vector [xx yy zz xy xz yz].');
  end
  Q = Q(:).';
end

% Determine number of nuclei
if isfield(Sys,'Nucs')
  Nucs = nucstring2list(Sys.Nucs);
  nNuclei = numel(Nucs);
else
  nNuclei = 0;
end

% Initialize output structure
NewSys = Sys;
if ~isfield(NewSys,'AFrame'), NewSys.AFrame = []; end
if ~isfield(NewSys,'Q'), NewSys.Q = zeros(nNuclei,3); end
if ~isfield(NewSys,'QFrame'), NewSys.QFrame = []; end
iNuc = nNuclei + 1;

% Simplest case: no prior nuclei
if nNuclei==0
  NewSys.Nucs = Nuc;
  NewSys.A = A;
  NewSys.AFrame = AFrame;
  NewSys.Q = Q;
  NewSys.QFrame = QFrame;
  NewSys = cleanemptyfields(NewSys);
  return
end

% Append isotope to Nucs field
NewSys.Nucs = [NewSys.Nucs ',' Nuc];

% Append multiplicity
if isfield(NewSys,'n')
  NewSys.n(iNuc) = 1;
end

% Append atom index (e.g. from orca2easyspin); unknown for added nucleus
if isfield(NewSys,'NucsIdx')
  NewSys.NucsIdx(iNuc) = NaN;
end

% Append A and AFrame
NewSys.A = appendtensor(NewSys.A,NewSys.AFrame,A,AFrame,nNuclei,'A');
fullA = size(NewSys.A,1)==3*iNuc;
symA = size(NewSys.A,2)==6;
if fullA || symA
  NewSys.AFrame = [];  % frames are already included in the matrices
else
  NewSys.AFrame(iNuc,:) = AFrame;
end

% Append Q and QFrame
if isfield(Sys,'Q') || any(Q(:)~=0)
  I = quadrupolespins([Nucs {Nuc}]);
  NewSys.Q = appendtensor(NewSys.Q,NewSys.QFrame,Q,QFrame,nNuclei,'Q',I);
  fullQ = size(NewSys.Q,1)==3*iNuc;
  symQ = size(NewSys.Q,2)==6;
  if fullQ || symQ
    NewSys.QFrame = [];  % frames are already included in the matrices
  else
    NewSys.QFrame(iNuc,:) = QFrame;
  end
end

NewSys = cleanemptyfields(NewSys);

end

%-------------------------------------------------------------------------------
function NewSys = cleanemptyfields(Sys)

NewSys = Sys;

irrelevantfield = @(f) isfield(Sys,f) && (isempty(Sys.(f))||all(Sys.(f)(:)==0));

fields = {'AFrame','Q','QFrame'};
for f = 1:numel(fields)
  if irrelevantfield(fields{f})
    NewSys = rmfield(NewSys,fields{f});
  end
end

end

%-------------------------------------------------------------------------------
function Afull = fullifyA(A,AFrame)
switch numel(A)
  case 1, A = A([1 1 1]);
  case 2, A = A([1 1 2]);
  case 3
  otherwise
    error('A must contain 1, 2, or 3 elements.');
end

if isempty(AFrame) || all(AFrame==0)
  Afull = diag(A);
  return
end

R_T2M = erot(AFrame).'; % tensor frame -> molecular frame
Afull = R_T2M*diag(A)*R_T2M.';

end

%-------------------------------------------------------------------------------
function Qfull = fullifyQ(Q,QFrame,I)
switch numel(Q)
  case 1
    eeQqh = Q;
    eta = 0;
    Qpv = eeQqh*qprefactor(I) * [-1+eta, -1-eta, 2];
  case 2
    eeQqh = Q(1);
    eta = Q(2);
    Qpv = eeQqh*qprefactor(I) * [-1+eta, -1-eta, 2];
  case 3
    Qpv = Q;
  otherwise
    error('Q must contain 1, 2, or 3 elements.');
end

if isempty(QFrame) || all(QFrame==0)
  Qfull = diag(Qpv);
  return
end

R_T2M = erot(QFrame).'; % tensor frame -> molecular frame
Qfull = R_T2M*diag(Qpv)*R_T2M.';

end

%-------------------------------------------------------------------------------
% Convert symmetric matrices given as rows [xx yy zz xy xz yz] (n x 6) to
% stacked full 3x3 matrices (3n x 3).
function M = sym2full(V)
n = size(V,1);
M = zeros(3*n,3);
for k = 1:n
  v = V(k,:);
  M(3*k-2:3*k,:) = [v(1) v(4) v(5); v(4) v(2) v(6); v(5) v(6) v(3)];
end
end

%-------------------------------------------------------------------------------
% Convert a symmetric 3x3 matrix to a row [xx yy zz xy xz yz].
function v = full2sym(M)
v = [M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3)];
end

%-------------------------------------------------------------------------------
function Tnew = appendtensor(T0,T0Frame,T,TFrame,nNuclei,AQ,I)

Atensor = AQ=='A';

if ~Atensor && isempty(T)
  T = zeros(1, 3);
end

if size(T0,1)==3*nNuclei
  nT0 = 9;
elseif numel(T0)==nNuclei
  nT0 = 1;
elseif size(T0,2)==2
  nT0 = 2;
elseif size(T0,2)==3
  nT0 = 3;
elseif size(T0,2)==6
  nT0 = 6;
else
  error('Size of existing %s is inconsistent with the number of nuclei.',AQ);
end
fullT0 = nT0==9;
symT0 = nT0==6;

nT = numel(T);
fullT = nT==9;
symT = nT==6;

if Atensor
  one2two = @(T,I)T(:,[1 1]);
  one2three = @(T,I)T(:,[1 1 1]);
  two2three = @(T,I)T(:,[1 1 2]);
  fullify = @(T,TFrame,I)fullifyA(T,TFrame);
  I0 = zeros(1,nNuclei);  % not used for A
  Inew = 0;
else
  one2two = @(T,I)[T(:) zeros(size(T(:)))];
  one2three = @(T,I)T(:).*qprefactor(I(:)) .* [-1 -1 2];
  two2three = @(T,I)T(:,1).*qprefactor(I(:)) .* [-1+T(:,2) -1-T(:,2) 2*ones(size(T,1),1)];
  fullify = @(T,TFrame,I)fullifyQ(T,TFrame,I);
  I0 = I(1:nNuclei);  % spins of existing nuclei
  Inew = I(end);  % spin of added nucleus
end

Tnew = T0;
newNuc = nNuclei+1;

% Full matrices of existing nuclei given as principal values and Euler angles
if isempty(T0Frame)
  T0Frame = zeros(nNuclei,3);
end
if nT0==1
  T0pv = T0(:);  % isotropic values can be given as row or column
else
  T0pv = T0;
end
pv2full = @(k)fullify(T0pv(k,:),T0Frame(k,:),I0(k));

% Append T, using the most general representation present:
% principal values < symmetric matrices [xx yy zz xy xz yz] < full matrices
if fullT0 || fullT
  if fullT0
    T0full = T0;
  elseif symT0
    T0full = sym2full(T0);
  else
    T0full = zeros(3*nNuclei,3);
    for k = 1:nNuclei
      T0full(3*k-2:3*k,:) = pv2full(k);
    end
  end
  if fullT
    Tfull = T;
  elseif symT
    Tfull = sym2full(T);
  else
    Tfull = fullify(T,TFrame,Inew);
  end
  Tnew = [T0full; Tfull];
elseif symT0 || symT
  if symT0
    T0sym = T0;
  else
    T0sym = zeros(nNuclei,6);
    for k = 1:nNuclei
      T0sym(k,:) = full2sym(pv2full(k));
    end
  end
  if symT
    Tsym = T;
  else
    Tsym = full2sym(fullify(T,TFrame,Inew));
  end
  Tnew = [T0sym; Tsym];
else
  if nT==1
    if nT0==1
      Tnew(newNuc) = T;
    elseif nT0==2
      Tnew(newNuc,:) = one2two(T,Inew);
    else
      Tnew(newNuc,:) = one2three(T,Inew);
    end
  elseif nT==2
    if nT0==1
      Tnew = [one2two(Tnew(:),I0); T];
    elseif nT0==2
      Tnew(newNuc,:) = T;
    else
      Tnew(newNuc,:) = two2three(T,Inew);
    end
  else % nT==3
    if nT0==1
      Tnew = Tnew(:);
      Tnew = [one2three(Tnew,I0); T];
    elseif nT0==2
      Tnew = [two2three(Tnew,I0); T];
    else
      Tnew = [Tnew; T];
    end
  end
end

end

%-------------------------------------------------------------------------------
% Nuclear spins to use for converting quadrupole couplings. For an isotope
% (e.g. '14N'), its spin is used. For an element (e.g. 'N'), the spin of the
% quadrupole reference isotope (most abundant isotope with I>=1) is used,
% or 1/2 if the element has no such isotope.
function I = quadrupolespins(NucList)
I = zeros(1,numel(NucList));
for k = 1:numel(NucList)
  if any(isstrprop(NucList{k},'digit'))
    I(k) = nucspin(NucList{k});
  else
    [~,qref] = referenceisotope(NucList{k});
    if isempty(qref)
      I(k) = 1/2;
    else
      I(k) = qref.I;
    end
  end
end
end

%-------------------------------------------------------------------------------
% Prefactor for converting e^2qQ/h to principal values of Q; zero for I<1,
% since such nuclei have no quadrupole coupling.
function f = qprefactor(I)
f = zeros(size(I));
f(I>=1) = 1./(4*I(I>=1).*(2*I(I>=1)-1));
end
