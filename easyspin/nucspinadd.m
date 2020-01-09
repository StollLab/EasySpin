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
%    matrices can be specified in Afull and Qfull.
%
%    Examples:
%     Sys = struct('S',1/2,'g',[2 2 2.2]);
%     Sys = nucspinadd(Sys,'Cu',[50 50 520]);
%     Sys = nucspinadd(Sys,'14N',[20 0 0; 0 30 0; 0 0 50]);

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
  errpor('nucspinadd does not work if Sys.nn is present.');
end

% Check Nuc
if ~ischar(Nuc)
  error('Second input (Nuc) must be a character array, such as ''14N''.');
end

% Supplement AFrame and QFrame
if isempty(AFrame), AFrame = [0 0 0]; end
if isempty(QFrame), QFrame = [0 0 0]; end

% Check A, AFrame, Q, QFrame
if ~any(numel(A)==[1 2 3 9])
  error('Wrong size of hyperfine tensor (3rd input argument).');
end

if numel(AFrame)~=3
  error('Wrong size of AFrame (4th input argument).');
end

if ~isempty(Q) && ~any(numel(Q)==[1 2 3 9])
  error('Wrong size of quadrupole tensor (5th input argument).');
end

if numel(QFrame)~=3
  error('Wrong size of QFrame (6th input argument).');
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

% Append A and AFrame
NewSys.A = appendtensor(NewSys.A,NewSys.AFrame,A,AFrame,nNuclei,'A');
fullA = numel(A)==9 || size(Sys.A,1)==3*nNuclei;
if ~fullA
  NewSys.AFrame(iNuc,:) = AFrame;
end

% Append Q and QFrame
if isfield(Sys,'Q') || any(Q(:)~=0)
  I = nucspin(Nucs);
  NewSys.Q = appendtensor(NewSys.Q,NewSys.QFrame,Q,QFrame,nNuclei,'Q',I);
  fullQ = numel(Q)==9 || size(NewSys.Q,1)==3*nNuclei;
  if ~fullQ
    NewSys.QFrame(iNuc,:) = QFrame;
  end
end

NewSys = cleanemptyfields(NewSys);

return

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

return

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

R_T2M = erot(AFrame); % tensor frame -> molecular frame
Afull = R_T2M*diag(A)*R_T2M.';

return

%-------------------------------------------------------------------------------
function Qfull = fullifyQ(Q,QFrame,I)
switch numel(Q)
  case 1
    eeQqh = Q;
    eta = 0;
    Qpv = eeQqh/(4*I*(2*I-1)) * [-1+eta, -1-eta, 2];
  case 2
    eeQqh = Q(1);
    eta = Q(2);
    Qpv = eeQqh/(4*I*(2*I-1)) * [-1+eta, -1-eta, 2];
  case 3
    Qpv = Q;
  otherwise
    error('Q must contain 1, 2, or 3 elements.');
end

if isempty(QFrame) || all(QFrame==0)
  Qfull = diag(Qpv);
  return
end

R_T2M = erot(QFrame); % tensor frame -> molecular frame
Qfull = R_T2M*diag(Qpv)*R_T2M.';

return

%-------------------------------------------------------------------------------
function Tnew = appendtensor(T0,T0Frame,T,TFrame,nNuclei,AQ,I)

Atensor = AQ=='A';

if size(T0,1)==3*nNuclei
  nT0 = 9;
elseif numel(T0)==nNuclei
  nT0 = 1;
elseif size(T0,2)==2
  nT0 = 2;
elseif size(T0,2)==3
  nT0 = 3;
end
fullT0 = nT0==9;

nT = numel(T);
fullT = nT==9;

if Atensor
  one2two = @(T)T(:,[1 1]);
  one2three = @(T)T(:,[1 1 1]);
  two2three = @(T)T(:,[1 1 2]);
else
  one2two = @(T)[T(:) zeros(size(T(:)))];
  one2three = @(T)T(:)./(4*I.*(2*I-1)) * [-1 -1 2];
  two2three = @(T)T(:,1)./(4*I.*(2*I-1)) * [-1+T(:,2) -1-T(:,2) 2];
end

Tnew = T0;
newNuc = nNuclei+1;

% Append T
if fullT0 || fullT
  if ~fullT0
    if Atensor
      Tnew = [fullifyA(T0,T0Frame); T];
    else
      Tnew = [fullifyQ(T0,T0Frame,I); T];
    end
  elseif ~fullT
    if Atensor
      Tnew = [T0; fullifyA(T,TFrame)];
    else
      Tnew = [T0; fullifyQ(T,TFrame,I)];
    end
  else
    Tnew = [T0; T];
  end
else
  if nT==1
    if nT0==1
      Tnew(newNuc) = T;
    elseif nT0==2
      Tnew(newNuc,:) = one2two(T);
    else
      Tnew(newNuc,:) = one2three(T);
    end
  elseif nT==2
    if nT0==1
      Tnew = [one2two(Tnew(:)); T];
    elseif nT0==2
      Tnew(newNuc,:) = T;
    else
      Tnew(newNuc,:) = two2three(T);
    end
  else % nT==3
    if nT0==1
      Tnew = Tnew(:);
      Tnew = [one2three(Tnew); T];
    elseif nT0==2
      Tnew = [two2three(Tnew); T];
    else
      Tnew = [Tnew; T];
    end
  end
end
return
