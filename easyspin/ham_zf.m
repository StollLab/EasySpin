% ham_zf  Electronic zero field interaction Hamiltonian 
%
%   F = ham_zf(SpinSystem)
%   F = ham_zf(SpinSystem,Electrons)
%   F = ham_zf(SpinSystem,Electrons,'sparse')
%
%   Returns the electronic zero-field interaction (ZFI)
%   Hamiltonian of the system SpinSystem, in units of MHz.
%
%   If the vector Electrons is given, the ZFI of only the
%   specified electrons is returned (1 is the first, 2 the
%   second, etc). Otherwise, all electrons are included.
%
%   If 'sparse' is given, the matrix is returned in sparse format.

function H = ham_zf(SpinSystem,idxElectrons,opt)

if nargin==0, help(mfilename); return; end

[Sys,err] = validatespinsys(SpinSystem);
error(err);

if nargin<2, idxElectrons = []; end
if nargin<3, opt = ''; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

if isempty(idxElectrons)
  idxElectrons = 1:Sys.nElectrons;
end

if any(idxElectrons>Sys.nElectrons) || any(idxElectrons<1)
  error('Electron spin index/indices (2nd argument) out of range!');
end

H = sparse(Sys.nStates,Sys.nStates);
Spins = Sys.Spins;

% Quadratic term S*D*S
%-------------------------------------------------------------------------------
for iSpin = idxElectrons
    
  % Prepare full 3x3 D matrix
  if Sys.fullD
    D = Sys.D(3*(iSpin-1)+(1:3),:);
  else
    D = diag(Sys.D(iSpin,:));
  end

  if ~any(D(:))
    continue
  end

  % Transform D tensor to molecular frame
  ang = Sys.DFrame(iSpin,:);
  if any(ang)
    R_M2D = erot(ang);  % mol frame -> D frame
    R_D2M = R_M2D.';    % D frame -> mol frame
    D = R_D2M*D*R_D2M.';
  end

  % Construct spin operator matrices
  for c = 3:-1:1
    Sxyz{c} = sop(Spins,[iSpin,c],'sparse');
  end

  % Construct SDS term
  for c1 = 1:3
    for c2 = 1:3
      H = H + D(c1,c2)*(Sxyz{c1}*Sxyz{c2});
    end
  end

end

% Fourth-order terms a and F
%-----------------------------------------------------------------------------
% Spin system fields: a and F

% Abragam/Bleaney, pp. 142, 437
% Bleaney/Trenam, Proc.Roy.Soc.A, 223(1152), 1-14, (1954)
% Doetschman/McCool, Chem.Phys. 8, 1-16 (1975)
% Scullane, J.Magn.Reson. 47, 383-397 (1982)
% Jain/Lehmann, phys.stat.sol.(b) 159, 495-544 (1990)

% These terms are no longer supported.
if isfield(Sys,'aF')
  error('Sys.aF is no longer supported. Use Sys.B4 and Sys.B4Frame instead.');
end

% High-order terms in extended Stevens operator format
%-----------------------------------------------------------------------------
% Spin system fields:
%   B1, B2, B3, ... (k = 1...12) -> processed to Sys.B
%   B1Frame, B2Frame, B3Frame, ... -> processed to Sys.BFrame

% Run over all ranks k
for k = 1:numel(Sys.B)
  Bk = Sys.B{k};
  if isempty(Bk), continue; end

  % Run over all desired electron spins
  for iSpin = idxElectrons

    % If D is used, skip rank-2 Stevens operator terms
    D_present = isfield(Sys,'D') && ~isempty(D) && any(D(:)~=0);
    if D_present && k==2, continue; end

    % Skip if all rank-k coefficients are zero
    if all(Bk(iSpin,:)==0), continue; end
    
    % Apply transformation if non-zero tilt angles are given
    % (Sys.BFrame is processed from Sys.B?Frame by validatespinsys)
    tiltang = Sys.BFrame{k}(iSpin,:);
    if any(tiltang)
      % Calculate transformation matrix for ISTOs
      Dk = wignerd(k,tiltang(1),tiltang(2),tiltang(3));
      Ck = isto2stev(k);
      % Transform to Stevens operator basis
      DB = Ck*Dk.'/Ck;
      DB = removenumericalnoise(DB);
      % Transform Bk
      Bk_M = Bk(iSpin,:)*DB;  % eigenframe of Bk -> molecular frame
    else
      Bk_M = Bk(iSpin,:);
    end
    
    % Build Hamiltonian
    q = k:-1:-k;
    for iq = find(Bk_M~=0)
      H = H + Bk_M(iq)*stev(Spins,[k,q(iq),iSpin],'sparse');
    end
    
  end % for all electron spins specified

end % for all tensor ranks

H = (H+H')/2;  % Hermitianize
if ~sparseResult
  H = full(H);
end

end
%===============================================================================

function A_ = removenumericalnoise(A)
reA = real(A);
imA = imag(A);
thr = 1e-14;
reA(abs(reA)<thr&reA~=0) = 0;
imA(abs(imA)<thr&imA~=0) = 0;
A_ = complex(reA,imA);
end
