% zfield  Electronic zero field interaction Hamiltonian 
%
%   F = zfield(SpinSystem)
%   F = zfield(SpinSystem,Electrons)
%   F = zfield(SpinSystem,Electrons,'sparse')
%
%   Returns the electronic zero-field interaction (ZFI)
%   Hamiltonian [MHz] of the system SpinSystem.
%
%   If the vector Electrons is given, the ZFI of only the
%   specified electrons is returned (1 is the first, 2 the
%   second, etc). Otherwise, all electrons are included.
%
%   If 'sparse' is given, the matrix is returned in sparse format.

function H = zfield(SpinSystem,idxElectrons,opt)

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

for iSpin = idxElectrons
    
  % Quadratic term S*D*S
  %----------------------------------------------------------
  % Prepare full 3x3 D matrix
  if Sys.fullD
    D = Sys.D(3*(iSpin-1)+(1:3),:);
  else
    D = diag(Sys.D(iSpin,:));
  end
  if any(D(:))
    % Apply rotation if DFrame is given
    if any(Sys.DFrame(iSpin,:))
      R_M2D = erot(Sys.DFrame(iSpin,:)); % mol frame -> D frame
      R_D2M = R_M2D.';                   % D frame -> mol frame
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
  %   B1Frame, B2Frame, B3Frame, ...
  %   BFrame
  
  % If D is used, skip corresponding Stevens operator terms
  D_present = any(Sys.D(iSpin,:));
  
  % Run over all ranks k
  for k = 1:numel(Sys.B)
    Bk = Sys.B{k};
    if isempty(Bk), continue; end
    if D_present && k==2, continue; end
    if all(Bk==0), continue; end
    
    % Apply transformation if non-zero tilt angles are given
    % (Sys.BFrame is processed from Sys.B?Frame by validatespinsys)
    tiltang = Sys.BFrame{k};
    if any(tiltang)
      % Calculate transformation matrix for ISTOs
      Dk = wignerd(k,tiltang(1),tiltang(2),tiltang(3));
      Ck = isto2stev(k);
      % Transform to Stevens operator basis
      DB = Ck*Dk.'/Ck;
      DB = removenumericalnoise(DB);
      % Transform Bk
      Bk_M = Bk*DB; % eigenframe of Bk -> molecular frame
    else
      Bk_M = Bk;
    end
    
    % Build Hamiltonian
    q = k:-1:-k;
    for iq = find(Bk_M(iSpin,:)~=0)
      H = H + Bk_M(iSpin,iq)*stev(Spins,[k,q(iq),iSpin],'sparse');
    end
    
  end % for all tensor ranks

end % for all electron spins specified

H = (H+H')/2; % Hermitianize
if ~sparseResult
  H = full(H);
end

return
%===============================================================================

function A_ = removenumericalnoise(A)
reA = real(A);
imA = imag(A);
thr = 1e-14;
reA(abs(reA)<thr&reA~=0) = 0;
imA(abs(imA)<thr&imA~=0) = 0;
A_ = complex(reA,imA);
