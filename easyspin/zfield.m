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

  if isfield(Sys,'aF') && any(Sys.aF(:))
    % work only for first electron spin
    if iSpin~=1
      continue
    end
    % not available if D frame is tilted (would necessitate rotation
    % of the a and F terms which is not implemented).
    if isfield(Sys,'DFrame')
      if ~isempty(Sys.DFrame) && any(Sys.DFrame(:))
        error('It''s not possible to use Sys.aF with a tilted D frame (Sys.DFrame).');
      end
    end
    S = Spins(1);
    n = S*(S+1);
    Sz = sop(Spins,[1,3],'sparse');
    O40 = (35*Sz^4-30*n*Sz^2+25*Sz^2-(6*n-3*n^2)*speye(length(Sz)));
    F = Sys.aF(2);
    if (F~=0)
      H = H + (F/180)*O40;
    end
    a = Sys.aF(1);
    if a~=0
      Sp = sop(Spins,[1,4],'sparse');
      Sm = sop(Spins,[1,5],'sparse');
      if ~isfield(Sys,'aFFrame'), Sys.aFFrame = 4; end
      if (Sys.aFFrame==3)
        % along threefold axis (see Abragam/Bleaney p.142, p.437)
        O43 = (Sz*(Sp^3+Sm^3)+(Sp^3+Sm^3)*Sz)/2;
        H = H - 2/3*(a/120)*(O40 + 10*sqrt(2)*O43);
      elseif (Sys.aFFrame==4)
        % along fourfold (tetragonal) axis (used by some)
        O44 = (Sp^4+Sm^4)/2;
        H = H + (a/120)*(O40 + 5*O44);
      else
        error('Unknown Sys.aFFrame value. Use 3 for trigonal and 4 for tetragonal (collinear with D).');
      end
    end
  end

  % High-order terms in extended Stevens operator format
  %-----------------------------------------------------------------------------
  % Spin system fields:
  %   B1, B2, B3, ... (k = 1...12) -> processed to Sys.B
  %   B1Frame, B2Frame, B3Frame, ...
  %   BFrame
  
  % If D and aF are used, skip corresponding Stevens operator terms
  D_present = any(Sys.D(iSpin,:));
  aF_present = any(Sys.aF(:));
  
  % Run over all ranks k
  for k = 1:numel(Sys.B)
    Bk = Sys.B{k};
    if isempty(Bk), continue; end
    if D_present && k==2, continue; end
    if aF_present && k==4, continue; end
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
