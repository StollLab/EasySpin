% zfield  Electronic zero field interaction Hamiltonian 
%
%   F = zfield(SpinSystem)
%   F = zfield(SpinSystem,Electrons)
%
%   Returns the electronic zero-field interaction (ZFI)
%   Hamiltonian [MHz] of the system SpinSystem.
%
%   If the vector Electrons is given, the ZFI of only the
%   specified electrons is returned (1 is the first, 2 the
%   second, etc). Otherwise, all electrons are included.

function H = zfield(SpinSystem,Electrons)

if (nargin==0), help(mfilename); return; end

[Sys,err] = validatespinsys(SpinSystem);
error(err);

if (nargin==1), Electrons = 1:Sys.nElectrons; end

if any(Electrons>Sys.nElectrons) || any(Electrons<1),
  error('Electron spin index/indices (2nd argument) out of range!');
end

H = sparse(Sys.nStates,Sys.nStates);
spvc = Sys.Spins;

for k = 1:length(Electrons)
  
  idx = Electrons(k);

  % S or I < 1 -> no internal interactions possible -> go to next spin
  if spvc(idx)<1, continue; end
  
  % Quadratic term S*D*S
  %----------------------------------------------------------
  % Prepare full D matrix
  if Sys.fullD
    D = Sys.D(3*(idx-1)+(1:3),:);
  else
    Ddiag = Sys.D(idx,:);
    if any(Ddiag)
      Rp = erot(Sys.Dpa(idx,:));
      D = Rp*diag(Ddiag)*Rp.';
    else
      D = 0;
    end
  end
  if any(D(:))
    % Construct spin operator matrices
    for c = 3:-1:1
      so{c} = sop(spvc,idx,c,'sparse');
    end
    % Construct SDS term
    for c1 = 1:3
      for c2 = 1:3
        H = H + D(c1,c2)*(so{c1}*so{c2});
      end
    end
  end

  % Fourth-order terms a and F
  %---------------------------------------------------------
  % Abragam/Bleaney p. 437
  % Bleaney/Trenam, Proc Roy Soc A, 223(1152), 1-14, (1954)
  % Doetschman/McCool, Chem Phys 8, 1-16 (1975)
  % Scullane, JMR 47, 383 (1982)
  % Jain/Lehmann, phys.stat.sol.(b) 159, 495 (1990)

  % work only for first electron spin
  if (idx~=1), continue; end

  if isfield(Sys,'aF')
    % not availabble if D frame is tilted
    % (this would necessitate rotation of the a and F terms
    % which is not implemented).
    if isfield(Sys,'Dpa')
      if ~isempty(Sys.Dpa) && any(Sys.Dpa)
        error('It''s not possible to use Sys.aF with a tilted D frame (Sys.Dpa).');
      end
    end
    S = spvc(1);
    n = S*(S+1);
    Sz = sop(spvc,1,3,'sparse');
    O40 = (35*Sz^4-30*n*Sz^2+25*Sz^2-(6*n-3*n^2)*speye(length(Sz)));
    a = Sys.aF(1);
    F = Sys.aF(2);
    if (a~=0)
      Sp = sop(spvc,1,4,'sparse');
      Sm = sop(spvc,1,5,'sparse');
      % along threefold axis (as in Abragam/Bleaney etc)
      if ~isfield(Sys,'aFrame'), Sys.aFrame = 3; end
      if (Sys.aFrame==3)
        O43 = (Sz*(Sp^3+Sm^3)+(Sp^3+Sm^3)*Sz)/2;
        H = H - a/180*(O40 + 10*sqrt(2)*O43);
      elseif (Sys.aFrame==4)
        % along fourfold (tetragonal) axis (used by some)
        O44 = (Sp^4+Sm^4)/2;
        H = H + a/120*(O40 + 5*O44);
      else
        error('Unknown Sys.aFrame value. Use 3 for trigonal and 4 for tetragonal (collinear with D).');
      end
    end
    if (F~=0)
      H = H + F/180*O40;
    end
  end

  % If D (and aF) are used, skip Stevens operator terms
  if any(Sys.D(idx,:)), continue; end


  % Extended Stevens operators
  %---------------------------------------------------------
  for k = 2:2:6
    if 2*spvc(idx)<k, break; end
    for q = 0:k
      fi = sprintf('B%d%d',k,abs(q));
      if isfield(Sys,fi)
        if any(Sys.D)
          warning('The extended Stevens operators are not used if Sys.D is given!!');
        else
          coeffs = Sys.(fi);
          if (numel(coeffs)<1) || (numel(coeffs)>1+(q>0))
            error('Wrong size of %s field in spin system structure!',fi);
          else
            H = H + coeffs(1)*stev(spvc,k,q);
            if (numel(coeffs)==2) && (q>0)
              H = H + coeffs(2)*stev(spvc,k,-q);
            end
          end
        end
      end
    end
  end

end % for all spins specified
H = full(H);
H = (H+H')/2; % Hermitianise

return
