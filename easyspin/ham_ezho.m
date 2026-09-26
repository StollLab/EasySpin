% ham_ezho  Multiple-order electron Zeeman interaction Hamiltonian
%
%   H = ham_ezho(SpinSystem,B)
%   H = ham_ezho(SpinSystem,B,Spins)
%   H = ham_ezho(SpinSystem,B,Spins,'sparse')
%   H = ham_ezho(SpinSystem,B,Spins,'sparse',lB)
%
%   [G0,G1,...] = ham_ezho(SpinSystem)
%   [G0,G1,...] = ham_ezho(SpinSystem,Spins)
%   [G0,G1,...] = ham_ezho(SpinSystem,Spins,'sparse')
%   [G_lB1,G_lB2,...] = ham_ezho(SpinSystem,Spins,'sparse',lB)
%   cG = ham_ezho(SpinSystem)
%   cG = ham_ezho(SpinSystem,[],Spins,'sparse',lB)
%
%   Returns the multiple-order electron Zeeman interaction Hamiltonian for the
%   electron spins Spins of the spin system SpinSystem. If a magnetic field B
%   is given, the Hamiltonian matrix H is returned. Otherwise, the tensors G0,
%   G1, G2, ... that contain the terms of order 0, 1, 2, ... in B are returned.
%
%   Input:
%   - SpinSystem: Spin system structure. The Hamiltonian is specified by the
%     fields Ham<lB><lS><l>, e.g. Ham112. Here, lB (0 to 8) is the order in the
%     magnetic field, lS (1 to 8) is the rank of the spin operators, and l is
%     the rank of the coupling, with |lB-lS| <= l <= lB+lS. Each field
%     contains the 2*l+1 real parameters (components m = l,...,-l) for each
%     electron spin (one row per electron spin) in units of MHz/mT^lB. A
%     column vector with one value per electron spin is taken as the m = 0
%     component. Fields with odd lB+lS (time-reversal symmetry) or odd l are
%     ignored, and lS = 0 is not supported.
%     lB = 0 corresponds to the zero-field terms of ham_zf, lB = 1 to the
%     Zeeman terms of ham_ez.
%   - B: Magnetic field vector in the molecular frame, in mT. If B is given,
%     the Hamiltonian matrix is returned. If B is omitted or empty, tensors
%     are returned.
%   - Spins: Vector of electron spin indices. If omitted or empty, all
%     electron spins are included.
%   - 'sparse': If given, results are returned in sparse format.
%   - lB: Order(s) in B. If given, only terms of these orders are returned.
%     Default: 0:8 for H, 0:highest for tensors. For tensors, lB can be 0 to 3.
%
%   Output:
%   - H: Hamiltonian matrix (MHz), for the given B.
%   - G0,G1,G2,G3: Tensors in units of MHz/mT^n, for n = 0, 1, 2, 3. They
%     are the coefficients of the terms of order n in B, so that
%       H = G0 + sum_i G1{i}*B(i) + sum_ij G2{i,j}*B(i)*B(j)
%             + sum_ijk G3{i,j,k}*B(i)*B(j)*B(k)
%     where i, j, k = 1, 2, 3 are the x, y, z components of B. G0 is a matrix,
%     G1 is a 1x3 cell array, G2 is a 3x3 cell array, and G3 is a 3x3x3 cell
%     array of matrices. G2 and G3 are symmetric with respect to their
%     indices. The Gn are the n-th derivatives of H divided by n!, e.g.
%     G3{1,2,3} = 1/6*d^3H/(dBx dBy dBz).
%     Tensors are available only up to order 3. Fields of higher order are
%     ignored. Without the input lB, nargout outputs give the orders
%     0,...,nargout-1. With lB, the number of outputs must be numel(lB).
%   - cG: Cell array {G0,G1,...} of tensors, with an entry for each lB
%     (default: 0:highest order in the spin system).
%
%   Uses the Hamiltonian as given in
%   McGavin, Tennant, Weil, J. Magn. Reson. 87, 92-109 (1990)
%
%   It is complete in the sense that it contains all usual electron spin
%   terms for a single spin, but it uses a non-standard convention. Nuclear
%   and hyperfine terms are not implemented.

function varargout = ham_ezho(SpinSystem, varargin)

if nargin==0, help(mfilename); return; end

if nargin>5, error('Wrong number of input arguments!'); end

% Tensors are returned if several outputs are requested or if no field is given
tensorOutput = nargout>1 || nargin<2 || isempty(varargin{1});

if tensorOutput
  % Zero-field Hamiltonian and field-dependent tensors up to 3rd order in B

  % Single output has an empty field input, multiple outputs do not
  if nargout>1
    args = varargin;  % Spins, opt, lB
  else
    args = varargin(2:end);  % B (empty), Spins, opt, lB
  end
  if numel(args)>3, error('Wrong number of input arguments!'); end
  Spins = [];
  opt = '';
  lb = [];
  if numel(args)>=1, Spins = args{1}; end
  if numel(args)>=2, opt = args{2}; end
  if numel(args)>=3, lb = args{3}; end
  if ~ischar(opt)
    error('Last input must be a string, ''sparse''.');
  end

  % Get highest order in B
  fields = fieldnames(SpinSystem);
  highest = 0;
  for n = 0:3
    if any(strncmp(fields,sprintf('Ham%i',n),4))
      highest = n;
    end
  end

  % Determine orders in B to return
  if isempty(lb)
    if nargout>1
      lb = 0:nargout-1;
    else
      lb = 0:highest;
    end
  elseif nargout>1 && nargout~=numel(lb)
    error('Number of outputs must match the number of requested orders in B!');
  end
  if ~all(ismember(lb,0:highest))
    error('Requested order in B is higher than that of the spin system (maximum is 3)!');
  end

  if any(lb==0)
    G0 = ham_ezho(SpinSystem,[0,0,0],Spins,opt,0);
  end
  xyz = 1:3;
  if any(lb==1)
    for n = 3:-1:1
      Field = zeros(1,3); Field(n) = 1;
      G1{n} = ham_ezho(SpinSystem,Field,Spins,opt,1);
    end
  end

  % For orders 2 and 3, the terms are homogeneous polynomials in B. Their
  % coefficients are extracted from the Hamiltonian at suitable test fields.
  if any(lb==2)
    % Diagonal elements from fields along x, y, z
    for n = 3:-1:1
      Field = zeros(1,3); Field(n) = 1;
      G2{n,n} = ham_ezho(SpinSystem,Field,Spins,opt,2);
    end
    % Off-diagonal elements from fields with one component zero
    for n = 1:3
      Field = ones(1,3); Field(n) = 0;
      mt = find(xyz~=n);
      G2{mt(1),mt(2)} = 1/2*(ham_ezho(SpinSystem,Field,Spins,opt,2)...
        -G2{mt(1),mt(1)}-G2{mt(2),mt(2)});
      G2{mt(2),mt(1)} = G2{mt(1),mt(2)};
    end
  end
  if any(lb==3)
    % Elements (n,n,n) from fields along x, y, z
    for n = 3:-1:1
      Field = zeros(1,3); Field(n) = 1;
      G3{n,n,n} = ham_ezho(SpinSystem,Field,Spins,opt,3);
    end
    % Elements (p,q,q) and permutations from fields with one component zero
    % and the other two equal to +1 or (+1,-1)
    for n = 1:3
      Field = ones(1,3); Field(n) = 0;
      lp = ham_ezho(SpinSystem,Field,Spins,opt,3);
      mt = perms(find(xyz~=n));
      for m = size(mt,1):-1:1
        Field(mt(m,2)) = -1;
        lm = ham_ezho(SpinSystem,Field,Spins,opt,3);
        me = 1/6*(lp+lm-2*G3{mt(m,1),mt(m,1),mt(m,1)});
        ind = perms([mt(m,1),mt(m,2),mt(m,2)]);
        for k = 1:size(ind,1)
          G3{ind(k,1),ind(k,2),ind(k,3)} = me;
        end
        Field(mt(m,2)) = 1;
      end
    end
    % Elements (x,y,z) and permutations from the field (1,1,1), after
    % subtracting all elements determined above
    Field = ones(1,3);
    ind = perms([1,2,3]);
    dif = 0*G3{1,1,1};
    for a = 1:3
      for b = 1:3
        for c = 1:3
          if ~isempty(G3{a,b,c})
            dif = dif + G3{a,b,c};
          end
        end
      end
    end
    me = 1/size(ind,1)*(ham_ezho(SpinSystem,Field,Spins,opt,3)-dif);
    for k = 1:size(ind,1)
      G3{ind(k,1),ind(k,2),ind(k,3)} = me;
    end
  end

  G = cell(1,numel(lb));
  for n = 1:numel(lb)
    switch lb(n)
      case 0, G{n} = G0;
      case 1, G{n} = G1;
      case 2, G{n} = G2;
      case 3, G{n} = G3;
    end
  end
  if nargout>1
    varargout = G;
  else
    varargout = {G};
  end

else  % full Hamiltonian is output

  Field = varargin{1};
  if nargin<3, Spins = []; else, Spins = varargin{2}; end
  if nargin<4, opt = ''; else, opt = varargin{3}; end
  if nargin<5, lBlist = 0:8; else, lBlist = varargin{4}; end

  if ~ischar(opt)
    error('Last input must be a string, ''sparse''.');
  end
  sparseResult = strcmp(opt,'sparse');

  if numel(Field)~=3
    error('Magnetic field vector (2nd input) must be a 3-element array.');
  end
  if ~all(ismember(lBlist,0:8))
    error('Orders in B (5th input) must be integers between 0 and 8.');
  end

  % Validate spin system
  [Sys,err] = validatespinsys(SpinSystem);
  error(err);

  % Vector of spin quantum numbers
  SpinVec = Sys.Spins;

  % Get number of electrons and states
  nElectrons = Sys.nElectrons;
  nStates = Sys.nStates;

  % No 'Spins' specified -> use all
  if isempty(Spins), Spins = 1:nElectrons; end

  % Validate third argument (Spins)
  if any(Spins<1) || any(Spins>nElectrons)
    error('Electron spin indices (3rd input argument) invalid!');
  end

  % Convert B from cartesian to spherical coordinates
  rB = norm(Field);
  [phiB,theta] = vec2ang(Field);
  ctheta = cos(theta);

  % Table of normalization factors alpha_l = (2*l-1)!!/l!, for l = 0,...,8
  % With these, the spherical tensor components of B are T_l,l = (-1)^l*2^(-l/2)*B_+^l,
  % using the same convention as for the spin operators T_l,m(S) below.
  alphapm1 = [1, 1, 3/2, 5/2, 35/8, 63/8, 231/16, 429/16, 6435/128];
  alphapm1 = 1./sqrt(alphapm1);

  % Table of conversion factors between spherical tensor operators T_l,m
  % and Stevens operators O_l,m (Alm(l,m+1) for m = 0,...,l), as used in stev():
  %   T_l,0 = O_l,0/Alm(l,1)
  %   T_l,m = (-1)^m/(sqrt(2)*Alm(l,m+1))*(O_l,m + i*O_l,-m),  m > 0
  %   T_l,-m = 1/(sqrt(2)*Alm(l,m+1))*(O_l,m - i*O_l,-m),  m > 0
  Alm(8,:) = [24*sqrt(1430),2*sqrt(1430),4*sqrt(143/7),2*sqrt(78/7), ...
    4*sqrt(130/7), 2*sqrt(10/7), 4*sqrt(15), 2*sqrt(2), 8*sqrt(2)];
  Alm(7,1:8) = [4*sqrt(429), 8*sqrt(429/7),4*sqrt(286/7),8*sqrt(143/7),...
    4*sqrt(13/7), 8*sqrt(13/7),4*sqrt(2/7),8];
  Alm(6,1:7) = [4*sqrt(231),2*sqrt(11),4*sqrt(22/5),2*sqrt(22/5),...
    4*sqrt(11/3), 2*sqrt(2/3), 4*sqrt(2)];
  Alm(5,1:6) = [6*sqrt(14), 2*sqrt(42/5),sqrt(6/5),12/sqrt(5), 2*sqrt(2/5), 4];
  Alm(4,1:5) = [2*sqrt(70), sqrt(7), sqrt(14), 1, 2*sqrt(2)];
  Alm(3,1:4) = [sqrt(10), 2*sqrt(5/3), sqrt(2/3),2];
  Alm(2,1:3) = [sqrt(6), 1/sqrt(2), sqrt(2)];
  Alm(1,1:2) = [1,1];

  fields = fieldnames(Sys);
  H = sparse(nStates,nStates);  % Hamiltonian
  for idx = 1:numel(Spins)
    iSpin = Spins(idx);

    % Run over all orders lB in B (Ham0.., Ham1.., Ham2.., ...)
    for k = 1:length(lBlist)
      lB = lBlist(k);
      strlB = sprintf('Ham%i',lB);
      paramtext = fields(strncmp(fields,strlB,4));
      if isempty(paramtext), continue; end

      % Spherical tensor components T_lB,mB (mB = lB,...,-lB) of B
      mBs = lB:-1:-lB;
      amB = abs(mBs);
      % Associated Legendre functions for m = 0,...,lB (they include the
      % Condon-Shortley phase). The extra phase factor is (-1)^|m| for m<0
      % and 1 for m>=0, since T_l,-m = (-1)^m*conj(T_l,m).
      LlBmB = legendre(lB,ctheta);
      pre = sqrt(factorial(lB-amB)./factorial(lB+amB)).*(-1).^(mBs.*(mBs<0));
      TlBmB = alphapm1(lB+1)*pre.*LlBmB(amB+1).'.*exp(1i*mBs*phiB);

      % Get ranks lS of spin operators for which parameters are given
      lStemp = zeros(1,numel(paramtext));
      for n = 1:numel(paramtext)
        lStemp(n) = str2double(paramtext{n}(5));
      end
      lStemp = unique(lStemp);

      % Run over all allowed ranks lS in S: lB+lS has to be even, due to
      % the time-reversal symmetry of the Hamiltonian
      for n = find(~mod(lB+lStemp,2))
        lS = lStemp(n);
        if lS==0
          error('Sys.Ham%i0* is not supported. lS must be between 1 and 8.',lB);
        end
        strlBlS = sprintf('Ham%i%i',lB,lS);
        Glb = rB^lB/sqrt(2);
        minl = abs(lB-lS);
        maxl = lB+lS;
        paramtext = fields(strncmp(fields,strlBlS,5));
        l_ = zeros(1,numel(paramtext));
        for n_ = 1:numel(paramtext)
          l_(n_) = str2double(paramtext{n_}(6:end));
        end
        l_ = l_(minl<=l_ & l_<=maxl & ~mod(l_,2));  % only even l

        for l = l_
          strlBlSl = sprintf('Ham%i%i%i',lB,lS,l);
          ZBlBlSlm = Sys.(strlBlSl)(iSpin,:);
          if ~any(ZBlBlSlm), continue; end

          % Construct the complex coefficients a^lB,lS,l_m of the coupling
          % from the real parameters ZB^lB,lS,l_m. Index iq = 1,...,2*l+1
          % corresponds to m = l,...,-l for both.
          alBlSlm = zeros(1,2*l+1);
          for m = l:-1:1
            lm = l+1+m;
            lp = l+1-m;
            alBlSlm(lp) = Glb*(-1)^m * (ZBlBlSlm(lp)-1i*ZBlBlSlm(lm));
            alBlSlm(lm) = Glb* (ZBlBlSlm(lp)+1i*ZBlBlSlm(lm));
          end
          alBlSlm(l+1) = Glb * sqrt(2) * ZBlBlSlm(l+1);

          for iq = find(alBlSlm~=0)
            m = l+1-iq;
            prefac = (-1)^m*sqrt(2*l+1);
            for mB = lB:-1:-lB
              if abs(TlBmB(lB-mB+1))<1e-10, continue; end
              mS = m-mB;  % Wigner 3j symbol is nonzero only for mB+mS = m
              if -lS<=mS && mS<=lS
                % Spin operator T_lS,mS from Stevens operators
                if mS == 0
                  TlSmS = stev(SpinVec,[lS,mS,iSpin],'sparse')/Alm(lS,abs(mS)+1);
                elseif mS > 0
                  TlSmS = (-1)^mS/(Alm(lS,mS+1)*sqrt(2))*...
                    (stev(SpinVec,[lS,mS,iSpin],'sparse')+1i*stev(SpinVec,[lS,-mS,iSpin],'sparse'));
                else
                  TlSmS = 1/(Alm(lS,abs(mS)+1)*sqrt(2))*...
                    (stev(SpinVec,[lS,-mS,iSpin],'sparse')-1i*stev(SpinVec,[lS,mS,iSpin],'sparse'));
                end
                H = H + alBlSlm(iq)*prefac*wigner3j(lB,lS,l,mB,mS,-m)*...
                  TlBmB(lB-mB+1)*TlSmS;
              end
            end % mB loop
          end % loop over iq
        end % loop over allowed l values
      end % loop over lS values
    end % loop over lB
  end % loop over spins

  if ~sparseResult
    H = full(H);
  end
  varargout = {H};
end

end
