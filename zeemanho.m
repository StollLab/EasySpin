
% zeemanho  higher order Zeeman interaction Hamiltonian 
%
%
%   H = zeemanho(SpinSystem, B)
%   H = zeemanho(SpinSystem, B, Spins)
%   H = zeemanho(SpinSystem, B, Spins, 'sparse')
%   H = zeemanho(SpinSystem, B, Spins, 'sparse', lB)
%   H = zeemanho(SpinSystem, B, [], '', lB)
%

%   [G0,G1] = zeemanho(SpinSystem)
%   [G0,G1,G2] = zeemanho(SpinSystem)
%   [G0,G1,G2,G3] = zeemanho(SpinSystem)
%   GlB = zeemanho(SpinSystem, Spins, 'sparse',lB)
%   [G0,...] = zeemanho(SpinSystem, Spins, 'sparse')
%   cG = zeemanho(SpinSystem)
%   cG = zeemanho(SpinSystem, [], Spins, 'sparse',lB)
%   (for cell output and spin selection, empty field for B required)
%
%   Returns the higher order Zeeman interaction Hamiltonian for
%   the spins 'Spins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla. If B is ommited zero field
%   is used (for single output).
%   - Spins: Vector of spin numbers. For one electron spin: 1
%     is the electron, >=2 are the nuclei. For two electron
%     spins: 1 and 2 electrons, >=3 nuclei, etc. If Spins is
%     omitted, all spins are included.
%   - 'sparse': If given, results returned in sparse format.
%   - lB: If given only terms of order lB in the magnetic field are
%   returned
%
%   Output:
%   - H: the Hamiltonian of the Zeeman interaction.
%   -[G0,G1,G2,G3]: components of the higher order Zeeman interaction Hamiltonian
%     for the selected spins, defined by the 0th, 1th, 2nd, and 3rd derivative of
%     of the parts of the Hamiltonian which contains the magnetic up to
%     this order.  Gn contain (d^n/dB^3) H, a tensor of rank n.
%     So  G3{1,2,3} contain d^3 /(dB_x dB_y dB_z) H. 
%   -cG cell containing the G
%
%   Uses the Hamiltonian as given in
%   MgGavin, Tennant and Weil, Jour. Mag. Res. 87,92-109 (1990)
%
%   it is complete in the sense that it contain all usual terms (in
%   principal also the nuclear and hyperfine for a single nuclei, but this
%   is not implemented) for a single spin, but in a non-popular convention.

function varargout = zeemanho(SpinSystem, varargin)


if (nargin==0), help(mfilename); return; end

if (nargin<1) || (nargin>5), error('Wrong number of input arguments!'); end

if (nargout>1 || nargin<2 ) 
% if no B value is given or several outputs are provided
  TensorOutput = 1;
elseif( nargin>2 && isempty(varargin{1}))
  % check if field is non-numeric
   TensorOutput = 1;
elseif( nargin==1 )
  TensorOutput = 1;
else
  TensorOutput = 0;
end

%zero-field Hamiltonian and field dependent tensors up to 3rd order in B
%are provided
if TensorOutput
  if nargout>4
    error('Incorrect number of outputs! Tensor output implemented only up to 3th order in B0');
  end
  if (nargout==1)
    if (nargin<3), Spins = []; else Spins = varargin{2}; end
    if (nargin<4), opt = ''; else opt = varargin{3}; end
    if (nargin<5)
      lb=[];
      selection = 0;
    else
      lb = varargin{4};
      selection = 1;
    end
    if (nargin>5), error('Wrong number of input arguments!'); end
  else
    if (nargin<2), Spins = []; else Spins = varargin{1}; end
    if (nargin<3), opt = ''; else opt = varargin{2}; end
    if (nargin<4)
      lb=[];
      selection = 0;
    else
      lb = varargin{3};
      selection = 1;
    end
    if (nargin>4), error('Wrong number of input arguments!'); end
  end
  
   %get highest order in B0
   fields = fieldnames(SpinSystem);
   for n = 0:3
    if any(strncmp(fields,['ZB',num2str(n,'%i')],3));
      highest = n;
    end
   end
   if selection
     if any(lb>highest)
       error('Requested order in B0 is higher than that of the SpinSystem!');
     end
   else
     lb = 0:highest;
   end
   if (nargout ~= length(lb) && nargout>2 )
     error('Incorrect number of outputs!');
   end

  if any(lb==0), GO = zeemanho(SpinSystem, [0,0,0], Spins, opt,0);end;
  xyz = 1:3;
  if any(lb==1)
    for n=3:-1:1
      Field = zeros(1,3);Field(n) = 1;
      G1{n} = zeemanho(SpinSystem, Field, Spins, opt,1);
    end
  end
  if any(lb==2)
    for n=3:-1:1
      Field = zeros(1,3);Field(n) = 1;
      G2{n,n} = zeemanho(SpinSystem, Field, Spins, opt,2);
    end
    for n=1:3
      Field = ones(1,3);Field(n) = 0;
      mt = find(xyz~=n);
      G2{mt(1),mt(2)} = 1/2*(zeemanho(SpinSystem, Field, Spins, opt,2)...
        -G2{mt(1),mt(1)}-G2{mt(2),mt(2)});
      G2{mt(2),mt(1)} = G2{mt(1),mt(2)};
    end    
  end
  if any(lb==3)
    for n=3:-1:1
      Field = zeros(1,3);Field(n) = 1;
      G3{n,n,n} = zeemanho(SpinSystem, Field, Spins, opt,3);
    end
    for n=1:3
      Field = ones(1,3);Field(n) = 0;
      lp = zeemanho(SpinSystem, Field, Spins, opt,3);
      mt = perms(find(xyz~=n));
      for m = length(mt):-1:1 
        Field(mt(m,2)) = -1;
        lm = zeemanho(SpinSystem, Field, Spins, opt,3);
        me = 1/6 * (lp+lm-2*G3{mt(m,1),mt(m,1),mt(m,1)});
        ind = perms([mt(m,1),mt(m,2),mt(m,2)]);
        for k =1:length(ind)
          G3{ind(k,1),ind(k,2),ind(k,3)} = me;
        end
        Field(mt(m,2)) = 1;
      end
    end
    Field = ones(1,3);
    ind = perms([1,2,3]);
    len = length(ind);
    si = size(G3{1,1,1});
    dif = zeros(si);
    for a =1:3, for b =1:3, for c = 1:3
          if (size(G3{a,b,c}) == si)
            dif = dif + G3{a,b,c};
          end
    end,end,end    
    me = 1/len*(zeemanho(SpinSystem, Field, Spins, opt,3)-dif);
    for k =1:len
      G3{ind(k,1),ind(k,2),ind(k,3)} = me;
    end
  end

    if(nargout == 1)
      for n=length(lb):-1:1
        switch lb(n)
          case 0, varargout{1}{n} = GO;
          case 1, varargout{1}{n} = G1;
          case 2, varargout{1}{n} = G2;
          case 3, varargout{1}{n} = G3;
        end
      end
    else
      for n=length(lb):-1:1
        switch lb(n)
          case 0, varargout{n} = GO;
          case 1, varargout{n} = G1;
          case 2, varargout{n} = G2;
          case 3, varargout{n} = G3;
        end
      end
    end
else  %full Hamiltonian is provided
  if nargin<2
    Field = zeros(1,3);
  else
    Field = varargin{1};
  end
  if (nargin<3), Spins = []; else Spins = varargin{2}; end
  if (nargin<4), opt = ''; else opt = varargin{3}; end
  if (nargin<5), lBlist=0:8; else lBlist = varargin{4}; end
  
  
  if ~ischar(opt)
    error('Last input must be a string, ''sparse''.');
  end
  sparseResult = strcmp(opt,'sparse');
  
  
  
  % Validate spin system
  [Sys,err] = validatespinsys(SpinSystem);
  error(err);
  
  % Vector of spin quantum numbers
  SpinVec = Sys.Spins;
  
  % No 'Spins' specified -> use all
  if isempty(Spins), Spins = 1:numel(SpinVec); end
  
  % Validate second argument (Spins)
  if any(Spins<1) || any(Spins>length(SpinVec))
    error('Spin indices (2nd input argument) invalid!');
  end

  % Get number of electrons, nuclei and states
  nElectrons = Sys.nElectrons;
  nStates = Sys.nStates;
  error(err);
  
  %convert B form cartesian to spherical
  [phiB, t, rB]= cart2sph(Field(1),Field(2),Field(3));
  ctheta= cos(pi/2-t); %dull definition of cart2sph
  
  %Table of conversion factors for spherical harmonics
  alphapm1 = 1./sqrt([1, 1,3/2,5/2,35/8,63/8,231/16,429,16,6435/128]);
  
  %Table of conversion factors for Stevens operator
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
  
  %constant G in MHz/mT
  Gconst =gfree*bmagn/(planck*1e9);
  
  hZ = sparse(nStates,nStates);
  for idx = 1:numel(Spins)
    iSpin = Spins(idx);
    if (iSpin<=nElectrons),  % If it's an electron...
      % Run over all ranks lb in B (ZB^lb,ls = ZB1l, ZB2l, ZB3l, ZB4l, ...)
      for k = 1:length(lBlist)
        lB =lBlist(k);
        lStemp = 1:min((2*SpinVec(iSpin)),8); %Stevens operators up to 12th order supported, conversion factors up to 8th order
        
        % run over all allowed ranks ls in S
        for n=find(mod(lB+lStemp,2)-1) %lb+ls has to be even, time-inversion sym of Hamiltonian
          fieldname = sprintf('ZB%d%d',lB,lStemp(n));
          lS = lStemp(n);
          if ~isfield(Sys,fieldname), continue; end
          ZBlBlS = Sys.(fieldname);
          if isempty(ZBlBlS), continue; end
          
          mB = lB:-1:-lB;
          amB =abs(mB);
          LlBmB = legendre(lB,ctheta);
          pre = sqrt(factorial(lB-amB)./factorial(lB+amB)).*...
            ((-1).^(heaviside(-mB).*mB)); % 1 for q<=0 (-1)^q for q>0
          TlBmB = alphapm1(lB+1)*pre.*LlBmB(amB+1).'.*exp(1i*mB*phiB);
          
          Glb = (Gconst*rB)^lB/sqrt(2);
          minl = abs(lB-lS);
          maxl = lB+lS;
          
          if ~(isfield(ZBlBlS,'l') && isfield(ZBlBlS,'vals'))
            error('%s is not properly defined. Fields l and vals are required!',fieldname);
          end
          
          len = length(ZBlBlS.l);
          
          for indl = find(minl<=ZBlBlS.l & ZBlBlS.l <=maxl & ~mod(ZBlBlS.l,2))
            l = ZBlBlS.l(indl);
            if iscell(ZBlBlS.vals)
              ZBlBlSlm = ZBlBlS.vals{indl};
            elseif len>1
              error('vals in %s has to be cell array when more then one l is given', fieldname);
            else
              ZBlBlSlm = ZBlBlS.vals;
            end
            
            if length(ZBlBlSlm)~= (2*l+1)
              error('%s is not properly defined. Cell vals{%d} must have %d entries!'...
                ,fieldname,l, 2*l+1);
            end
            if ~any(ZBlBlSlm), continue; end
            %construc a^lB lS_l m from ZB^lS lM
            for m = l:-1:1
              lm = l+1+m;
              lp = l+1-m;
              alBlSlm(lp) = Glb*(-1)^m * (ZBlBlSlm(lp)-1i*ZBlBlSlm(lm));
              alBlSlm(lm) = Glb* (ZBlBlSlm(lp)+1i*ZBlBlSlm(lm));
            end
            alBlSlm(l+1) = Glb * sqrt(2) * ZBlBlSlm(l+1);
            
            for iq = find(alBlSlm~=0)
              %identify m
              m = l+1-iq;
              pre = (-1)^m*sqrt(2*l+1);
              for mB=lB:-1:-lB
                if abs(TlBmB(lB-mB+1))<1e-10, continue; end
                mS = m-mB;% Wigner3j is ~= 0 for m = mB+mS
                if -lS<=mS && mS<=lS
                  if mS == 0
                    TlSmS = stev(SpinVec, lS,mS, iSpin)/Alm(lS,abs(mS)+1);
                  elseif mS > 0
                    TlSmS = (-1)^mS/(Alm(lS,mS+1)*sqrt(2))*...
                      (stev(SpinVec, lS,mS, iSpin)+1i*stev(SpinVec, lS,-mS, iSpin));
                  else
                    TlSmS = 1/(Alm(lS,abs(mS)+1)*sqrt(2))*...
                      (stev(SpinVec, lS,-mS, iSpin)-1i*stev(SpinVec, lS,mS, iSpin));
                  end
                  hZ = hZ + alBlSlm(iq)*pre*wigner3j(lB,lS,l,mB,mS,-m)*...
                    TlBmB(lB-mB+1)*TlSmS;
                end
              end % mB loop
            end % loop over iq
            clear alBlSlm
          end % loop over allowed l values
        end % loop over ls values
      end % loop over lb
    end % electron condition
  end %loop over spin centers
  
  if ~sparseResult
    hZ = full(hZ);
  end
  varargout = {hZ};
end

