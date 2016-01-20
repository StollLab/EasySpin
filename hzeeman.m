% hzeeman  higher order Zeeman interaction Hamiltonian 
%
%   H = hzeeman(SpinSystem)
%   H = hzeeman(SpinSystem, B)
%   H = hzeeman(SpinSystem, B, Spins)
%   H = hzeeman(SpinSystem, B, Spins, 'sparse')
%
%   Returns the higher order Zeeman interaction Hamiltonian for
%   the spins 'Spins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla. If B is ommited zero field
%   is used.
%   - Spins: Vector of spin numbers. For one electron spin: 1
%     is the electron, >=2 are the nuclei. For two electron
%     spins: 1 and 2 electrons, >=3 nuclei, etc. If Spins is
%     omitted, all spins are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - H: the Hamiltonian of the Zeeman interaction.
%
%   Uses the Hamiltonian as given in
%   MgGavin, Tennant and Weil, Jour. Mag. Res. 87,92-109 (1990)
%
%   it is complete in the sense that it contain all usual terms (in
%   principal also the nuclear and hyperfine for a single nuclei, but this
%   is not implemented) for a single spin, but in a non-popular convention.

function hZ = hzeeman(SpinSystem, varargin)

if (nargin==0), help(mfilename); return; end

if (nargin<1) || (nargin>4), error('Wrong number of input arguments!'); end
if nargin<2
  Field = zeros(1,3);
else
  Field = varargin{1};
end
if (nargin<3), Spins = []; else Spins = varargin{2}; end
if (nargin<4), opt = ''; else opt = varargin{3}; end


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
    for lB = 0:8
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
end
