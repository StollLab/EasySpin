% curry  Computation of magnetometry data (magnetic moment, susceptibility)
%
%   curry(Sys,Exp)
%   curry(Sys,Exp,Opt)
%   ... = curry(...)
%
%    Calculates the magnetic moment and the molar static magnetic
%    susceptibility of a powder sample for given values of
%    applied magnetic field and temperature.
%
%    Input:
%      Sys    spin system
%        TIP          temperature-independent susceptibility
%
%      Exp    experimental parameter
%        Field        list of field values (mT)
%        Temperature  list of temperatures (K)
%
%      Opt    calculation options
%        Output    string of keywords defining the outputs
%                  'mu'        single-center magnetic moment along field
%                  'mumol'     molar magnetic moment (magnetization) along field
%                  'muBM'      single-center magnetic moment along field,
%                                as multiple of Bohr magnetons
%                  'mueff'     effective magnetic moment (unitless)
%                  'chi'       single-center magnetic susceptibility,
%                                component along field
%                  'chimol'    molar magnetic susceptibility,
%                                component along field
%                  'chimolT'   chimol times temperature
%                  '1/chimol'  inverse of chimol
%                  The default is 'muBM chimol'.
%        Units     'SI' (for SI units, default) or 'CGS' (for CGS-emu units)
%        Method    calculation method, 'operator' (default) or 'energies'
%        GridSize  grid size for powder average
%        Spins     electron spin indices, for spin-selective calculation
%
%
%    Output depends on the settings in Opt.Output. If Opt.Output is not given,
%    two outputs are provided:
%      muzBM    magnetic moment along zL axis, as multiple of Bohr magnetons
%                 (value is the same in SI and CGS-emu)
%      chimol   molar susceptibility, zLzL component, in SI units (m^3 mol^-1)
%
%    The size of outputs is nB x nT, where nB is the number of field values in
%    Exp.Field and nT is the number of temperature values in Exp.Temperature.
%
%   If no output argument is given, the computed data are plotted.

function varargout = curry(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version
error(chkmlver);

if nargin<2, Exp = struct; end
if nargin<3, Opt = struct; end
if nargin>3
  error('Up to three input arguments are possible. You gave %d.',nargin);
end

if ~isstruct(Opt)
  error('Opt (third input argument) must be a structure or a string.');
end

if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

logmsg(1,['=begin=curry======' datestr(now) '=================']);

doPlot = (nargout==0);

% Spin system
%-------------------------------------------------------------------------------
logmsg(1,'Spin system');
if iscell(Sys)
  error('curry does not support calculations with multiple components.');
end
if ~isfield(Sys,'TIP')
  Sys.TIP = 0;
end

% Experimental parameters
%-------------------------------------------------------------------------------
logmsg(1,'Experimental parameters');

% Field
if ~isfield(Exp,'Field')
  Exp.Field = 0;
  disp('Exp.Field is missing, assuming zero field.');
end
B = Exp.Field/1e3; % magnetic field, T
nFields = numel(B);
logmsg(1,'  number of field values: %d',nFields);

% Temperature
if ~isfield(Exp,'Temperature')
  error('Exp.Temperature is missing.');
end

T = reshape(Exp.Temperature,1,[]); % temperature, K
nTemperatures = numel(T);
logmsg(1,'  number of temperature values: %d',nTemperatures);

if any(T<0)
  error('Negative temperatures are not possible.')
end
zeroTemp = any(T==0);

% Other fields
if ~isfield(Exp,'CrystalOrientation')
  Exp.CrystalOrientation = [];
end
if ~isfield(Exp,'CrystalSymmetry')
  Exp.CrystalSymmetry = '';
end
if ~isfield(Exp,'MolFrame')
  Exp.MolFrame = [];
end

doPowderSimulation = isempty(Exp.CrystalOrientation);

if doPowderSimulation
  logmsg(1,'Powder calculation');
else
  logmsg(1,'Crystal calculation');
end

% Options
%-------------------------------------------------------------------------------
logmsg(1,'Options');
if ~isfield(Opt,'GridSize')
  Opt.GridSize = 10;
end
logmsg(1,'  number of knots: %d',Opt.GridSize);
if ~isfield(Opt,'GridSymmetry')
  Opt.GridSymmetry = []; % needed for p_symandgrid
end
if ~isfield(Opt,'GridFrame')
  Opt.GridFrame = []; % needed for p_symandgrid
end
if ~isfield(Opt,'deltaB')
  Opt.deltaB = 1; % mT
end
dB = Opt.deltaB*1e-3; % mT -> T

% Parse output quantity list in Opt.Output
if ~isfield(Opt,'Output')
  switch nargout
    case 0
      Opt.Output = '';
    case 1
      Opt.Output = 'muBM';
    case 2
      Opt.Output = 'muBM chimol';
    otherwise
      error('Incorrect number of outputs (1 or 2 expected if Opt.Output is not given.');
  end
end
logmsg(1,'  output: %s',Opt.Output);

calculateMu = (nargout==0);
calculateChi = (nargout==0) || (nargout>1);
calculateMuVec = false;
keywords = strread(Opt.Output,'%s'); %#ok
for k = 1:numel(keywords)
  switch keywords{k}
    case 'mu', calculateMu = true;
    case 'mumol', calculateMu = true;
    case 'muBM', calculateMu = true;
    case 'mueff', calculateChi = true;
    case 'muvec', calculateMu = true; calculateMuVec = true;
    case 'chi', calculateChi = true;
    case 'chimol', calculateChi = true;
    case 'chimolT', calculateChi = true;
    case '1/chimol', calculateChi = true;
    otherwise
      error('''%s'' keyword in Opt.Output is not known.',keywords{k});
  end
end
logmsg(1,'  number of outputs: %d',numel(keywords));

if calculateMuVec && doPlot
  error('Cannot plot results when calculating full magentic moment vector.');
end

% Parse unit system in Opt.Units
if ~isfield(Opt,'Units')
  Opt.Units = 'SI';
end
switch upper(Opt.Units)
  case 'SI', useCGSunits = false;
  case 'CGS', useCGSunits = true;
  otherwise
    error('Unknown units ''%s'' in Opt.Units. Use either ''SI'' or ''CGS''.',...
      Opt.Units);
end

% Parse calculation method in Opt.Method
if ~isfield(Opt,'Method')
  Opt.Method = 'operator';
end
switch lower(Opt.Method)
  case 'operator'
    useOperatorMethod = true;
  case 'partitionfunction'
    useOperatorMethod = false;
  otherwise
    error('Opt.Method can be ''operator'' or ''partitionfunction''!');
end
logmsg(1,'  calculation method: %s',Opt.Method);

% Spin indices for spin-selective magnetic moment calculation
if ~isfield(Opt,'Spins')
  Opt.Spins = [];
elseif ~isempty(Opt.Spins)
  if any(Opt.Spins<1) || any(Opt.Spins>numel(Sys.S))
    error('Opt.Spins must contain electron spin indices between 1 and %d.',numel(Sys.S));
  end
  if ~useOperatorMethod
    error('To use Opt.Spins, you must use Opt.Method = ''operator''.');    
  end
end


% Set up Hamiltonian and magnetic dipole moment operators
%-------------------------------------------------------------------------------
% zero-field Hamiltonian F (MHz)
% magnetic dipole moment operators muOpxM,muOpyM,muOpzM (MHz/mT)
% all are in the molecular frame
[H0,GxM,GyM,GzM] = sham(Sys);
if ~isempty(Opt.Spins)
  [GxM,GyM,GzM] = zeeman(Sys,Opt.Spins);
end

% zero-field spin Hamiltonian
H0 = H0*1e6*planck; % MHz -> J

% magnetic-dipole operators, in molecular frame
c = 1e6*1e3*planck; % conversion factor, MHz/mT -> J/T
muOpxM = -GxM*c; % MHz/mT -> J/T
muOpyM = -GyM*c; % MHz/mT -> J/T
muOpzM = -GzM*c; % MHz/mT -> J/T

% Set up sample orientations
%-------------------------------------------------
Exp.PowderSimulation = doPowderSimulation; % for communication with p_*
[Exp,Opt] = p_symandgrid(Sys,Exp,Opt);

% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,~,~] = p_crystalorientations(Exp,Opt);
Exp.OriWeights = Exp.OriWeights/4/pi;

beta = 1./T/boltzm;

% Initialize output arrays
muz = zeros(nFields,nTemperatures);
if calculateMuVec
  mux = zeros(nFields,nTemperatures);
  muy = zeros(nFields,nTemperatures);
end
if calculateChi
  chizz = zeros(nFields,nTemperatures);
end

% Calculation of magnetic moment and susceptibility
%-------------------------------------------------------------------------------
% All calculations are done in the molecular frame.

% Define utility functions
calcmu = @(Vecs,muOp,pop) (real(diag(Vecs'*muOp*Vecs)).'*pop)./sum(pop,1);
getmuproj = @(nM) nM(1)*muOpxM + nM(2)*muOpyM + nM(3)*muOpzM;

% Orientation loop
for iOri = 1:nOrientations
  [xL_M,yL_M,zL_M] = erot(Orientations(iOri,:),'rows');
  
  % Projection of magnetic moment operator onto lab-frame axes
  % (field direction is along z axis of lab frame, zL)
  muOpzL = getmuproj(zL_M); % J/T
  if calculateMuVec
    muOpxL = getmuproj(xL_M); % J/T
    muOpyL = getmuproj(yL_M); % J/T
  end
  
  for iB = 1:nFields
    
    if useOperatorMethod
      
      % Calculate mu with magnetic-moment operators, chi as numerical derivative
      %-------------------------------------------------------------------------
      [V,E] = eig(H0 - B(iB)*muOpzL);
      E = diag(E); % J
      populations = exp(-(E-E(1))*beta);
      if zeroTemp
        populations(1) = 1;
      end
      
      % Calculate expectation value of zL-component of magnetic moment
      % for each spin state and do population average
      muz_ = calcmu(V,muOpzL,populations);
      if calculateMuVec
        mux_ = calcmu(V,muOpxL,populations);
        muy_ = calcmu(V,muOpyL,populations);
      end
      
      if calculateChi
        
        % Solve eigenproblem at slightly higher field
        B_ = B(iB) + dB;
        [V,E] = eig(H0 - B_*muOpzL);
        E = diag(E); % J
        populations = exp(-(E-E(1))*beta);
        if zeroTemp
          populations(1) = 1;
        end
        
        % Calculate population-weighted average.
        muz2_ = calcmu(V,muOpzL,populations);
        
        % Calculate zz component of susceptibility tensor as numerical
        % derivative of muz.
        chizz_ = (muz2_-muz_)/dB;
        
      end
      
    else
      
      % Calculate mu and chi using logarithm of partition function
      %-------------------------------------------------------------------------
      lnZ = @(E,Emin) log(sum(exp(-(E-Emin)*beta),1)); % log of partition function
      
      E1 = eig(H0 - (B(iB)-dB)*muOpzL);
      E3 = eig(H0 - (B(iB)+dB)*muOpzL);
      Emin = min([E1;E3]);
      lnZ1 = lnZ(E1,Emin);
      lnZ3 = lnZ(E3,Emin);
      if calculateChi
        E2 = eig(H0 - B(iB)*muOpzL);
        lnZ2 = lnZ(E2,Emin);
      end
      
      % Calculate mu via symmetric difference quotient of ln(Z).
      if calculateMu
        muz_ = (lnZ3-lnZ1)/(2*dB)./beta;
      end
      
      % Calculate chi via symmetric second difference quotient of ln(Z).
      if calculateChi
        chizz_ = (lnZ3-2*lnZ2+lnZ1)/(dB^2)./beta;
      end
      
    end
    
    % Accumulate
    if calculateMu
      muz(iB,:) = muz(iB,:) + Exp.OriWeights(iOri)*muz_;
    end
    if calculateMuVec
      mux(iB,:) = mux(iB,:) + Exp.OriWeights(iOri)*mux_;
      muy(iB,:) = muy(iB,:) + Exp.OriWeights(iOri)*muy_;
    end
    if calculateChi
      chizz(iB,:) = chizz(iB,:) + Exp.OriWeights(iOri)*chizz_;
    end
    
  end % loop over field values
  
end % loop over orientations


% Unit conversions
%-------------------------------------------------------------------------------
if calculateMu
  muz_SI = muz; % single-center magnetic moment, SI units
  if useCGSunits
    muz_CGS = muz_SI/1e-3; % single-center magnetic moment, CGS-emu units
  end
end
if calculateMuVec
  mux_SI = mux; % single-center magnetic moment, SI units
  muy_SI = muy; % single-center magnetic moment, SI units
  if useCGSunits
    mux_CGS = mux_SI/1e-3; % single-center magnetic moment, CGS-emu units
    muy_CGS = muy_SI/1e-3; % single-center magnetic moment, CGS-emu units
  end
end
if calculateChi
  chizz_SI = chizz*mu0 + Sys.TIP/avogadro;   % single-center SI, add TIP
  if useCGSunits
    chizz_CGS = chizz_SI/(4*pi*1e-6);   % SI -> CGS-emu unit conversion
  end
end


%-------------------------------------------------------------------------------
% Graphical plotting
%-------------------------------------------------------------------------------
if doPlot
  DataDimensionality = (nFields>1) + (nTemperatures>1);
  if DataDimensionality==0
    % do nothing
    clf
    
  elseif DataDimensionality==1
    if nFields==1
      x = T;
      xLabel = 'temperature (K)';
    else
      x = B;
      xLabel = 'magnetic field (T)';
    end
    
    clf
    
    % magnetic moment, z component
    subplot(2,2,1)
    plot(x,muz_SI/bmagn); % SI and CGS-emu value are numerically identical
    ylabel('\mu_z/\mu_B,  \mu_{mol,z}/(N_A\mu_B)')
    axis tight
    grid on
    
    % molar magnetic moment, z component
    subplot(2,2,2)
    if useCGSunits
      plot(x,muz_CGS*avogadro);
      ylabel('\mu_{mol,z}  (erg G^{-1} mol^{-1})')
    else
      plot(x,muz_SI*avogadro);
      ylabel('\mu_{mol,z}  (J T^{-1} mol^{-1})')
    end
    axis tight
    grid on
    
    % molar magnetic susceptibility, zz component
    subplot(2,2,3)
    if useCGSunits
      plot(x,chizz_CGS*avogadro);
      ylabel('\chi_{mol,zz}  (cm^3 mol^{-1})');
    else
      plot(x,chizz_SI*avogadro);
      ylabel('\chi_{mol,zz}  (m^3 mol^{-1})');
    end
    axis tight
    grid on
    
    % molar magnetic susceptibility, zz component, times temperature
    subplot(2,2,4)
    if useCGSunits
      plot(x,chizz_CGS*avogadro.*T);
      ylabel('\chi_{mol,zz}T  (K cm^3 mol^{-1})');
    else
      plot(x,chizz_SI*avogadro.*T);
      ylabel('\chi_{mol,zz}T  (K m^3 mol^{-1})');
    end
    axis tight
    grid on
    
    for i = 1:4
      subplot(2,2,i);
      xlabel(xLabel);
      xlim([min(x) max(x)]);
    end
    
  else % 2D data
    cla
    Plot2D = (nFields>10);
    chizzT_SI = chizz_SI.*repmat(T,nFields,1);
    if Plot2D
      subplot(2,2,1);
      surf(T,B,muz_SI/bmagn); % SI and CGS-emu values are numerically identical
      shading flat
      zlabel(' \mu_z/\mu_B,  \mu_{mol,z}/(N_A\mu_B)')
      
      subplot(2,2,2);
      if useCGSunits
        surf(T,B,muz_CGS*avogadro);
        zlabel(' \mu_{mol,z}  (J T^{-1} mol^{-1})')
      else
        surf(T,B,muz_SI*avogadro);
        zlabel(' \mu_{mol,z}  (erg G^{-1} mol^{-1})')
      end
      shading flat
      
      subplot(2,2,3);
      if useCGSunits
        surf(T,B,chizz_CGS);
        zlabel('\chi_{mol,zz}  (cm^3 mol^{-1})');
      else
        surf(T,B,chizz_SI);
        zlabel('\chi_{mol,zz}  (m^3 mol^{-1})');
      end
      shading flat
      
      subplot(2,2,4);
      if useCGSunits
        surf(T,B,chizzT_SI);
        zlabel('\chi_{mol,zz}T  (K cm^3 mol^{-1})');
      else
        surf(T,B,chizzT_SI);
        zlabel('\chi_{mol,zz}T  (K m^3 mol^{-1})');
      end
      shading flat
      
      for i = 1:4
        subplot(2,2,i);
        xlabel('T (K)');
        ylabel('B (T)')
        xlim([min(T) max(T)]);
        ylim([min(B) max(B)]);
      end
      
    else
      subplot(2,2,1)
      plot(T,muz_SI/bmagn); % SI and CGS-emu value are numerically identical
      axis tight
      ylabel('\mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
      
      subplot(2,2,2)
      if useCGSunits
        plot(T,muz_CGS*avogadro);
        ylabel('\mu_{mol,z} (erg G^{-1} mol^{-1})')
      else
        plot(T,muz_SI*avogadro);
        ylabel('\mu_{mol,z} (J T^{-1} mol^{-1})')
      end
      axis tight
      
      subplot(2,2,3)
      if useCGSunits
        plot(T,chizz_CGS);
        ylabel('\chi_{mol,zz} (cm^3 mol^{-1})');
      else
        plot(T,chizz_SI);
        ylabel('\chi_{mol,zz} (m^3 mol^{-1})');
      end
      
      subplot(2,2,4)
      if useCGSunits
        plot(T,chizzT_CGS);
        ylabel('\chi_{mol,zz}T (K cm^3 mol^{-1})');
      else
        plot(T,chizzT_SI);
        ylabel('\chi_{mol,zz}T (K m^3 mol^{-1})');
      end
      
      for i = 1:4
        subplot(2,2,i);
        xlabel('temperature (K)');
        xlim([min(T) max(T)]);
      end
      
    end % if Plot2D
  end % DataDimensionality
end % if doPlot

%-------------------------------------------------------------------------------
% Calculate and assign outputs
%-------------------------------------------------------------------------------

% Assign unit-ful values
if useCGSunits
  if calculateMu, muz = muz_CGS; end
  if calculateMuVec, mux = mux_CGS; muy = muy_CGS; end
  if calculateChi, chizz = chizz_CGS; end
  muB = bmagn/1e-3; % Bohr magneton in CGS-emu units (erg/G = abA cm^2)
  kB = boltzm/1e-7; % Boltzmann constant in CGS units (erg/K)
  c = 3*kB/muB^2; % conversion factor needed for 'mueff'
else
  if calculateMu, muz = muz_SI; end
  if calculateMuVec, mux = mux_SI; muy = muy_SI; end
  if calculateChi, chizz = chizz_SI; end
  muB = bmagn; % Bohr magneton in SI units (J/T = A m^2)
  kB = boltzm; % Boltmann constant in SI units (J/K)
  c = 3*kB/muB^2/mu0; % conversion factor needed for 'mueff'
end

% Calculate and assign outputs
for n = numel(keywords):-1:1
  
  switch keywords{n}
    case 'mu'
      outval = muz;
    case 'mumol'
      outval = muz*avogadro;
    case 'muBM'
      outval = muz/muB;
    case 'muvec'
      outval = [mux; muy; muz];
    case 'chi'
      outval = chizz;
    case 'chimol'
      outval = chizz*avogadro;
    case 'chimolT'
      outval = chizz*avogadro.*repmat(T,nFields,1);
    case '1/chimol'
      outval = 1./(chizz*avogadro);
    case 'mueff'
      outval = sqrt(chizz.*repmat(T,nFields,1)*c);
    case 'chitensor'
      error('chi tensor currently not implemented.');
    otherwise
      error('Keyword %s in Opt.Output is unknown.',keywords{n});
  end
  
  varargout{n} = outval;
  
end

logmsg(1,'=end=curry========%s=================\n',datestr(now));

clear global EasySpinLogLevel
