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
%                  'mu'       single-center magnetic moment along field
%                  'mumol'    molar magnetic moment (magnetization) along field
%                  'muBM'     single-center magnetic moment along field,
%                               as multiple of Bohr magnetons
%                  'mueff'    effective magnetic moment (unitless)
%                  'chi'      single-center magnetic susceptibility,
%                               component along field
%                  'chimol'   molar magnetic susceptibility,
%                               component along field
%                  'chimolT'  chimol times temperature
%                  '1/chimol' inverse of chimol
%                  The default is 'muBM chimol'.
%        Units     'SI' (for SI units, default) or 'CGS' (for CGS-emu units)
%        Method    calculation method, 'operator' (default) or 'energies'
%        nKnots    number of knots for powder average
%
%
%    Output (if Opt.Output is not given):
%      muzBM    magnetic moment along zL axis, as multiple of Bohr magnetons
%                 (value is the same in SI and CGS-emu)
%      chizz    molar susceptibility, zLzL component, in SI units (m^3 mol^-1)
%
%    zL is the direction of the applied static magnetic field.
%
%    The size of outputs is nB x nT, where nB is the number of
%    field values in Exp.Field and nT is the number of temperature values
%    in Exp.Temperature.
%
%   If no output argument is given, the computed data are plotted.

function varargout = curry(Sys,Exp,Opt)

if (nargin==0), help(mfilename); return; end

if (nargin<2), Exp = struct; end
if (nargin<3), Opt = struct; end
if (nargin>3)
  error('Up to three input arguments are possible. You gave %d.',nargin);
end

if ~isstruct(Opt)
  if ischar(Opt)
    Opt = struct('Output',Opt);
  else
    error('Opt (third input argument) must be a structure or a string.');
  end
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
if ~isfield(Opt,'nKnots')
  Opt.nKnots = 10;
end
logmsg(1,'  number of knots: %d',Opt.nKnots);
if ~isfield(Opt,'Symmetry')
  Opt.Symmetry = []; % needed for p_symandgrid
end

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

% Parse output quantity list in Opt.Output
calculateMu = (nargout==0);
calculateChi = (nargout==0) || (nargout>1);
OneColumn = false;
keywords = strread(Opt.Output,'%s');
rmv = false(size(keywords));
for k = 1:numel(keywords)
  switch keywords{k}
    case 'mu', calculateMu = true;
    case 'mumol', calculateMu = true;
    case 'muBM', calculateMu = true;
    case 'mueff', calculateChi = true;
    case 'chi', calculateChi = true;
    case 'chimol', calculateChi = true;
    case 'chimolT', calculateChi = true;
    case '1/chimol', calculateChi = true;
    case 'onecolumn', rmv(k) = true; OneColumn = true;
    otherwise
      error('%s keyword in Opt.Output is not known.',keywords{k});
  end
end
keywords(rmv) = [];
logmsg(1,'  number of outputs: %d',numel(keywords));

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
  case 'energies'
    useOperatorMethod = false;
    if calculateChi && any(B(:))
      error('Energy method can calculate susceptibility only in the absence of magnetic field!');
    end
  otherwise
    error('Opt.Method can be ''operator'' or ''energies''!');
end
logmsg(1,'  calculation method: %s',Opt.Method);

% Set up Hamiltonian and magnetic dipole moment
%-------------------------------------------------
% zero-field Hamiltonian F (MHz)
% magnetic dipole moment operators muOpxM,muOpyM,muOpzM (MHz/mT)
% all are in the molecular frame
[H0,GxM,GyM,GzM] = sham(Sys);

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
% This sets Orientations, nOrientations, nSites and AverageOverChi
[Orientations,nOrientations,nSites,AvgOverChi] = p_crystalorientations(Exp,Opt);
Exp.OriWeights = Exp.OriWeights/4/pi;

beta = 1./T/boltzm;

% Initialize output arrays
muz = zeros(nFields,nTemperatures);
if calculateChi
  chizz = zeros(nFields,nTemperatures);
end

% Calculation loop
%-------------------------------------------------
% All calculations are done in the molecular frame.
if useOperatorMethod
  for iOri = 1:nOrientations
    [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
    
    % Projection of magnetic moment operator onto lab-frame axes
    % (field direction is along z axis of lab frame, zLab)
    muOpzL = zLab_M(1)*muOpxM + zLab_M(2)*muOpyM + zLab_M(3)*muOpzM; % J/T
    
    for iB = 1:nFields
      
      [V,E] = eig(H0 - B(iB)*muOpzL);
      E = diag(E); % J
      populations = exp(-(E-E(1))*beta);
      if zeroTemp
        populations(1) = 1;
      end
      
      % Calculate expectation value of zL-component of magnetic moment
      % for each spin state.
      muz_expect = real(diag(V'*muOpzL*V));
      
      % Calculate population-weighted average.
      muz_avg = (muz_expect.'*populations)./sum(populations,1);
      
      % Accumulate for powder average
      muz(iB,:) = muz(iB,:) + Exp.OriWeights(iOri)*muz_avg;
      
      if calculateChi
        
        % Determine step size for numerical derivative
        dB = eps^(1/3)*max(B(iB),1); % optimal step size for numerical derivative
        B1 = B(iB) + dB; dB = B1 - B(iB); % prevent round-off errors
        
        [V,E] = eig(H0 - (B(iB)+dB)*muOpzL);
        E = diag(E); % J
        populations = exp(-(E-E(1))*beta);
        if zeroTemp
          populations(1) = 1;
        end
        
        % Calculate population-weighted average.
        muz_expect = real(diag(V'*muOpzL*V));
        muz_avg2 = (muz_expect.'*populations)./sum(populations,1);
        
        % Calculate zz component of susceptibility tensor as numerical
        % derivative of muz. Accumulate for powder average.
        chizz_ = (muz_avg2-muz_avg)/dB;
        chizz(iB,:) = chizz(iB,:) + Exp.OriWeights(iOri)*chizz_;
        
      end
      
    end % loop over field values
  end % loop over orientations
  
else
  
  % Calculate M and chi using logarithm of partition function
  %----------------------------------------------------------------------
  
  % Pre-calculate log(Z) for zero field (used for the calculation of chi
  % at zero field)
  if calculateChi
    E0 = eig(H0);
    mEchi = E0(1);
    lE = length(E0);
    betachi = repmat(beta,lE,1);
    logZ0 = log(sum(exp(-repmat(E0-mEchi,1,nTemperatures).*betachi)));
  end
  
  for iOri = 1:nOrientations
    
    [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
    
    % Projection of magnetic moment operator onto lab z axis
    % (field direction along z axis of lab frame, zLab)
    muOpzL = zLab_M(1)*muOpxM + zLab_M(2)*muOpyM + zLab_M(3)*muOpzM;
    
    if calculateMu
      for iB = nFields:-1:1
        dB = eps^(1/3)*max(B(iB),1); % optimal step size for numerical derivative
        B1 = B(iB) + dB; dB = B1 - B(iB); % prevent round-off errors
        B2 = B(iB) - dB;
        E1 = eig(H0 - B1*muOpzL);
        E2 = eig(H0 - B2*muOpzL);
        Emin = min([E1;E2]);
        logZ1 = log(sum(exp(-(E1-Emin)*beta)));
        logZ2 = log(sum(exp(-(E2-Emin)*beta)));
        dlogZdB = (logZ1-logZ2)/(2*dB);
        mz(:,iB)= boltzm * Exp.Temperature(:) .* dlogZdB(:);
      end
      % Accumulation for powder average
      muz = muz + Exp.OriWeights(iOri)*mz.';
    end
    
    if calculateChi
      dB = eps^(1/3);
      E1 = eig(H0 + dB*muOpzL);
      logZ1 = log(sum(exp(-repmat(E1-mEchi,1,nTemperatures).*betachi)));
      chizz_ = boltzm *Exp.Temperature.' .* (logZ1-logZ0)*2/dB^2;
      % Accumulation for powder average
      chizz = chizz + Exp.OriWeights(iOri)*chizz_;
    end
    
  end  % loop over orientations
  
end

% Unit conversions
if calculateMu
  muz_SI = muz; % single-center magnetic moment, SI units
  if useCGSunits
    muz_CGS = muz_SI/1e-3; % single-center magnetic moment, CGS-emu units
  end
end
if calculateChi
  chizz_SI = chizz*mu0 + Sys.TIP/avogadro;   % single-center SI, add TIP
  if useCGSunits
    chizz_CGS = chizz_SI/(4*pi*1e-6);   % SI -> CGS unit conversion
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
    subplot(2,2,1)
    plot(x,muz_SI/bmagn); % SI and CGS value are numerically identical
    ylabel('\mu_z/\mu_B,  \mu_{mol,z}/(N_A\mu_B)')
    axis tight
    
    subplot(2,2,2)
    if useCGSunits
      plot(x,muz_CGS*avogadro);
      ylabel('\mu_{mol,z}  (erg G^{-1} mol^{-1})')
    else
      plot(x,muz_SI*avogadro);
      ylabel('\mu_{mol,z}  (J T^{-1} mol^{-1})')
    end
    axis tight
    
    subplot(2,2,3)
    if useCGSunits
      plot(x,chizz_CGS);
      ylabel('\chi_{mol,zz}  (cm^3 mol^{-1})');
    else
      plot(x,chizz_SI);
      ylabel('\chi_{mol,zz}  (m^3 mol^{-1})');
    end
    
    subplot(2,2,4)
    if useCGSunits
      plot(x,chizz_CGS.*T);
      ylabel('\chi_{mol,zz}T  (K cm^3 mol^{-1})');
    else
      plot(x,chizz_SI.*T);
      ylabel('\chi_{mol,zz}T  (K m^3 mol^{-1})');
    end
    
    for i = 1:4
      subplot(2,2,i);
      xlabel(xLabel);
      xlim([min(x) max(x)]);
    end
    
  else % 2D data
    cla
    Plot2D = (nFields>10);
    chizzT_SI = chizz_SI.*repmat(T(:).',nFields,1);
    if Plot2D
      subplot(2,2,1);
      surf(T,B,muz_SI/bmagn); % SI and CGS values are numerically identical
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
      plot(T,muz_SI/bmagn); % SI and CGS value are numerically identical
      axis tight
      ylabel('\mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
      
      subplot(2,2,2)
      if useCGSunits
        plot(T,muz_CGS*avogadro);
      else
        plot(T,muz_SI*avogadro);
      end
      axis tight
      ylabel('\mu_{mol,z} (J T^{-1}mol^{-1})')
      
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
      else
        plot(T,chizzT_SI);
      end
      ylabel('\chi_{mol,zz}T (K cm^3 mol^{-1})');
      
      for i = 1:4
        subplot(2,2,i);
        xlabel('temperature (K)');
        xlim([min(T) max(T)]);
      end
      
    end % if Plot2D
  end % DataDimensionality
end % if doPlot

%-------------------------------------------------------------------------------
% Assign output arguments
%-------------------------------------------------------------------------------

% Assign unit-ful values
if useCGSunits
  if calculateMu, muz = muz_CGS; end
  if calculateChi, chizz = chizz_CGS; end
  muB = bmagn/1e-3; % Bohr magneton in CGS units
else
  if calculateMu, muz = muz_SI; end
  if calculateChi, chizz = chizz_SI; end
  muB = bmagn; % Bohr magneton in SI units
end

% Units-independent output assignment
for n = numel(keywords):-1:1
  switch keywords{n}
    case 'mu'
      outval = muz;
    case 'mumol'
      outval = muz*avogadro;
    case 'muBM'
      outval = muz/muB; % dimensionless; numerically identical in SI and CGS-emu
    case 'mueff'
      if useCGSunits
        c = 1/0.1250494086;
        outval = sqrt(chizz*avogadro.*repmat(T(:).',nFields,1)*c);
      else
        c = 3*boltzm/avogadro/bmagn^2/mu0;
        outval = sqrt(chizz*avogadro.*repmat(T(:).',nFields,1)*c);
      end
    case 'chi'
      outval = chizz;
    case 'chimol'
      outval = chizz*avogadro;
    case 'chimolT'
      outval = chizz*avogadro.*repmat(T(:).',nFields,1);
    case '1/chimol'
      outval = 1./(chizz*avogadro);
    otherwise
      error('Keyword %s in Opt.Output is unknown.',keywords{n});
  end
  
  if OneColumn
    outdim = nFields*nTemperatures;
    for m = nTemperatures:-1:1
      templine((m-1)*nFields+1:m*nFields) = outval(:,m);
    end
    varargout{1}((n-1)*outdim+1:n*outdim) = templine;
  else
    varargout{n} = outval;
  end
  
end

logmsg(1,'=end=curry========%s=================\n',datestr(now));

clear global EasySpinLogLevel
