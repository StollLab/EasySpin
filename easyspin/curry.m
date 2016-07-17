% curry  Computation of magnetometry data
%         (magnetic moment, susceptibility)
%
%   curry(Sys,Exp)
%   curry(Sys,Exp,Opt)
%   muz = curry(...)
%   [muz,chizz] = curry(...)
%
%    Calculates the magnetic moment and the molar static magnetic
%    susceptibility of a powder sample for given values of
%    applied magnetic field and temperature.
%
%    Input:
%      Sys    spin system
%      Exp    experimental parameter
%        Field        list of field values (mT)
%        Temperature  list of temperatures (K)
%      Opt    calculation options
%        nKnots       number of knots for powder average
%
%    Output:
%      muz      magnetic moment along zL axis
%               in units of Bohr magnetons
%      chizz    molar susceptibility, zLzL component
%               in SI units (m^3 mol^-1)
%    
%    zL is the direction of the applied static magnetic field
%
%    The size of muz and chizz is nB x nT, where nB is the number of
%    field values in Exp.Field and nT is the number of temperature values
%    in Exp.Temperature.
%
%   If no output argument is given, the computed data are plotted.

function varargout = curry(Sys,Exp,Opt)

if (nargin==0), help(mfilename); return; end

if (nargin<2), Exp = struct; end
if (nargin<3), Opt = struct; end

if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

logmsg(1,['=begin=curry======' datestr(now) '=================']);

calculateChi = (nargout==0) || (nargout>1);
doPlot = (nargout==0);

% Spin system
%-------------------------------------------------
logmsg(1,'Spin system');
if iscell(Sys)
  error('curry does not support calculations with multiple components.');
end

% Experimental parameters
%-------------------------------------------------
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
T = Exp.Temperature(:).'; % temperature, K
nTemperatures = numel(T);
logmsg(1,'  number of temperature values: %d',nTemperatures);

if any(T<0)
  error('Negative temperatures are not possible.')
end
zeroTemp = T==0;

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

PowderSimulation = isempty(Exp.CrystalOrientation);

if PowderSimulation
  logmsg(1,'Powder calculation');
else
  logmsg(1,'Crystal calculation');
end

% Options
%-------------------------------------------------
logmsg(1,'Options');
if ~isfield(Opt,'nKnots')
  Opt.nKnots = 10;
end
logmsg(1,'  number of knots: %d',Opt.nKnots);
if ~isfield(Opt,'Symmetry')
  Opt.Symmetry = []; % needed for p_symandgrid
end

% Set up Hamiltonian and magnetic dipole moment
%-------------------------------------------------
% zero-field Hamiltonian F (MHz)
% magnetic dipole moment operators muOpxM,muOpyM,muOpzM (MHz/mT)
% all are in the molecular frame
[H0,GxM,GyM,GzM] = sham(Sys);

% zero-field spin Hamiltonian
H0 = H0*1e6*planck; % MHz -> J

% magnetic dipole operators, in molecular frame
muOpxM = -GxM*1e6*1e3*planck; % MHz/mT -> J/T
muOpyM = -GyM*1e6*1e3*planck; % MHz/mT -> J/T
muOpzM = -GzM*1e6*1e3*planck; % MHz/mT -> J/T

% Set up sample orientations
%-------------------------------------------------
Exp.PowderSimulation = PowderSimulation; % for communication with p_*
[Exp,Opt] = p_symandgrid(Sys,Exp,Opt);
% Process crystal orientations, crystal symmetry, and frame transforms
% This sets Orientations, nOrientations, nSites and AverageOverChi
[Orientations,nOrientations,nSites,AverageOverChi] = p_crystalorientations(Exp,Opt);
Exp.OriWeights = Exp.OriWeights/4/pi;

beta = 1./T/boltzm; 

% Initialize output arrays
muz = zeros(nFields,nTemperatures);
if calculateChi
  chizz = zeros(nFields,nTemperatures);
end

% Calculation loop
%-------------------------------------------------
% all calculations are done in the molecular frame

for iOri = 1:nOrientations
  [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
  % projection of magnetic moment operator onto lab axes
  % (field direction along z axis of lab frame, zLab)
  muOpzL = zLab_M(1)*muOpxM + zLab_M(2)*muOpyM + zLab_M(3)*muOpzM; % J/T
  
  for iB = 1:numel(Exp.Field)
    
    [V,E] = eig(H0 - B(iB)*muOpzL);
    E = diag(E) - E(1); % J
    populations = exp(-E*beta);
    if zeroTemp, populations(1) = 1; end
    
    % expectation values for each spin state
    muz_expect = real(diag(V'*muOpzL*V));
    
    % population-weighted average
    muz_avg = (muz_expect.'*populations)./sum(populations,1);
    
    % accumulation for powder average
    muz(iB,:) = muz(iB,:) + Exp.OriWeights(iOri)*muz_avg;
    
    if calculateChi
      dB = eps^(1/3)*max(B(iB),1); % optimal step size for numerical derivative
      B1 = B(iB) + dB; dB = B1 - B(iB); % prevent round-off errors
      
      [V,E] = eig(H0 - (B(iB)+dB)*muOpzL);
      E = diag(E) - E(1); % J
      populations = exp(-E*beta);
      if zeroTemp, populations(1) = 1; end
      
      % population-weighted average
      muz_expect = real(diag(V'*muOpzL*V));
      muz_avg2 = (muz_expect.'*populations)./sum(populations,1);
      
      chizz_ = (muz_avg2-muz_avg)/dB;
      chizz(iB,:) = chizz(iB,:) + Exp.OriWeights(iOri)*chizz_;
    end
    
  end
end

% Unit conversions
if  calculateChi
  chizz_SI = chizz*mu0*avogadro;   % single molecule SI -> molar SI
end

%chizz_cgs = chizz_SI/(4*pi*1e-6);   % SI -> CGS-emu unit conversion

% Graphical plotting
%-----------------------------------------------------
if doPlot
  DataDimensionality = (nFields~=1) + (nTemperatures~=1);
  if DataDimensionality==0
    % do nothing
    clf
    
  elseif DataDimensionality==1
    if (nFields==1)
      x = T;
      xLab = 'temperature (K)';
    else
      x = B;
      xLab = 'magnetic field (T)';
    end
    
    cla
    subplot(2,2,1)
    plot(x,muz/bmagn);
    axis tight
    ylabel('\mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
    
    subplot(2,2,2)
    plot(x,muz*avogadro);
    axis tight
    ylabel('\mu_{mol,z} (J T^{-1} mol^{-1})')
    
    subplot(2,2,3)
    plot(x,chizz_SI);
    ylabel('\chi_{mol,zz} (m^3 mol^{-1})');
    
    subplot(2,2,4)
    plot(x,chizz_SI.*T);
    ylabel('\chi_{mol,zz}T (K m^3 mol^{-1})');
    
    for i=1:4
      subplot(2,2,i);
      xlabel(xLab);
      xlim([min(x) max(x)]);
    end
    
  else
    cla
    Plot2D = (nFields>10);
    chizzT_SI = chizz_SI.*repmat(T(:).',nFields,1);
    if Plot2D
      subplot(2,2,1);
      surf(T,B,muz/bmagn);
      shading flat
      zlabel(' \mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
      subplot(2,2,2);
      surf(T,B,muz*avogadro);
      shading flat
      zlabel(' \mu_{mol,z} (J T^{-1} mol^{-1})')
      subplot(2,2,3);
      surf(T,B,chizz_SI);
      shading flat
      zlabel('\chi_{mol,zz} (m^3 mol^{-1})');
      subplot(2,2,4);
      surf(T,B,chizzT_SI);
      shading flat
      zlabel('\chi_{mol,zz}T (K m^3 mol^{-1})');
      for i=1:4
        subplot(2,2,i);
        xlabel('T (K)');
        ylabel('B (T)')
        xlim([min(T) max(T)]);
        ylim([min(B) max(B)]);
      end
    else
      subplot(2,2,1)
      plot(T,muz/bmagn);
      axis tight
      ylabel('\mu_z (\mu_B), \mu_{mol,z} (N_A\mu_B)')
      
      subplot(2,2,2)
      plot(T,muz*avogadro);
      axis tight
      ylabel('\mu_{mol,z} (J T^{-1}mol^{-1})')
      
      subplot(2,2,3)
      plot(T,chizz_SI);
      ylabel('\chi_{mol,zz} (m^3 mol^{-1})');
      
      subplot(2,2,4)
      plot(T,chizzT_SI);
      ylabel('\chi_{mol,zz}T (K m^3 mol^{-1})');
      
      for i=1:4
        subplot(2,2,i);
        xlabel('temperature (K)');
        xlim([min(T) max(T)]);
      end
      
    end
  end
end

% Assign output arguments
%---------------------------------------------------------------
switch (nargout)
  case 0
  case 1
    varargout = {muz/bmagn};
  case 2
    varargout = {muz/bmagn,chizz_SI};
  case 3
    varargout = {muz/bmagn,chizz_SI,chizzT_SI};
  otherwise
    varargout = cell(1,nargout);
    varargout(1:3) = {muz/bmagn,chizz_SI,chizzT_SI};
end

logmsg(1,'=end=curry========%s=================\n',datestr(now));

clear global EasySpinLogLevel
