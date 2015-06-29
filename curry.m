% curry  Computation of magnetometry data
%         (magnetic moment, susceptibility)
%
%   curry(Sys,Exp)
%   curry(Sys,Exp,Opt)
%   muz = curry(...)
%   [muz,chizz] = curry(...)
%   [muz,chizz,chizzT] = curry(...)
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
%               in cgs units (cm^3 mol^-1)
%      chizzT   chizz*T (cm^3 mol^-1 K)
%    
%    zL is the direction of the applied static magnetic field
%
%    The size of muz, chizz, chizzT is nB x nT, where nB is the number of
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
T = Exp.Temperature; % temperature, K
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
p_symandgrid;
p_crystalorientations;
Weights = Weights/4/pi;

beta = 1./T/boltzm; 

% Initialize output arrays
muz = zeros(nFields,nTemperatures);
if calculateChi
  chizz = zeros(nFields,nTemperatures);
  chizzT = zeros(nFields,nTemperatures);
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
    muz(iB,:) = muz(iB,:) + Weights(iOri)*muz_avg;
    
    if calculateChi
      h = eps^(1/3)*max(B(iB),1); % optimal step size for numerical derivative
      B1 = B(iB) + h; h = B1 - B(iB); % prevent round-off errors
      
      [V,E] = eig(H0 - (B(iB)+h)*muOpzL);
      E = diag(E) - E(1); % J
      populations = exp(-E*beta);
      if zeroTemp, populations(1) = 1; end
      
      % population-weighted average
      muz_expect = real(diag(V'*muOpzL*V));
      muz_avg2 = (muz_expect.'*populations)./sum(populations,1);
      
      chi_ = (muz_avg2-muz_avg)/h;
      chizz(iB,:) = chizz(iB,:) + Weights(iOri)*chi_;
      chizzT(iB,:) = chizz(iB,:).*T;
    end
    
  end
end

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
    ylabel('\mu_z (\mu_B), \mu_{z,mol} (N_A\mu_B)')
    title('magnetic moment, natural units');
    
    subplot(2,2,2)
    plot(x,muz*avogadro);
    axis tight
    ylabel('\mu_{z,mol} (J T^{-1} mol^{-1})')
    title('magnetic moment, SI units');
    
    subplot(2,2,3)
    plot(x,chizz*avogadro/10);
    ylabel('\chi_{mol} (cm^3 mol^{-1})');
    title('magnetic susceptibility, CGS units');
    
    subplot(2,2,4)
    plot(x,chizzT*avogadro/10);
    ylabel('\chi_{mol}T (cm^3 mol^{-1} K)');
    title('magnetic susceptibility, SI units');
    
    for i=1:4
      subplot(2,2,i);
      xlabel(xLab);
      xlim([min(x) max(x)]);
    end
    
  else
    cla
    Plot2D = (nFields>10);
    if Plot2D
      subplot(2,2,1);
      surf(T,B,muz/bmagn);
      xlabel('T (K)'); 
      ylabel('B (T)');
      zlabel(' \mu (\mu_B), \mu_{mol} (N_A\mu_B)')
      subplot(2,2,2);
      surf(T,B,muz*avogadro);
      xlabel('T (K)');
      ylabel('B (T)');
      zlabel(' \mu_{mol} (J T^{-1} mol^{-1})')
      subplot(2,2,3);
      surf(T,B,chizz*avogadro/10);
      xlabel('T (K)');
      ylabel('B (T)')
      zlabel('\chi_{mol} (cm^3 mol^{-1})');
      subplot(2,2,4);
      surf(T,B,chizzT*avogadro/10);
      xlabel('T (K)');
      ylabel('B (T)')
      zlabel('\chi_{mol}T (cm^3 mol^{-1} K)');
      for i=1:4
        subplot(2,2,i);
        xlim([min(T) max(T)]);
        ylim([min(B) max(B)]);
      end
    else
      subplot(2,2,1)
      plot(T,muz/bmagn);
      axis tight
      ylabel('\mu (\mu_B), \mu_{mol} (N_A\mu_B)')
      title('magnetic moment, natural units');
      
      subplot(2,2,2)
      plot(T,muz*avogadro);
      axis tight
      ylabel('\mu_{mol} (J T^{-1}mol^{-1})')
      title('magnetic moment, SI units');
      
      subplot(2,2,3)
      plot(T,chizz*avogadro/10);
      ylabel('\chi_{mol} (cm^3 mol^{-1})');
      title('magnetic susceptibility, CGS units');
      
      subplot(2,2,4)
      plot(T,chizzT*avogadro/10);
      ylabel('\chi_{mol}T (cm^3 mol^{-1} K)');
      title('magnetic susceptibility, SI units');
      
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
    varargout = {muz/bmagn,chizz*avogadro/10};
  case 3
    varargout = {muz/bmagn,chizz*avogadro/10,chizzT*avogadro/10};
  otherwise
    varargout = cell(1,nargout);
    varargout(1:3) = {muz/bmagn,chizz*avogadro/10,chizzT*avogadro/10};
end

logmsg(1,'=end=curry========%s=================\n',datestr(now));

clear global EasySpinLogLevel
