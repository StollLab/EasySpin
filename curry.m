% curry  Computation of magnetometry data
%         (magnetic moment, susceptibility)
%
%   curry(Sys,Exp)
%   curry(Sys,Exp,Opt)
%   mu = curry(...)
%   [mu,chi] = curry(...)
%   [mu,chi,chiT] = curry(...)
%
%    Calculates the molar magnetic moment and the molar
%    magnetic susceptibility for given fields and temperatures.
%
%  Input:
%    Sys    spin system
%    Exp    experimental parameter
%      Field:  list of field values (mT)
%      Temperature:  list of temperatures (K)
%      CrystalOrientation: crystal orientation; if absent, a powder
%         average is calculated
%      CrystalSymmetry: space group of crystal (for crystal calculations)
%    Opt    calculation options
%      nKnots       number of knots for powder average
%
%  Output:
%    mu     magnetic moment, in units of Bohr magnetons
%    chi    molar susceptibility, in cgs units (cm^3 mol^-1)
%    chiT   chi*T (cm^3 mol^-1 K)
%    
%    The size of mu, chi, chiI is nB x nT, where nB is the number of
%    field values in Exp.Field and nT is the number of temperatures in
%    Exp.Temperature.
%
%   If no output argument is given, the computed data are plotted.

function varargout = curry(Sys,Exp,Opt)

if (nargin==0), help(mfilename); return; end

if (nargin<2), Exp = struct; end
if (nargin<3), Opt = struct; end

calculateChi = (nargout==0) || (nargout>1);
doPlot = (nargout==0);

% Experimental parameters
%-------------------------------------------------
if ~isfield(Exp,'Field')
  Exp.Field = 0;
  disp('Exp.Field is missing, assuming zero field.');
end
B = Exp.Field/1e3; % magnetic field, T
nFields = numel(B);

if ~isfield(Exp,'Temperature')
  error('Exp.Temperature is missing.');
end
T = Exp.Temperature; % temperature, K
nTemperatures = numel(T);

if any(T<0)
  error('Negative temperatures are not possible.')
end
zeroTemp = T==0;

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
Exp.PowderSimulation = PowderSimulation; % for communication with resf*

% Options
%-------------------------------------------------
if ~isfield(Opt,'nKnots')
  Opt.nKnots = 20;
end
if ~isfield(Opt,'Symmetry')
  Opt.Symmetry = [];
end

% Set up Hamiltonian
%-------------------------------------------------
% zero-field Hamiltonian, and magnetic dipole moment operators
% F: MHz;  GxM,GyM,GzM: MHz/mT
useSparse = false;
if useSparse
  [H0,GxM,GyM,GzM] = sham(Sys,[],'sparse');
else
  [H0,GxM,GyM,GzM] = sham(Sys);
end
H0 = H0*1e6*planck; % MHz -> J
% magnetic dipole operators, in molecular frame
muOpxM = -GxM*1e6*1e3*planck; % MHz/mT -> J/T
muOpyM = -GyM*1e6*1e3*planck; % MHz/mT -> J/T
muOpzM = -GzM*1e6*1e3*planck; % MHz/mT -> J/T

% Set up sample orientations
%-------------------------------------------------
p_symandgrid;
p_crystalorientations;
Weights = Weights/4/pi;

%{
doPowderAverage = ~isfield(Exp,'CrystalOrientation') || isempty(Exp.CrystalOrientation);
if doPowderAverage
  [zL,weights] = sphgrid('Ci',Opt.nKnots,'c');
  weights = weights/(4*pi);
else
  R_C2L = erot(Exp.CrystalOrientation);
  zL = R_C2L(3,:).'; % field is along z axis of lab frame
  weights = 1/nSites;
end
nOrientations = numel(weights);
%}
beta = 1./T/boltzm; 
mu = zeros(nFields,nTemperatures);
if calculateChi
  chi = zeros(nFields,nTemperatures);
  chiT = zeros(nFields,nTemperatures);
end

% Calculation loop
%-------------------------------------------------
for iOri = 1:nOrientations
  [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
  % projection of magnetic moment operator onto field direction
  % (field direction along z axis of lab frame)
  muOpn = zLab_M(1)*muOpxM + zLab_M(2)*muOpyM + zLab_M(3)*muOpzM; % J/T
  
  for iB = 1:numel(Exp.Field)
    if useSparse
      [V,E] = eigs(H0 - B(iB)*muOpn);
    else
      [V,E] = eig(H0 - B(iB)*muOpn);
    end
    E = diag(E); % J
    E = E - E(1);
    mun_ = real(diag(V'*muOpn*V));
    populations = exp(-E*beta);
    if zeroTemp, populations(1) = 1; end
    mun1 = (mun_.'*populations)./sum(populations,1);
    
    mu(iB,:) = mu(iB,:) + Weights(iOri)*mun1;

    if calculateChi
      h = eps^(1/3)*max(B(iB),1); % optimal step size for numerical derivative
      B1 = B(iB) + h; h = B1 - B(iB); % prevent round-off errors
      if useSparse
        [V,E] = eigs(H0 - (B(iB)+h)*muOpn);
      else
        [V,E] = eig(H0 - (B(iB)+h)*muOpn);
      end
      E = diag(E); % J
      E = E - E(1);
      mun_ = real(diag(V'*muOpn*V));
      populations = exp(-E*beta);
      if zeroTemp, populations(1) = 1; end
      mun2 = (mun_.'*populations)./sum(populations,1);
      chi_ = (mun2-mun1)/h;
      chi(iB,:) = chi(iB,:) + Weights(iOri)*chi_;
      chiT(iB,:) = chi(iB,:).*T;
    end
    
  end
end

% Graphical plotting
%-----------------------------------------------------
if doPlot
  if (nFields==1) && (nTemperatures==1)
    % do nothing
    clf
  elseif (nFields==1)
    subplot(2,2,1)
    plot(T,mu/bmagn);
    axis tight
    xlabel('T (K)');
    ylabel('\mu (\mu_B), \mu_{mol} (N_A\mu_B)')
    
    subplot(2,2,2)
    plot(T,mu*avogadro);
    axis tight
    xlabel('T (K)');
    ylabel('\mu_{mol} (J T^{-1} mol^{-1})')
    
    subplot(2,2,3)
    plot(T,chi*avogadro/10);
    xlabel('T (K)');
    ylabel('\chi_{mol} (cm^3 mol^{-1})');
    
    subplot(2,2,4)
    plot(T,chiT*avogadro/10);
    xlabel('T (K)');
    ylabel('\chi_{mol}T (cm^3 mol^{-1} K)');
    
  elseif (nTemperatures==1)
    cla
    subplot(2,2,1)
    plot(B,mu/bmagn);
    axis tight
    xlabel('magnetic field (T)');
    ylabel('\mu (\mu_B), \mu_{mol} (N_A\mu_B)')
    
    subplot(2,2,2)
    plot(B,mu*avogadro);
    axis tight
    xlabel('magnetic field (T)');
    ylabel('\mu_{mol} (J T^{-1}mol^{-1})')
    
    subplot(2,2,3)
    plot(B,chi*avogadro/10);
    xlabel('magnetic field (T)');
    ylabel('\chi_{mol} (cm^3 mol^{-1})');
    
    subplot(2,2,4)
    plot(B,chiT*avogadro/10);
    xlabel('magnetic field (T)');
    ylabel('\chi_{mol}T (cm^3 mol^{-1} K)');
    for i=1:4
      subplot(2,2,i);
      xlim([min(B) max(B)]);
    end
    
  else
    cla
    Plot2D = (nTemperatures>10);
    if Plot2D
      subplot(2,2,1);
      surf(T,B,mu/bmagn);
      xlabel('T (K)'); 
      ylabel('B (T)');
      zlabel(' \mu (\mu_B), \mu_{mol} (N_A\mu_B)')
      subplot(2,2,2);
      surf(T,B,mu*avogadro);
      xlabel('T (K)');
      ylabel('B (T)');
      zlabel(' \mu_{mol} (J T^{-1} mol^{-1})')
      subplot(2,2,3);
      surf(T,B,chi*avogadro/10);
      xlabel('T (K)');
      ylabel('B (T)')
      zlabel('\chi_{mol} (cm^3 mol^{-1})');
      subplot(2,2,4);
      surf(T,B,chiT*avogadro/10);
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
      plot(B,mu/bmagn);
      axis tight
      xlabel('magnetic field (T)');
      ylabel('\mu (\mu_B), \mu_{mol} (N_A\mu_B)')
      
      subplot(2,2,2)
      plot(B,mu*avogadro);
      axis tight
      xlabel('magnetic field (T)');
      ylabel('\mu_{mol} (J T^{-1}mol^{-1})')
      
      subplot(2,2,3)
      plot(B,chi*avogadro/10);
      xlabel('magnetic field (T)');
      ylabel('\chi_{mol} (cm^3 mol^{-1})');
      
      subplot(2,2,4)
      plot(B,chiT*avogadro/10);
      xlabel('magnetic field (T)');
      ylabel('\chi_{mol}T (cm^3 mol^{-1} K)');
      for i=1:4
        subplot(2,2,i);
        xlim([min(B) max(B)]);
      end
    end
  end
end

% Assign output arguments
%---------------------------------------------------------------
switch (nargout)
  case 0
  case 1
    varargout = {mu/bmagn};
  case 2
    varargout = {mu/bmagn,chi*avogadro/10};
  case 3
    varargout = {mu/bmagn,chi*avogadro/10,chiT*avogadro/10};
  otherwise
    varargout = cell(1,nargout);
    varargout(1:3) = {mu/bmagn,chi*avogadro/10,chiT*avogadro/10};
end
