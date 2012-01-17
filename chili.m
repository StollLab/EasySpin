% chili    Simulation of cw EPR spectra in the slow motional regime
%
%   chili(Sys,Exp,Opt)
%   spc = chili(...)
%   [B,spc] = chili(...)
%
%   Computes the slow-motion cw EPR spectrum of systems with
%   one electron and one nuclear spin.
%
%   Input:
%     Sys                 spin system structure
%
%     Sys.tcorr           rotational correlation time (in seconds)
%     Sys.logtcorr        log10 of rotational correlation time (in seconds)
%     Sys.Diff            diffusion rate (s^-1)
%     Sys.logDiff         log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Sys.Diffpa          Euler angles of the diffusion tensor (default [0 0 0])
%     Sys.lw              vector with FWHM residual broadenings [mT]
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
%     Sys.lwpp            peak-to-peak line widths [mT], same format as Sys.lw
%     Sys.lambda          ordering potential coefficients
%                         [lambda20 lambda22 lambda40 lambda42 lambda44]
%     Sys.Exchange        Heisenberg exchange frequency (in MHz)
%
%     Exp.mwFreq          spectrometer frequency, in GHz
%     Exp.CenterSweep     [centerField sweepWidth], in mT
%     Exp.Range           [minField maxField], in mT
%
%          Exp.Range is only used if Exp.CenterSweep is not given.
%          If both Exp.CenterField and Exp.Range are omitted, the
%          magnetic field range is determined automatically.
%
%     Exp.nPoints         number of points (default 1024)
%     Exp.Harmonic        detection harmonic: 0, 1, 2 (default 1)
%     Exp.MOMD            0: single orientation, 1: powder (MOMD)
%     Exp.psi             "director tilt" orientation
%
%     Opt.LLKM            basis size: [evenLmax oddLmax Kmax Mmax]
%     Opt.Verbosity       0: no display, 1: show info
%     Opt.nKnots          number of knots for powder simulation (MOMD)
%
%   Output:
%     B      magnetic field axis vector, in mT
%     spc    simulated spectrum, arbitrary units
%
%     If no output arguments are specified, chili plots the
%     simulated spectrum.

function varargout = chili(Sys,Exp,Opt)

if (nargin==0), help(mfilename); return; end

error(chkmlver);
if (nargin<2) || (nargin>3), error('Wrong number of input arguments!'); end
if (nargout<0), error('Not enough output arguments.'); end
if (nargout>2), error('Too many output arguments.'); end

if (nargin<3), Opt = struct('unused',NaN); end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

% --------License ------------------------------------------------
LicErr = 'Could not determine license.';
Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------

global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;


%==================================================================
% Loop over components, if necessary. 
if ~isfield(Sys,'singlecomponent')
  
  if isstruct(Sys), Sys = {Sys}; end

  spec = 0;
  for iSys = 1:numel(Sys)
    Sys{iSys}.singlecomponent = 1;
    [xField,spec_] = chili(Sys{iSys},Exp,Opt);
    if isfield(Sys{iSys},'weight')
      spec = spec + spec_*Sys{iSys}.weight;
    else
      spec = spec + spec_;
    end
  end
  
  % Output and plotting
  switch (nargout)
    case 0,
      cla
      if (xField(2)<10000)
        plot(xField,spec);
        xlabel('magnetic field [mT]');
      else
        plot(xField/1e3,spec);
        xlabel('magnetic field [T]');
      end
      axis tight
      ylabel('intensity [a.u.]');
      title(sprintf('%0.8g GHz',Exp.mwFreq));
    case 1,
      varargout = {spec};
    case 2,
      varargout = {xField,spec};
  end

  return

end




logmsg(1,'-- slow motion regime simulation ----------------------------------');

% Spin system
%-------------------------------------------------------------------
if isfield(Sys,'S')
  if (numel(Sys.S)~=1) || (Sys.S~=1/2)
    error('Sys.S is invalid. Only S=1/2 systems are supported.');
  end
end

% add non-interacting nucleus if none given
% (kernel cannot handle the absence of a nucleus)
if ~isfield(Sys,'Nucs'), Sys.Nucs = '1H'; Sys.A = [0 0 0]; end

out = isotopologues(Sys.Nucs);
if (out.nIso>1)
  error('chili does not support isotope mixtures. Please specify pure isotopes in Sys.Nucs.');
end

[Sys,err] = validatespinsys(Sys);
error(err);

mT2MHz = mt2mhz(1,mean(Sys.g));

if (Sys.nNuclei>1)
  [val,maxNuc] = max(max(abs(Sys.A),[],2));
  if (maxNuc~=1)
    error('Nucleus with largest hyperfine interaction must be first!');
  end
end

if Sys.fullg
  error('chili does not support 3x3 g matrices in Sys.g.');
end

% Dynamics
%-------------------------------------------------------------------
if ~isfield(Sys,'lambda'), Sys.lambda = 0; end
if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end
if ~isfield(Sys,'Diffpa'), Sys.Diffpa = [0 0 0]; end

if isfield(Sys,'logtcorr'), Dynamics.logtcorr = Sys.logtcorr; end
if isfield(Sys,'tcorr'), Dynamics.tcorr = Sys.tcorr; end
if isfield(Sys,'logDiff'), Dynamics.logDiff = Sys.logDiff; end
if isfield(Sys,'Diff'), Dynamics.Diff = Sys.Diff; end
if isfield(Sys,'lwpp'), Dynamics.lwpp = Sys.lwpp; end
if isfield(Sys,'lw'), Dynamics.lw = Sys.lw; end
if isfield(Sys,'lambda'), Dynamics.lambda = Sys.lambda; end
if isfield(Sys,'Exchange'), Dynamics.Exchange = Sys.Exchange; end


% Experimental settings
%-------------------------------------------------------------------
if ~isfield(Exp,'nPoints'), Exp.nPoints = 1024; end
if ~isfield(Exp,'Harmonic'), Exp.Harmonic = 1; end
if ~isfield(Exp,'mwPhase'), Exp.mwPhase = 0; end
if ~isfield(Exp,'MOMD'), Exp.MOMD = 0; end
if ~isfield(Exp,'psi'), Exp.psi = 0; end
if ~isfield(Exp,'Temperature'), Exp.Temperature = NaN; end
if ~isfield(Exp,'ModAmp'), Exp.ModAmp = 0; end
if ~isfield(Exp,'Mode'), Exp.Mode = 'perpendicular'; end


% Number of points
if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
  error('Problem with Exp.nPoints. Needs to be a number >= 2.')
end


% Resonator mode
if isfield(Exp,'Detection')
  error('Exp.Detection is obsolete. Use Exp.Mode instead.');
end
switch Exp.Mode
  case 'perpendicular', ParallelMode = 0;
  case 'parallel', ParallelMode = 1;
  otherwise, error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
end
logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.Mode);


% Temperature
if ~isnan(Exp.Temperature)
  if (num(Exp.Temperature)~=1) || isinf(Exp.Temperature) || (Exp.Temperature<0)
    error('If given, Exp.Temperature must be a positive value.')
  end
end

% Powder
if (Exp.MOMD) && isempty(Sys.lambda)
  logmsg(0,'  No ordering potential given, skipping MOMD.');
  Exp.MOMD = 0;
end

% Field modulation
if (Exp.ModAmp>0)
  logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
  if (Exp.Harmonic<1)
    error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
  end
end


% Microwave frequency
if ~isfield(Exp,'mwFreq')
  error('Please supply the microwave frequency in Exp.mwFreq.');
end
if (numel(Exp.mwFreq)~=1) || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
  error('Uninterpretable microwave frequency in Exp.mwFreq.');
end
logmsg(1,'  mw frequency %0.8g GHz',Exp.mwFreq);


% Magnetic field range
if isfield(Exp,'CenterSweep')
  if isfield(Exp,'Range')
    logmsg(0,'Using Exp.CenterSweep and ignoring Exp.Range.');
  end
else
  if isfield(Exp,'Range')
    Exp.CenterSweep = [mean(Exp.Range) diff(Exp.Range)];
  else
    logmsg(1,'  automatic determination of magnetic field range');
    I = nucspin(Sys.Nucs).';
    if numel(I)>0
      Amax = max(abs(Sys.A),[],2);
      hf = sum(I.*Amax)*1e6; % MHz -> Hz
    else
      hf = 0;
    end
    minB = planck*(Exp.mwFreq*1e9 - hf)/bmagn/max(Sys.g)/1e-3;
    maxB = planck*(Exp.mwFreq*1e9 + hf)/bmagn/min(Sys.g)/1e-3;
    Exp.CenterSweep = [(maxB+minB)/2, 1.25*(maxB-minB)];
  end
end

Exp.CenterField = Exp.CenterSweep(1);
Exp.Sweep = Exp.CenterSweep(2);
Exp.Range = Exp.CenterField + [-1 1]/2*Exp.Sweep;

if any(Exp.Range<0) || diff(Exp.Range)<=0
  error('Invalid field range! Check Exp.CenterSweep or Exp.Range.');
end

logmsg(1,'  field range [mT]: min %g, max %g, center %g, width %g',...
  Exp.Range(1),Exp.Range(2),Exp.CenterField,Exp.Sweep);

% Options
%-------------------------------------------------------------------
if isempty(Opt), Opt = struct('unused',NaN); end
Opt.maxOffset = Exp.Sweep/2*mT2MHz*1e6; % Hz
if ~isfield(Opt,'Rescale'), Opt.Rescale = 1; end % rescale A before Lanczos
if ~isfield(Opt,'Threshold'), Opt.Threshold = 1e-6; end
if ~isfield(Opt,'Diagnostic'), Opt.Diagnostic = 0; end
if ~isfield(Opt,'Method'), Opt.Method = 'L'; end
if ~isfield(Opt,'Lentz'), Opt.Lentz = 1; end
if ~isfield(Opt,'IncludeNZI'), Opt.IncludeNZI = 1; end
if ~isfield(Opt,'Output'), Opt.Output = 'summed'; end
switch Opt.Output
  case 'summed', Opt.SeparateTransitions = 0;
  case 'separate', Opt.SeparateTransitions = 1;
  otherwise, error('Wrong setting in Options.Output.');
end

if ~isfield(Opt,'nKnots'), Opt.nKnots = [5 0]; end
if numel(Opt.nKnots)<1, Opt.nKnots(1) = 5; end
if numel(Opt.nKnots)<2, Opt.nKnots(2) = 0; end

if ~isfield(Opt,'LLKM')
  Opt.LLKM = [14 7 6 2];
end
Basis.LLKM = Opt.LLKM;

switch Opt.Method
  case 'L'
    logmsg(1,'  tridiagonalization: Lanczos algorithm');
    if Opt.Lentz==1 % Lentz method
      logmsg(1,'  continued fraction evaluation: left-to-right');
    else
      logmsg(1,'  continued fraction evaluation: right-to-left');
    end
  case 'C'
    logmsg(1,'  tridiagonalization: conjugate gradients');
    logmsg(1,'  continued fraction: right-to-left evaluation');
  case 'R'
    logmsg(1,'  Reference method: biconjugate gradients, stabilized');
  otherwise
    error('Unknown method in Options.Method. Must be ''L'', ''R'' or ''C''.');
end

if ~isfield(Opt,'Allocation')
  Opt.Allocation = 1e6; % maxElements, used in chili_liouvmatrix
end
if Opt.Allocation<1e3
  error('Opt.Allocation is too small.');
end
logmsg(2,'  allocation: %d max elements',Opt.Allocation(1));

% Process
%-------------------------------------------------------
Sys = processspinsys(Sys,Exp.CenterField);
if ~Opt.IncludeNZI, Sys.NZI0 = 0; end
[Dynamics,err] = processdynamics(Dynamics);
error(err);

logmsg(1,'Computing Xlk coefficients...');
Dynamics.xlk = chili_xlk(Dynamics);
Dynamics.maxL = size(Dynamics.xlk,1)-1;

Basis = processbasis(Sys,Basis,Dynamics);

nu = linspace(-1,1,Exp.nPoints)*Opt.maxOffset;  % Hz
%z = complex(1/(2*pi*Dynamics.T2),-nu+mwFreq); % Hz
z0 = complex(1/(Dynamics.T2),2*pi*(-nu+Exp.mwFreq*1e9)); % Hz
xField = nu/1e6/mT2MHz + Exp.CenterField;


% Set up quantum numbers for basis
%-------------------------------------------------------
logmsg(1,'Examining basis...');
[Basis.Size,Indices] = chili_basiscount(Sys,Basis,Exp.psi~=0);
logmsg(1,'  Le,Lo,K,M:  %d,%d,%d,%d',...
  Basis.LLKM(1),Basis.LLKM(2),Basis.LLKM(3),Basis.LLKM(4));
logmsg(1,'  basis size: %d',Basis.Size);


% Set up list of orientations
%=====================================================================
if (Exp.MOMD)
  if Opt.nKnots(1)==1
    psi = 0;
    GeomWeights = 4*pi;
  else
    [dummy,psi,GeomWeights] = sphgrid('Dinfh',Opt.nKnots(1));
  end
  logmsg(1,'  MOMD simulation with %d orientations',numel(psi));
else
  psi = Exp.psi;
  GeomWeights = 4*pi;
  logmsg(2,'  Single orientation simulation');
end
nOrientations = numel(psi);


% Loop over all orientations
%=====================================================================

for iOri = 1:nOrientations
  
  % Director tilt
  %-------------------------------------------------------
  Dynamics.psi = psi(iOri);
  Dynamics.d2psi = wignerd(2,[0 Dynamics.psi 0]);
  logmsg(2,'orientation %d of %d: psi = %g° (weight %g)',...
    iOri,nOrientations,psi(iOri)*180/pi,GeomWeights(iOri));

  % Starting vector
  %-------------------------------------------------------
  logmsg(1,'Computing starting vector(s)...');
  StartingVector = chili_startingvector(Sys,Basis,Dynamics,Opt);
  BasisSize = size(StartingVector,1);
  nVectors = size(StartingVector,2);
  logmsg(1,'  basis size: %d, vectors: %d',BasisSize,nVectors);
  logmsg(1,'  total non-zero elements: %d (%0.2f%%)',...
    nnz(StartingVector),100*nnz(StartingVector)/BasisSize);

  % Liouville matrix
  %-------------------------------------------------------
  logmsg(1,'Computing Liouville matrix...');
  [r,c,Vals,nRows,nElm] = ...
    chili_liouvmatrix(Sys,Basis.v,Dynamics,Opt.Allocation);
  if (nRows~=BasisSize)
    Msg = sprintf('Matrix size (%d) inconsistent with basis size (%d). Please report.',nRows,BasisSize);
    error(Msg);
  end

  z = z0;
  
  if (Opt.Rescale)
    scale = -min(imag(Vals));
    Vals = Vals/scale;
    z = z/scale;
  end
  
  idx = 1:nElm;
  A = sparse(r(idx)+1,c(idx)+1,Vals(idx),nRows,nRows);
  %StartingVector = sparse(StartingVector);

  maxDvalLim = 1e3;
  maxDval = max(abs(real(Vals)));
  logmsg(1,'  maxabs diffusion matrix: %g',maxDval);
  if maxDval>maxDvalLim
    error(sprintf('Numerical instability, values in diffusion matrix are too large (%g)!',maxDval));
  end

  logmsg(1,'  non-zero elements: %d (%0.2f%%)',nnz(A),100*nnz(A)/length(A)^2);

  
  %==============================================================
  % Computation of the spectral function
  %==============================================================
  logmsg(1,'Computing spectrum...');
  nDim = length(A);

  switch Opt.Method
    case 'L'
      for iVec = 1:nVectors
        [alpha,beta,minerr] = chili_lanczos(A,StartingVector(:,iVec),z,Opt);
        minerr = minerr(end);
        if (minerr<Opt.Threshold)
          thisspec(iVec,:) = chili_contfracspec(z,alpha,beta);
          logmsg(1,'  vector %d: converged to within %g at iteration %d/%d',...
            iVec,Opt.Threshold,numel(alpha),nDim);
        else
          thisspec = ones(size(z));
          logmsg(0,'Tridiagonalization did not converge to within %g after %d steps! Increase Options.LLKM (current settings [%d,%d,%d,%d])',...
            Opt.Threshold,nDim,Opt.LLKM');
        end
      end

    case 'C'
      CGshift = 1e-6 + 1e-6i;
      [xx,alpha,beta,err,StepsDone] = chili_conjgrad(A,StartingVector,CGshift);

      logmsg(1,'  step %d/%d: CG converged to within %g',...
        StepsDone,nDim,err);

      thisspec = chili_contfracspec(z,alpha,beta);

    case 'R'
      for iz = 1:numel(z)
        u = bicgstab(A+z(iz)*speye(size(A)),StartingVector,Opt.Threshold,nDim);
        thisspec(iz) = real(u'*StartingVector);
      end

  end

  % Phasing
  if (Exp.mwPhase~=0)
    cph = cos(Exp.mwPhase);
    sph = sin(Exp.mwPhase);
    for iTrans = 1:nVectors
      spec{iTrans}(iOri,:) = cph*real(thisspec(iTrans,:))-sph*imag(thisspec(iTrans,:));
    end
  else
    for iTrans = 1:nVectors
      spec{iTrans}(iOri,:) = real(thisspec(iTrans,:));
    end
  end

end % orientation loop
%==============================================================



%==============================================================
% Accumulation of spectra
%==============================================================
% Accumulation
logmsg(1,'Spectra accumulation (%d spectra)',nOrientations);
totalspec = zeros(nVectors,Exp.nPoints);
for iTrans = 1:nVectors
  for iOri = 1:nOrientations
    totalspec(iTrans,:) = totalspec(iTrans,:) + spec{iTrans}(iOri,:)*GeomWeights(iOri);
  end
end
spec = totalspec;

%==============================================================




%==============================================================
% Postconvolution
%==============================================================
if (numel(Sys.I)>1)
  logmsg(1,'Postconvolution...');
  shfSys.g = mean(Sys.g);
  shfSys.A = mean(Sys.A(2:end,:),2);
  if isfield(Sys,'n')
    shfSys.n = Sys.n(2:end);
  end
  shfSys.Nucs = nucstringmake(Sys.Nucs(2:end));
  shfExp.Range = Exp.Range;
  shfExp.mwFreq = mt2mhz(mean(Exp.Range),shfSys.g)/1e3; % GHz
  dx = diff(Exp.Range)/(Exp.nPoints-1);
  shfSys.lw = dx/12;
  shfExp.Harmonic = 0;
  shfExp.nPoints = Exp.nPoints;
  spec_shf = garlic(shfSys,shfExp);
  spec_shf = spec_shf/sum(spec_shf);
  spec = conv(spec,spec_shf);
  spec = spec(fix(numel(spec_shf)/2)+(1:Exp.nPoints));
end
%==============================================================




%==============================================================
% Basis set analysis
%==============================================================
Opt.BasisAnalysis = 0;
if (Opt.BasisAnalysis)
  logmsg(1,'-------------------------------------------------------------------');
  logmsg(1,'Basis set analysis');
  zz = linspace(z(1),z(end),12);
  u_sum = 0;
  for iz = 1:numel(zz)
    u = bicgstab(A+z(iz)*speye(size(A)),StartingVector,1e-7,180);
    u_sum = u_sum + abs(u)/abs(StartingVector'*u);
  end
  u_sum = u_sum/max(u_sum);

  Thr = [1e-3 1e-4 1e-5 1e-6 1e-8];
  for iThr = 1:numel(Thr)
    inc = (u_sum>Thr(iThr));
    LL = Indices.L(inc);
    Le = max(LL(mod(LL,2)==0));
    Lo = max(LL(mod(LL,2)~=0));
    jK = max(Indices.jK(inc));
    K = max(Indices.K(inc));
    M = max(Indices.M(inc));
    fprintf('%0.2e: %3d %3d %3d %3d %3d, size %d\n',Thr(iThr),Le,Lo,jK,K,M,sum(inc));
  end
end
%==============================================================


% Temperature: include Boltzmann equilibrium polarization
if isfinite(Exp.Temperature)
  e = exp(-planck*Exp.mwFreq*1e9/boltzm/Exp.Temperature);
  Population = [1 e];
  Population = Population/sum(Population);
  Polarization = Population(1) - Population(2);
  spec = spec*Polarization;
end

% Parallel mode: no intensities
if ParallelMode
  spec = spec*0;
end


%==============================================================
% Convolutional broadening
%---------------------------------------------------------
%
GaussianFWHM = Sys.lw(1);
dx = xField(2) - xField(1);
if (GaussianFWHM>0)
  logmsg(1,'Convoluting with Gaussian (FWHM %g mT)...',GaussianFWHM);
  spec = convspec(spec,dx,GaussianFWHM,0,1);
end

% Lorentzian broadening is already included in the slow-motion
% simulation.

outspec = spec;
%==============================================================


%==============================================================
% Field modulation, derivatives, etc
%--------------------------------------------------------------
if (Exp.ModAmp>0)
  logmsg(1,'  applying field modulation');
  outspec = fieldmod(xField,outspec,Exp.ModAmp,Exp.Harmonic);
else
  if (Exp.Harmonic>0), outspec = deriv(xField,outspec.').'; end
  if (Exp.Harmonic>1), outspec = deriv(xField,outspec.').'; end
end
%==============================================================



%==============================================================
%  Final processing
%==============================================================

switch (nargout)
case 0,
  cla
  if (xField(2)<10000)
    plot(xField,outspec);
    xlabel('magnetic field [mT]');
  else
    plot(xField/1e3,outspec);
    xlabel('magnetic field [T]');
  end
  axis tight
  ylabel('intensity [a.u.]');
  title(sprintf('%0.8g GHz',Exp.mwFreq));
case 1,
  varargout = {outspec};
case 2,
  varargout = {xField,outspec};
end
%==============================================================


logmsg(1,'-------------------------------------------------------------------');

clear global EasySpinLogLevel

return
%====================================================================
%====================================================================
%====================================================================

%--------------------------------------------------------------------
function Sys = processspinsys(Sy,Field)

Sys = Sy;
Sys.I = nucspin(Sys.Nucs);
Sys.gn = nucgval(Sys.Nucs);
nNucs = numel(Sys.I);

% Transformation form A/g tensor frame to diffusion frame
RDiff = wignerd(2,Sys.Diffpa);

% g tensor
%------------------------------------
%Compute spherical tensor component coefficients in eigenframe
[T0,T2] = cart2spher_transforms;
Rg = wignerd(2,Sys.gpa);
T2 = RDiff*Rg*T2;
Sys.g2 = T2*Sys.g(:);
Sys.g0 = T0*Sys.g(:);

Sys.g_axial = Sys.g(1)==Sys.g(2);

if (nNucs>0)
  [Sys.A0, Sys.A2] = isto(Sys.A(1,:));
  Sys.A_axial = Sys.A2(1)==0;
  Sys.gn0 = -sqrt(1/3)*(3*Sys.gn(1));
end

if (nNucs>0)
  if ~isfield(Sys,'Apa'), Sys.Apa = [0 0 0]; end
  RA = wignerd(2,Sys.Apa(1,:));
  Sys.A2 = RDiff*RA*Sys.A2;
end

% Convert all tensorial coefficients to units of Hz
% - Electron Zeeman (rank 0 and 2)
Sys.EZ0 = bmagn*Field/1e3*Sys.g0/planck;
Sys.EZ2 = bmagn*Field/1e3*Sys.g2/planck;
if (nNucs>0)
  % Nuclear Zeeman (rank 0)
  Sys.NZ0 = nmagn*Field/1e3*Sys.gn0/(planck);
  % Hyperfine (rank 0 and 2)
  Sys.HF0 = Sys.A0*1e6;
  Sys.HF2 = Sys.A2*1e6;
else
  Sys.NZ0 = 0;
  Sys.HF0 = 0;
  Sys.HF2 = 0;
end

% Frequency (Hz) -> angular frequency
Sys.EZ0 = 2*pi*Sys.EZ0;
Sys.EZ2 = 2*pi*Sys.EZ2;
Sys.NZ0 = 2*pi*Sys.NZ0;
Sys.HF0 = 2*pi*Sys.HF0;
Sys.HF2 = 2*pi*Sys.HF2;

return

%--------------------------------------------------------------------
% isto   Computes spherical tensor components of g or A matrix,
%        given its three principal values x = [gx gy gz] or x = [Ax Ay Az].
%        F0 is the zero-rank component (one number(, F2 is a 5-vector
%        with the second-rank components.
function [F0,F2] = isto(x)
[T0,T2] = cart2spher_transforms;
F0 = T0*x(:);
F2 = T2*x(:);
return

function [T0,T2] = cart2spher_transforms
T0 = -sqrt(1/3)*[1 1 1];
T2 = [1/2 -1/2 0; 0 0 0; sqrt(2/3)*[-1/2 -1/2 1]; 0 0 0; 1/2 -1/2 0];
return

%--------------------------------------------------------------------
function Basis = processbasis(Sys,Bas,Dyn)

Basis = Bas;

Basis.evenLmax = Basis.LLKM(1);
Basis.oddLmax = Basis.LLKM(2);
Basis.Kmax = Basis.LLKM(3);
Basis.Mmax = Basis.LLKM(4);

if (Basis.oddLmax>Basis.evenLmax)
  Basis.oddLmax = Basis.evenLmax;
end

% Set jKmin = +1 if tensorial coefficients are all real. This is the
% case when g and A are collinear and the orientation of the diffusion
% tensor is described by the angles (0,beta,0)
if all(isreal(Sys.EZ2)) && all(isreal(Sys.HF2))
  Basis.jKmin = +1;
else
  Basis.jKmin = -1;
end

if (Sys.nNuclei>0)
  pImax = 2*Sys.I(1);
  if ~isfield(Basis,'pImax')
    Basis.pImax = pImax;
  end
  Basis.pImax = min(Basis.pImax,pImax);
else
  Basis.pImax = 0;
end

if ~isfield(Basis,'pSmin')
  Basis.pSmin = 0;
end

% Use only even K if there is no magnetic or diffusion tilt.
if all(Sys.EZ2([2 4])==0) && all(Sys.HF2([2 4])==0)
  Basis.deltaK = 2;
else
  Basis.deltaK = 1;
end

% Use only even L values (deltaL=2) and no K values (Kmx=0)
% in case of axial magnetic tensors, axial potential, 
% and no magnetic/diffusion tilt
if Sys.g_axial && Sys.A_axial && (Basis.deltaK==2) && (max(Dyn.KK)==0)
  Basis.deltaL = 2;
  Basis.Kmax = 0;
else
  Basis.deltaL = 1;
end


Basis.v = [Basis.evenLmax Basis.oddLmax Basis.Kmax Basis.Mmax, ...
    Basis.jKmin Basis.pSmin Basis.pImax Basis.deltaL Basis.deltaK];

return
% 
%========================================================================
function [Dyn,err] = processdynamics(D)

Dyn = D;
err = '';

% diffusion tensor, correlation time
%------------------------------------------------------------------------
% convert everything (tcorr, logcorr, logDiff) to Diff
if isfield(Dyn,'Diff')
  % Diff given
elseif isfield(Dyn,'logDiff')
  Dyn.Diff = 10.^Dyn.logDiff;
elseif isfield(Dyn,'tcorr')
  Dyn.Diff = 1/6./Dyn.tcorr;
elseif isfield(Dyn,'logtcorr')
  if Dyn.logtcorr>=0, error('Sys.logtcorr must be negative.'); end
  Dyn.Diff = 1/6./10.^Dyn.logtcorr;
else
  err = sprintf('You must specify a rotational correlation time or a diffusion tensor\n(Sys.tcorr, Sys.logtcorr, Sys.Diff or Sys.logDiff).');
  return
end

if any(Dyn.Diff<0)
  error('Negative diffusion rate or correlation times are not possible.');
elseif any(Dyn.Diff>1e12)
  fprintf('Diffusion rate very fast. Simulation might not converge.\n');
elseif any(Dyn.Diff<1e3)
  fprintf('Diffusion rate very slow. Simulation might not converge.\n');
end

% expand to rhombic tensor
switch numel(Dyn.Diff)
  case 1, Dyn.Diff = Dyn.Diff([1 1 1]);
  case 2, Dyn.Diff = Dyn.Diff([1 1 2]);
  case 3, % Diff already rhombic
  otherwise
    err = 'Sys.Diff must have 1, 2 or 3 elements (isotropic, axial, rhombic).';
    return
end

if isfield(Dyn,'lw')
  if numel(Dyn.lw)>1
    LorentzFWHM = Dyn.lw(2)*28 * 1e6; % mT -> MHz -> Hz
  else
    LorentzFWHM = 0;
  end
  if (LorentzFWHM~=0)
    % Lorentzian T2 from FWHM in freq domain 1/T2 = pi*FWHM
    Dyn.T2 = 1/LorentzFWHM/pi;
  else
    Dyn.T2 = inf;
  end
end

% Restricting potential
%------------------------------------------------------------------
if ~isfield(Dyn,'lambda'), Dyn.lambda = [0 0 0 0 0]; end
if numel(Dyn.lambda)<5, Dyn.lambda(5) = 0; end
if numel(Dyn.lambda)>5, err = 'Too many potential coefficients!'; return; end

% Process lambda list
Dyn.LL = [2 2 4 4 4];
Dyn.KK = [0 2 0 2 4];
LL = Dyn.LL;
KK = Dyn.KK;
LL(Dyn.lambda==0) = [];
KK(Dyn.lambda==0) = [];
Dyn.maxL = max(LL);
Dyn.maxK = max(KK);

% Heisenberg exchange
%------------------------------------------------------------------
if ~isfield(Dyn,'Exchange'), Dyn.Exchange=0; end
Dyn.Exchange = Dyn.Exchange*2*pi*1e6;

return
