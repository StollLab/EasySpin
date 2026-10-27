% chili    Simulation of cw EPR spectra in the slow motional regime
%
%   chili(Sys,Exp,Opt)
%   spc = chili(...)
%   [B,spc] = chili(...)
%   [nu,spc] = chili(...)
%   [B,spc,info] = chili(...)
%   [nu,spc,info] = chili(...)
%
%   Computes a slow-motion cw EPR spectrum.
%
%   Sys: spin system structure
%
%     Sys.tcorr       rotational correlation time (in seconds)
%     Sys.logtcorr    log10 of rotational correlation time (in seconds)
%     Sys.Diff        diffusion rate (s^-1)
%     Sys.logDiff     log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Sys.DiffFrame   Euler angles describing the orientation of the
%                     diffusion tensor in the molecular frame (default [0 0 0])
%     Sys.lw          vector with FWHM residual broadenings
%                     1 element:  GaussianFWHM
%                     2 elements: [GaussianFWHM LorentzianFWHM]
%                     units: mT for field sweeps, MHz for frequency sweeps
%     Sys.lwpp        peak-to-peak line widths, same format as Sys.lw
%     Sys.Exchange    spin exchange rate (microsecond^-1)
%     Sys.Potential   orientational potential coefficients
%                       [L1 M1 K1 lambda1; L2 M2 K2 lambda2; ...]
%
%    Exp: experimental parameter settings
%      mwFreq         microwave frequency, in GHz (for field sweeps)
%      Range          sweep range, [sweepmin sweepmax], in mT (for field sweeps)
%      CenterSweep    sweep range, [center sweep], in mT (for field sweeps)
%      Field          static field, in mT (for frequency sweeps)
%      mwRange        sweep range, [sweepmin sweepmax], in GHz (for frequency sweeps)
%      mwCenterSweep  sweep range, [center sweep], in GHz (for frequency sweeps)
%      nPoints        number of points
%      Harmonic       detection harmonic: 0, 1, 2
%      ModAmp         peak-to-peak modulation amplitude, in mT (field sweeps only)
%      mwPhase        detection phase (0 = absorption, pi/2 = dispersion)
%      Temperature    temperature, in K; if omitted: high-temperature limit
%
%   Opt: simulation options
%      LLMK           basis set parameters, [evenLmax oddLmax Mmax Kmax]
%      evenK          whether to use only even K values (true/false)
%      highField      whether to use the high-field approximation (true/false)  
%      pImax          maximum nuclear coherence order for basis
%      GridSize       grid size for powder simulation
%      PostConvNucs   nuclei to include perturbationally via post-convolution
%      Verbosity      0: no display, 1: show info
%      GridSymmetry   grid symmetry to use for powder simulation
%      separate       subspectra output, '' (default) or 'components'
%
%   Output:
%     B               magnetic field axis vector, in mT (for field sweeps)
%     nu              frequency axis vector, in GHz (for frequency sweeps)
%     spc             simulated spectrum, arbitrary units
%     info            structure containing details about the calculation
%
%     If no output arguments are specified, chili plots the simulated spectrum.

function varargout = chili(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker('internal'));

% Check Matlab version
warning(chkmlver);

if nargin<2 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>3, error('Too many output arguments.'); end

if nargin<3, Opt = struct; end
if isempty(Opt), Opt = struct; end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % print level
end
logmsg(Opt.Verbosity);

logmsg(1,'-- slow motion regime simulation ----------------------------------');

% Experimental settings
%-------------------------------------------------------------------------------
DefaultExp.nPoints = 1024;
DefaultExp.Harmonic = [];
DefaultExp.mwPhase = 0;
DefaultExp.Temperature = NaN;  % high-temperature limit
DefaultExp.ModAmp = 0;
DefaultExp.mwMode = 'perpendicular';
DefaultExp.lightBeam = '';
DefaultExp.Ordering = [];
DefaultExp.SampleFrame = [];
DefaultExp.SampleRotation = [];
DefaultExp.MolFrame = [];
DefaultExp.CrystalSymmetry = '';

Exp = adddefaults(Exp,DefaultExp);

% Photoselection is not supported
if ~isempty(Exp.lightBeam)
  error('Photoselection (via Exp.lightBeam) is not supported.')
end

% Crystal samples and sample rotations are not supported
if ~isempty(Exp.MolFrame) || ~isempty(Exp.CrystalSymmetry)
  error('Crystal samples (via Exp.MolFrame and Exp.CrystalSymmetry) are not supported.');
end
if ~isempty(Exp.SampleRotation)
  error('Sample rotations (via Exp.SampleRotation) are not supported.');
end
if ~isempty(Exp.SampleFrame) && ~isequal(size(Exp.SampleFrame),[1 3])
  error('Exp.SampleFrame must contain a single orientation (three Euler angles).');
end

% Temperature (NaN: high-temperature limit)
if ~p_temperature(Exp)
  Exp.Temperature = NaN;
end

logmsg(1,'Experiment:');

% Field or frequency sweep, mw frequency, field, number of points
allowZeroField = false;
FieldSweep = p_sweeptype(Exp,allowZeroField);
Exp.FrequencySweep = ~FieldSweep;

% Sweep range from CenterSweep or Range (CenterSweep has precedence),
% stored in Exp.Range or Exp.mwRange ([] for automatic range, which is
% available for field sweeps only)
if FieldSweep
  Exp.Range = p_sweeprange(Exp,false,false);
  Exp.CenterSweep = [];
  SweepAutoRange = isempty(Exp.Range);
else
  Exp.mwRange = p_sweeprange(Exp,true,false);
  Exp.mwCenterSweep = [];
  if isempty(Exp.mwRange)
    error('Either Exp.mwRange or Exp.mwCenterSweep need to be given.');
  end
  SweepAutoRange = false;
end

% Microwave phase: unlike pepper and garlic, Exp.mwPhase is not flipped for
% frequency sweeps, since chili_single applies it to the complex spectrum,
% which already has the correct dispersion sign along the frequency axis

% Resonator mode
switch Exp.mwMode
  case 'perpendicular', ParallelMode = false;
  case 'parallel', ParallelMode = true;
  otherwise, error('Exp.mwMode must be either ''perpendicular'' or ''parallel''.');
end
logmsg(1,'  %s mode',Exp.mwMode);
if ParallelMode
  error('chili does not support parallel-mode spectra.');
end

logmsg(1,'  %d points',Exp.nPoints);

% Director orientation: a single one given in Exp.SampleFrame (aligned
% sample), or orientationally disordered (Exp.SampleFrame empty). This is
% stored before p_sampletype, which sets an empty Exp.SampleFrame to [0 0 0].
Opt.alignedSample = ~isempty(Exp.SampleFrame);

% Sample type and partial ordering of the director (Exp.Ordering)
[Exp,Opt] = p_sampletype(Exp,Opt);

% Options
%-------------------------------------------------------------------------------
% Documented
% (Opt.LLMK default depends on the spin system, see chili_single)
if ~isfield(Opt,'highField'), Opt.highField = false; end
if ~isfield(Opt,'pImax'), Opt.pImax = []; end
if ~isfield(Opt,'pImaxall'), Opt.pImaxall = []; end
if ~isfield(Opt,'GridSize'), Opt.GridSize = [19 0]; end
if ~isfield(Opt,'LiouvMethod'), Opt.LiouvMethod = ''; end
if ~isfield(Opt,'FieldSweepMethod'), Opt.FieldSweepMethod = []; end
if ~isfield(Opt,'PostConvNucs'), Opt.PostConvNucs = ''; end
if ~isfield(Opt,'Solver'), Opt.Solver = ''; end
if ~isfield(Opt,'GridSymmetry'), Opt.GridSymmetry = 'Dinfh'; end
% Opt.Verbosity

% Common typo
if isfield(Opt,'LLKM')
  error('Opt.LLKM is not used by chili(), provide Opt.LLMK instead.');
end

% Obsolete options
if isfield(Opt,'nKnots')
  error('Options.nKnots is obsolete. Use Options.GridSize instead, e.g. Options.GridSize = 91.');
end
if isfield(Opt,'Symmetry')
  error('Options.Symmetry is obsolete. Use Options.GridSymmetry instead, e.g. Options.GridSymmetry = ''D2h''.');
end
if isfield(Opt,'Output')
  error('Options.Output is obsolete. Use Opt.separate instead.');
end

% Undocumented
if ~isfield(Opt,'Rescale'), Opt.Rescale = true; end
if ~isfield(Opt,'Threshold'), Opt.Threshold = 1e-6; end
if ~isfield(Opt,'LanczosBreakdownTol'), Opt.LanczosBreakdownTol = 1e-5; end
if ~isfield(Opt,'Lentz'), Opt.Lentz = true; end
if ~isfield(Opt,'IncludeNZI'), Opt.IncludeNZI = true; end
if ~isfield(Opt,'pqOrder'), Opt.pqOrder = false; end
if ~isfield(Opt,'GridFrame'), Opt.GridFrame = []; end
if ~isfield(Opt,'Diagnostics'), Opt.Diagnostics = ''; end
if ~isfield(Opt,'useLMKbasis'), Opt.useLMKbasis = false; end
if ~isfield(Opt,'useStartvecSelectionRules'), Opt.useStartvecSelectionRules = true; end
if ~isfield(Opt,'PeqTol'), Opt.PeqTol = []; end
if ~isfield(Opt,'jKmin'), Opt.jKmin = []; end
if ~isfield(Opt,'evenK'), Opt.evenK = []; end
if ~isfield(Opt,'BasisAnalysis'), Opt.BasisAnalysis = false; end

if isfield(Opt,'deltaK')
  error('Opt.deltaK is obsolete. Use Opt.evenK instead.');
end

if ~ischar(Opt.Diagnostics) && ~isempty(Opt.Diagnostics) && ~isvarname(Opt.Diagnostics)
  error('If given, Opt.Diagnosics must be a valid Matlab variable name.');
end

if numel(Opt.GridSize)<1, Opt.GridSize(1) = 19; end
if numel(Opt.GridSize)<2, Opt.GridSize(2) = 0; end
if Opt.GridSize(2)~=0
  error('chili cannot interpolate orientations. Set Opt.GridSize(2) to zero.');
end

% Basis settings
%-------------------------------------------------------------------------------
if isfield(Opt,'LLMK')
  if ~isnumeric(Opt.LLMK) || numel(Opt.LLMK)~=4 || any(Opt.LLMK<0) || any(mod(Opt.LLMK,1))
    error('Opt.LLMK must be a 4-element array with non-negative integers.');
  end
  if mod(Opt.LLMK(1),2)
    error('Opt.LLMK(1) must be an even number.');
  end
  if Opt.LLMK(2)>0 && ~mod(Opt.LLMK(2),2)
    error('Opt.LLMK(2) must be an odd integer.');
  end
  maxL = max(Opt.LLMK(1:2));
  if any(Opt.LLMK(3)>maxL)
    error('Opt.LLMK(3)=%d, but must be smaller or equal to the maximum L (%d).',Opt.LLMK(3),maxL);
  end
  if any(Opt.LLMK(4)>maxL)
    error('Opt.LLMK(4)=%d, but must be smaller or equal to the maximum L (%d).',Opt.LLMK(4),maxL);
  end
end

% high-field approximation
if isfield(Opt,'pSmin')
  error('Opt.pSmin is not supported - use Opt.highField instead.');
end

% pImax (maximum nuclear coherence order)
if ~isempty(Opt.pImax)
  if any(Opt.pImax<0)
    error('Every element in Opt.pImax must be 0 or larger.');
  end
end
if ~isempty(Opt.pImaxall)
  if numel(Opt.pImaxall)~=1
    error('Opt.pImaxall must be a single non-negative number.');
  end
end

% M-pS-pI symmetry (see Meirovitch J.Chem.Phys. 77 3915 (1982), eq. (A47))
if ~isfield(Opt,'MpSymm')
  Opt.MpSymm = false;
end


% Field sweep and linear solver
%-------------------------------------------------------------------------------
logmsg(1,'Methods:');

% Field sweep method: set default
if FieldSweep
  if isempty(Opt.FieldSweepMethod)
    Opt.FieldSweepMethod = 'approxlin';
  end
  [~,err] = parseoption(Opt,'FieldSweepMethod',{'explicit','approxinv','approxlin'});
  error(err);
end
explicitFieldSweep = FieldSweep && strcmp(Opt.FieldSweepMethod,'explicit');
logmsg(1,'  field sweep method: %s',Opt.FieldSweepMethod);

% Set solver if not given
if isempty(Opt.Solver)
  if strcmp(Opt.FieldSweepMethod,'explicit')
    Opt.Solver = '\';
  else
    Opt.Solver = 'L';
  end
end

if explicitFieldSweep && Opt.Solver~='\'
  error('For an explicit field sweep, use Opt.Solver=''\''.');
end

switch Opt.Solver
  case 'L'
    if Opt.Lentz==1 % Lentz method
      SolverString = 'Lanczos tridiagonalization, left-to-right continued fraction evaluation';
    else
      SolverString = 'Lanczos tridiagonalization, right-to-left continued fraction evaluation';
    end
  case 'A'
    SolverString = 'Arnoldi iteration';
  case '\'
    SolverString = 'backslash';
  case 'E'
    SolverString = 'eigenvalue method, sum of Lorentzians';
  case 'C'
    SolverString = 'conjugate gradients tridiagonalization, right-to-left continued fraction evaluation';
  case 'B'
    SolverString = 'biconjugate gradients, stabilized';
  otherwise
    error('Unknown method in Options.Solver. Must be ''L'', ''A'', ''\'', ''E'', ''C'', or ''B''.');
end
logmsg(1,'  solver: %s = %s',Opt.Solver,SolverString);

% Isotopologues and separate output
%-------------------------------------------------------------------------------
if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end

% Process Opt.separate
if ~isfield(Opt,'separate'), Opt.separate = ''; end
[separateOutput,err] = parseoption(Opt,'separate',{'','components','transitions','orientations','sites'});
error(err);
separateComponentSpectra = separateOutput==2;
separateTransitionSpectra = separateOutput==3;
separateOrientationSpectra = separateOutput==4;
separateSiteSpectra = separateOutput==5;
if separateTransitionSpectra || separateSiteSpectra || separateOrientationSpectra
  separateComponentSpectra = true;
end
if separateTransitionSpectra
  error('chili does not support Opt.separate=''transitions''.')
end
if separateSiteSpectra
  error('chili does not support Opt.separate=''sites''.')
end
if separateOrientationSpectra
  error('chili does not support Opt.separate=''orientations''.')
end


%===============================================================================
% Loop over components and isotopologues
%===============================================================================
thirdOutput = nargout>=3;
[xAxis,spec,info] = compisoloop(@chili_single,Sys,Exp,Opt,SweepAutoRange,thirdOutput,separateComponentSpectra);

% Output and plotting
switch nargout
  case 0
    cwepr_plot(xAxis,spec,Exp);
  case 1
    varargout = {spec};
  case 2
    varargout = {xAxis,spec};
  case 3
    varargout = {xAxis,spec,info};
end

end


%===============================================================================
% Single-isotopologue spectrum
%===============================================================================
% Exp and Opt have been parsed by chili(); Exp.Range is empty for automatic
% range determination.
function [xAxis,outspec,info] = chili_single(Sys,Exp,Opt)

FieldSweep = ~Exp.FrequencySweep;
useTemperature = ~isnan(Exp.Temperature);
partiallyOrderedSample = Opt.partiallyOrderedSample;
saveDiagnostics = ~isempty(Opt.Diagnostics);
explicitFieldSweep = FieldSweep && strcmp(Opt.FieldSweepMethod,'explicit');
useLanczosSolver = Opt.Solver=='L';

logmsg(1,'-- component spectrum simulation ----------------------------------');

% Spin system
%-------------------------------------------------------------------------------
if ~isfield(Sys,'Nucs'), Sys.Nucs = ''; end

[Sys,err] = validatespinsys(Sys);
error(err);

% Check for limitations in spin system
if Sys.MO_present, error('chili does not support Sys.Ham* parameters.'); end
if any(Sys.L(:)), error('chili does not support Sys.L.'); end
if isfield(Sys,'nn') && any(Sys.nn(:)~=0)
  error('chili does not support nuclear-nuclear couplings (Sys.nn).');
end
if any(Sys.lwAnisoG(:)) || ~isempty(Sys.StrainData.Q)
  error('chili does not support anisotropic line widths (Sys.lwAnisoG) or strains (Sys.StrainPars).');
end

% Dynamics
%-------------------------------------------------------------------------------
% Add defaults
if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end
if ~isfield(Sys,'Potential'), Sys.Potential = []; end
if isfield(Sys,'lwpp'), Dynamics.lwpp = Sys.lwpp; end
if isfield(Sys,'lw'), Dynamics.lw = Sys.lw; end


% Heisenberg exchange
%-------------------------------------------------------------------------------
if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end
Dynamics.Exchange = Sys.Exchange;
Dynamics.Exchange = Dynamics.Exchange*2*pi*1e6; % microseconds^-1 -> rad s^-1

% Check and parse diffusion rate information
%-------------------------------------------------------------------------------
DiffFields = {'logtcorr','tcorr','Diff','logDiff'};
hasDiffFields = isfield(Sys,DiffFields);
if sum(hasDiffFields)>1
  error('Only one of Sys.tcorr, Sys.logtcorr, Sys.Diff, and Sys.logDiff is allowed.');
elseif ~any(hasDiffFields)
  error('One of Sys.tcorr, Sys.logtcorr, Sys.Diff, and Sys.logDiff is required.');
end
for k = 1:numel(DiffFields)
  if hasDiffFields(k)
    Dynamics.(DiffFields{k}) = Sys.(DiffFields{k});
  end
end

% Orientational potential
%-------------------------------------------------------------------------------
% Error on obsolete field Sys.lambda, include explicit upgrade information
if isfield(Sys,'lambda') && ~isempty(Sys.lambda)
  if numel(Sys.lambda)>4
    error('Sys.lambda must be a vector with at most 4 elements.');
  end
  lam = Sys.lambda;
  LMK = [2 0 0; 2 0 2; 4 0 0; 4 0 2];
  str = '    Sys.Potential = [';
  for p = 1:numel(lam)
    str = [str sprintf('%d %d %d %g',LMK(p,1),LMK(p,2),LMK(p,3),lam(p))];  %#ok
    if p~=numel(lam), str = [str '; ']; end  %#ok
  end
  str = [str ']'];
  error('\n  Sys.lambda is obsolete.\n  Use the following instead:\n\n%s;    % L M K lambda\n',str);
end

% Extract and organize information about potential
if ~isempty(Sys.Potential)
  if size(Sys.Potential,2)~=4
    error('Sys.Potential needs 4 entries per row (L, M, K, lambda).');
  end
  Potential.L = Sys.Potential(:,1);
  Potential.M = Sys.Potential(:,2);
  Potential.K = Sys.Potential(:,3);
  Potential.lambda = Sys.Potential(:,4);
  rmv = Potential.lambda==0;
  rmv = rmv | (Potential.L==0 & Potential.M==0 & Potential.K==0);
  if any(rmv)
    Potential.L(rmv) = [];
    Potential.M(rmv) = [];
    Potential.K(rmv) = [];
    Potential.lambda(rmv) = [];
  end
else
  Potential.L = [];
  Potential.M = [];
  Potential.K = [];
  Potential.lambda = [];
end
usePotential = ~isempty(Potential.lambda);

% Validate inputs for orientational potential
if usePotential
  if any(Potential.L<0)
    error('L values of potential coefficients must be nonnegative.');
  end
  if any(abs(Potential.K)>Potential.L)
    error('L and K values of potential coefficients do not satisfy |K|<=L.');
  end
  if any(abs(Potential.M)>Potential.L)
    error('L and M values of potential coefficients do not satisfy |M|<=L.');
  end
  if any(Potential.K<0)
    error('Only potential terms with nonnegative values of K are allowed. Terms with negative K required to render the potential real-valued are supplemented automatically.');
  end
  if any(Potential.M(Potential.K==0)<0)
    error('For potential terms with K=0, M must be nonnegative. Terms with negative M required to render the potential real-valued are supplemented automatically.');
  end
  zeroMK = Potential.K==0 & Potential.M==0;
  if any(~isreal(Potential.lambda(zeroMK)))
    error('Potential coefficients for M=K=0 must be real-valued.');
  end
end

% Check for old-style potential (L=2,4; M=0; K=0,2; real-valued lambda)
if usePotential
  oldStylePotential = ...
     all(Potential.L==2 | Potential.L==4) && ...
     all(Potential.M==0) && ...
     all(Potential.K==0 | Potential.K==2) && ...
     all(isreal(Potential.lambda));
end

% Sweep range
%-------------------------------------------------------------------------------
if FieldSweep
  Range = Exp.Range;
  if isempty(Range)
    logmsg(1,'  automatic estimation of sweep range');
    [minB,maxB] = p_resfieldspan(Sys,Exp.mwFreq);
    Stretch = 1.25;
    Range = (maxB+minB)/2 + [-1 1]/2*Stretch*max(maxB-minB,5);
  end
  CenterField = mean(Range);
  SweepWidth = diff(Range);
  Exp.Range = CenterField + [-1 1]/2*SweepWidth;
else
  CenterFreq = mean(Exp.mwRange);
  SweepWidth = diff(Exp.mwRange);
  Exp.mwRange = CenterFreq + [-1 1]/2*SweepWidth;
  CenterField = Exp.Field;
end

if FieldSweep
  logmsg(1,'  field range (mT): min %g, max %g, center %g, width %g',...
    Exp.Range(1),Exp.Range(2),CenterField,SweepWidth);
else
  logmsg(1,'  frequency range (GHz): min %g, max %g, center %g, width %g',...
    Exp.mwRange(1),Exp.mwRange(2),CenterFreq,SweepWidth);
end

% Detection harmonic and field modulation amplitude
% (T2 relaxation always broadens the slow-motion spectrum)
[Exp,derivsLeft] = p_harmonic(Exp,FieldSweep,true);

% Determine whether to average over an orientational grid
% (without potential, no orientational averaging is necessary - it's identical to a
% single-orientation sim)
integrateOverGrid = ~Opt.alignedSample && usePotential;

% Options
%-------------------------------------------------------------------------------
if ~isfield(Opt,'LLMK')
  if usePotential
    error('Sys.Potential is given. Please provide basis set information (Opt.LLMK etc).');
  else
    Opt.LLMK = [14 7 2 6];
  end
end

% Determine default method for constructing Liouvillian
if isempty(Opt.LiouvMethod)
  if (Sys.nElectrons==1) && (Sys.S==1/2) && (Sys.nNuclei<=2) && ...
      (~usePotential || oldStylePotential)
    Opt.LiouvMethod = 'fast';
  else
    Opt.LiouvMethod = 'general';
  end
end

[LiouvMethod,err] = parseoption(Opt,'LiouvMethod',{'fast','general'});
error(err);
generalLiouvillian = LiouvMethod==2;

if generalLiouvillian
  if any(Sys.Exchange~=0)
    error('Opt.LiouvMethod=''general'' does not support spin exchange (Sys.Exchange).');
  end
else
  if Sys.nElectrons>1
    error('Opt.LiouvMethod=''fast'' does not work with more than one electron spin.');
  end
  if Sys.S~=1/2
    error('Opt.LiouvMethod=''fast'' does not work for S>1/2.');
  end
  if Sys.nNuclei>2
    error('Opt.LiouvMethod=''fast'' does not work more than 2 nuclei in the spin system.');
  end
  if usePotential
    if ~oldStylePotential
      error('Opt.LiouvMethod=''fast'' does not work with this orientational potential.');
    end
  end
end

% Post-convolution nuclei
doPostConvolution = ~isempty(Opt.PostConvNucs);
if doPostConvolution
  Opt.PostConvNucs = sort(unique(Opt.PostConvNucs));
  if any(Opt.PostConvNucs<1) || any(Opt.PostConvNucs>Sys.nNuclei)
    error('Opt.PostConvNucs must contain indices of nuclei (1 to %d).',Sys.nNuclei);
  end
  nPostConvNucs = numel(Opt.PostConvNucs);
  if (Sys.nNuclei-nPostConvNucs>2) && ~generalLiouvillian
    error('Cannot have more than two nuclei for the Stochastic Liouville equation with this Opt.LiouvMethod.');
  end
  fullSys = Sys;
  Sys = nucspinrmv(Sys,Opt.PostConvNucs);
  Sys.processed = false;
  [Sys,err] = validatespinsys(Sys);
  error(err);
end

if any(Sys.n~=1)
  error('chili cannot handle systems with nuclei with Sys.n > 1 only if these nuclei are treated using post-convolution (Opt.PostConvNucs).');
end

% Basis settings (Opt.LLMK has been validated by chili())
%-------------------------------------------------------------------------------
Basis.LLMK = Opt.LLMK;
Basis.jKmin = Opt.jKmin;
Basis.evenK = Opt.evenK;

% high-field approximation
if Opt.highField
  Basis.pSmin = +1;
else
  Basis.pSmin = -1;
end

% pImax (maximum nuclear coherence order)
if ~isempty(Opt.pImax)
  if numel(Opt.pImax)~=Sys.nNuclei && numel(Opt.pImax)~=1
    error('Opt.pImax must contain one entry for every nucleus or just a single number.');
  end
end
Basis.pImax = Opt.pImax;
Basis.pImaxall = Opt.pImaxall;

% M-pS-pI symmetry
Basis.MpSymm = Opt.MpSymm;


% Set allocation block size, used in chili_lm
%-------------------------------------------------------------------------------
if ~generalLiouvillian
  blockSize = 1e6;
  minBlockSize = 1e3;
  if ~isfield(Opt,'AllocationBlockSize')
    Opt.AllocationBlockSize = blockSize;
  end
  Opt.AllocationBlockSize = Opt.AllocationBlockSize(1);
  if Opt.AllocationBlockSize<minBlockSize
    error('Opt.AllocationBlockSize = %d is too small. Increase its value to at least %d.',Opt.AllocationBlockSize,minBlockSize);
  end
  logmsg(2,'  allocating memory in blocks of %d non-zero elements',Opt.AllocationBlockSize);
end

% Precalculate spin operator matrices
%-------------------------------------------------------------------------------
if generalLiouvillian
  logmsg(1,'  construction of Liouvillian: general code');
  
  % calculate spin operators
  for iSpin = numel(Sys.Spins):-1:1
    SpinOps{iSpin,1} = sop(Sys.Spins,[iSpin,1],'sparse'); % Sx
    SpinOps{iSpin,2} = sop(Sys.Spins,[iSpin,2],'sparse'); % Sy
    SpinOps{iSpin,3} = sop(Sys.Spins,[iSpin,3],'sparse'); % Sz
  end

  SdetOp = sparse(0);
  for e = 1:Sys.nElectrons
    SdetOp = SdetOp + SpinOps{e,1}; % Sx
  end
  
else
  logmsg(1,'  Liouvillian construction: fast S=1/2 code');
  % no need to calculate spin operators for the fast code
  SpinOps = [];
end

% Calculate ISTOs and symmetry properties
%-------------------------------------------------------------------------------
includeNQI = false;
[T,F,Sys,Symmetry,isFieldDep] = magint(Sys,SpinOps,CenterField,...
                                       Opt.IncludeNZI,includeNQI,...
                                       explicitFieldSweep);

if saveDiagnostics
  diagnostics.T = T;
  diagnostics.F = F;
end
                                     
noAnisotropiesPresent = all(F.F1(:)==0) && all(F.F2(:)==0);
if noAnisotropiesPresent
  error('This is an isotropic spin system. chili cannot calculate a slow-motion spectrum.');
end

% Process diffusion tensor and linewidth
%-------------------------------------------------------------------------------
[Dynamics,err] = processdynamics(Dynamics,FieldSweep);
error(err);


% Basis
%-------------------------------------------------------------------------------
Basis = processbasis(Basis,max(Potential.K),Sys.I,Symmetry);


% Set up g values
%-------------------------------------------------------------------------------
% Set reference g value (for frequency-to-field conversion and Boltzmann)
if Sys.fullg
  idx = 1:3;
  for iElectron = 1:Sys.nElectrons
    mean_g(iElectron) = mean(eig(Sys.g(idx,:)));  %#ok
    idx = idx + 3;
  end
  gavg = mean(mean_g);
else
  gavg = mean(mean(Sys.g));
end

% Set field and frequency axes and values
%-------------------------------------------------------------------------------
% Determine B0 and nu based on sweep domain and method
if FieldSweep
  B_ = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);
  switch Opt.FieldSweepMethod
    case 'explicit'
      B0 = B_;
      nu = Exp.mwFreq;
    case 'approxinv'
      B0 = CenterField;
      nu = Exp.mwFreq*B0./B_;
    case 'approxlin'
      B0 = CenterField;
      nu = Exp.mwFreq - unitconvert(B_-B0,'mT->MHz',gavg)/1e3; % GHz
    otherwise
      error('Unknown setting ''%s'' in Opt.FieldSweepMethod.',Opt.FieldSweepMethod);
  end
else % frequency sweep
  B0 = CenterField;
  nu = linspace(Exp.mwRange(1),Exp.mwRange(2),Exp.nPoints);
end

% Adjust units
B0 = B0/1e3; % mT -> T
omega0 = complex(2*pi*nu*1e9,1/Dynamics.T2); % GHz -> rad s^-1

% Set x axis
if FieldSweep
  xAxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);  % field axis, mT
else
  xAxis = nu; % frequency axis, GHz
end


% Set up list of orientations
%===============================================================================
if integrateOverGrid
  if Opt.GridSize(1)==1
    phi = 0;
    theta = 0;
    gridWeights = 4*pi;
    Opt.GridSymmetry = 'O3';
  else
    grid = sphgrid(Opt.GridSymmetry,Opt.GridSize(1));
    Vecs = grid.vecs;
    gridWeights = grid.weights;
    % Transform vector to reference frame representation and convert to polar angles.
    if isempty(Opt.GridFrame)
      [phi,theta] = vec2ang(Vecs);
    else
      [phi,theta] = vec2ang(Opt.GridFrame*Vecs);
    end
  end
  logmsg(1,'  grid simulation with %d orientations, grid symmetry %s (director frame relative to lab frame)',numel(phi),Opt.GridSymmetry);
else
  if Opt.alignedSample
    phi = -Exp.SampleFrame(3);
    theta = -Exp.SampleFrame(2);
  else
    phi = 0;
    theta = 0;
  end
  gridWeights = 4*pi;
  logmsg(2,'  single-orientation simulation');
end
nOrientations = numel(phi);
Basis.DirTilt = any(theta~=0);

% Partial ordering for protein/macromolecule
if partiallyOrderedSample
  orderingWeights = Exp.Ordering(0,-theta,-phi);
  if any(orderingWeights<0), error('User-supplied orientation distribution gives negative values!'); end
  if all(orderingWeights==0), error('User-supplied orientation distribution is all-zero.'); end
  logmsg(2,'  orientational potential');
else
  orderingWeights = ones(1,nOrientations);
end

Weights = gridWeights.*orderingWeights;
Weights = 4*pi*Weights/sum(Weights);


% Basis set preparations
%-------------------------------------------------------------------------------
logmsg(1,'Basis set:');
logmsg(1,'  truncation settings:');
logmsg(1,'    orientational: Leven max %d, Lodd max %d, Mmax %d, Kmax %d, evenK %d, jKmin %+d',...
  Basis.LLMK(1),Basis.LLMK(2),Basis.LLMK(3),Basis.LLMK(4),Basis.evenK,Basis.jKmin);
logmsg(1,'    spin: pSmin %+d, pImax %s, pImaxall %d',Basis.pSmin,num2str(Basis.pImax),Basis.pImaxall);
logmsg(1,'    M-p symmetry: %d',Basis.MpSymm);

if generalLiouvillian
  
  % Set up basis
  if Opt.useLMKbasis
    Basis = generateoribasis(Basis,'LMK');
  else
    Basis = generateoribasis(Basis,'LjKKM');
  end
  nOriBasis = numel(Basis.L);
  nSpinBasis = Sys.nStates^2;
  logmsg(1,'  orientational basis: %d functions',nOriBasis);
  
  % Get (p,q) quantum numbers for transitions, and index vector for reordering
  % basis states from m1-m2 order (standard) to p-q order (used in fast code)
  [idxpq,pq] = pqorder(Sys.Spins);
  % Reorder detection operator if needed
  SdetOp = SdetOp(:);
  if Opt.pqOrder
    SdetOp = SdetOp(idxpq);
  end
  
  % Removing unwanted spin functions (in mm ordering)
  keep = true;
  % (1) keep only transitions with pS>=pSmin, for each electron
  for ie = 1:Sys.nElectrons
    pS = pq(:,2*ie-1);
    keep = keep & (pS>=Basis.pSmin);
  end
  % (2a) keep only transitions with |pI|<=pImax, for each nucleus
  % (2b) keep only transitions with sum(pI)<=pImaxall
  pIsum = 0;
  for in = 1:Sys.nNuclei
    pI = pq(:,2*Sys.nElectrons+2*in-1);
    pIsum = pIsum + pI;
    keep = keep & (abs(pI)<=Basis.pImax(in));
  end
  keep = keep & (abs(pIsum)<=Basis.pImaxall);
  if Opt.pqOrder
    keep = keep(idxpq);
  end
  logmsg(1,'  spin basis: %d functions (%0.1f%% of %d)',...
    sum(keep),sum(keep)/nSpinBasis*100,nSpinBasis);
  keep = repmat(keep,nOriBasis,1);
  
  % Apply M=p-1 symmetry (Meirovitch J.Chem.Phys. 77(7), 3915-3938, Eq. (A47))
  if Opt.MpSymm
    M = Basis.M;
    psum = sum(pq(:,1:2:end),2);
    keep_Mp = (psum - M.')==1; % keep only basis states with pS+pI-M == 1
    keep = keep & keep_Mp(:);
    logmsg(1,'  applying M-p symmetry: keeping %d of %d functions',sum(keep),numel(keep));
  end

  logmsg(1,'  total basis size: %d functions',sum(keep));
  
else
  
  Basis = chili_basisbuild(Basis,Sys);
  logmsg(1,'  total basis size: %d functions',numel(Basis.L));
  
end
if saveDiagnostics
  diagnostics.basis = Basis;
end

% Precalculate 3j symbols
%-------------------------------------------------------------------------------
if generalLiouvillian
  
  logmsg(1,'Precalculating 3j symbols up to Leven=%d and Lodd=%d...',...
    Basis.LLMK(1),Basis.LLMK(2));
  computeRankOne = any(F.F1(:)~=0);
  [jjj0,jjj1,jjj2] = jjjsymbol(Basis.LLMK(1),Basis.LLMK(2),computeRankOne);
  
end

% Calculate diffusion operator matrix
%-------------------------------------------------------------------------------
% Pre-calculate diffusion operator Wigner expansion coefficient
% (needed for both Opt.LiouvMethod='fast' and 'general')
if usePotential
  logmsg(1,'Calculating Wigner expansion coefficients for potential-dependent part of diffusion matrix');
  XLMK = chili_xlmk(Potential,Dynamics.R);
else
  XLMK = {};
end

if generalLiouvillian

  logmsg(1,'Calculating the diffusion matrix');
  
  % Calculate diffusion superoperator in spatial basis
  Gamma = diffsuperop(Basis,Dynamics.R,XLMK,Potential);
  % Expand to full product basis
  Gamma = kroneye(Gamma,Sys.nStates^2);
  Gamma = Gamma(keep,keep); % prune
  
  maxerr = @(x) full(max(abs(x(:))));
  imagerr = @(x) maxerr(imag(x))/maxerr(real(x));
  logmsg(1,'  imag/real = %f',imagerr(Gamma));
  
else
  
  % Gamma is calculated simultaneously with H
  
end


% Starting vector
%-------------------------------------------------------------------------------
logmsg(1,'Calculating starting vector...');
if generalLiouvillian
  if ~isfield(Opt,'StartVec') || isempty(Opt.StartVec)
    % Calculate sqrt(Peq) vector
    Opt_.useSelectionRules = Opt.useStartvecSelectionRules;
    Opt_.PeqTolerances = Opt.PeqTol;
    [sqrtPeq,nInt] = chili_eqpopvec(Basis,Potential,Opt_);
    % Set up in full product basis, then prune
    StartVector = kron(sqrtPeq,SdetOp(:)/norm(SdetOp(:)));
    StartVector = StartVector(keep);
    normPeqVec = norm(sqrtPeq)^2;
    logmsg(1,'  norm of Peq vector: %g',normPeqVec);
    if normPeqVec<0.99
      fprintf('The norm of the equilibrium population vector in this truncated basis is %0.3g.\n It should be close to 1. The basis might be too small.\n',normPeqVec);
    end
  else
    logmsg(1,'  using provided vector');
    if numel(Opt.StartVec)==sum(keep)
      StartVector = Opt.StartVec;
      nInt = [];
    else
      error('Opt.StartVec must have %d elements.',sum(keep));
    end
  end
  StartVector = StartVector/norm(StartVector);
else
  if usePotential
    % Organize potential as expected by chili_startingvector
    LMK = [2 0 0; 2 0 2; 4 0 0; 4 0 2]; % standard terms and order for fast code
    PotLMK = [Potential.L Potential.M Potential.K];
    lambda_shortlist = [0; 0; 0; 0];
    for p = 1:4
      [found,idx] = ismember(LMK(p,:),PotLMK,'rows');
      if found, lambda_shortlist(p) = Potential.lambda(idx); end
    end
  else
    lambda_shortlist = [];
  end
  StartVector = chili_startingvector(Basis,lambda_shortlist);
  nInt = [];
end

if saveDiagnostics
  diagnostics.sv = StartVector;
end

BasisSize = size(StartVector,1);
logmsg(1,'  vector size: %dx1',BasisSize);
logmsg(1,'  non-zero elements: %d/%d (%0.2f%%)',...
  nnz(StartVector),BasisSize,100*nnz(StartVector)/BasisSize);
logmsg(1,'  maxabs %g, norm %g',full(max(abs(StartVector))),norm(StartVector));
if ~isempty(nInt)
  logmsg(1,'  evaluated integrals: 1D %d, 2D %d, 3D %d',nInt(1),nInt(2),nInt(3));
end

% Prepare transformation matrix for K-symmetrization
if Opt.useLMKbasis && useLanczosSolver
  TT = ksymmetrizer(Basis); % in spatial basis
  TT = kron(TT,eye(Sys.nStates^2)); % expand to full basis, incl. spin
  TT = TT(keep,keep); % prune
end


% Loop over all orientations
%===============================================================================
spec = 0;
for iOri = 1:nOrientations
  
  logmsg(1,'Orientation %d/%d: phi = %gdeg, theta = %gdeg (weight %g)',...
    iOri,nOrientations,phi(iOri)*180/pi,theta(iOri)*180/pi,Weights(iOri));

  % Liouvillian matrix
  %-----------------------------------------------------------------------------
  logmsg(1,'Computing Liouvillian matrix...');
  if generalLiouvillian
    [Q0B,Q1B,Q2B,Q0G,Q1G,Q2G] = rbos(T,F,[phi(iOri),theta(iOri),0],isFieldDep);
    
    if Opt.pqOrder
      Q0B = Q0B(idxpq,idxpq);
      Q0G = Q0G(idxpq,idxpq);
      for k = 1:numel(Q1B)
        Q1B{k} = Q1B{k}(idxpq,idxpq);
        Q1G{k} = Q1G{k}(idxpq,idxpq);
      end
      for k = 1:numel(Q2B)
        Q2B{k} = Q2B{k}(idxpq,idxpq);
        Q2G{k} = Q2G{k}(idxpq,idxpq);
      end
    end
    if explicitFieldSweep
      if Opt.useLMKbasis
        HB = liouvhamiltonian_LMK(Basis,Q0B,Q1B,Q2B,jjj0,jjj1,jjj2);
        HG = liouvhamiltonian_LMK(Basis,Q0G,Q1G,Q2G,jjj0,jjj1,jjj2);
      else
        HB = liouvhamiltonian(Basis,Q0B,Q1B,Q2B,jjj0,jjj1,jjj2);
        HG = liouvhamiltonian(Basis,Q0G,Q1G,Q2G,jjj0,jjj1,jjj2);
      end
      HB = HB(keep,keep);
      HG = HG(keep,keep);
    else
      Q0 = Q0B+Q0G;
      if any(F.F1(:))
        for i = 1:3
          for j = 1:3
            Q1{i,j} = Q1B{i,j} + Q1G{i,j};  %#ok
          end
        end
      else
        Q1 = {};
      end
      for i = 1:5
        for j = 1:5
          Q2{i,j} = Q2B{i,j} + Q2G{i,j};  %#ok
        end
      end
      if Opt.useLMKbasis
        H = liouvhamiltonian_LMK(Basis,Q0,Q1,Q2,jjj0,jjj1,jjj2);
      else
        H = liouvhamiltonian(Basis,Q0,Q1,Q2,jjj0,jjj1,jjj2);
      end
      H = H(keep,keep);
      if saveDiagnostics
        diagnostics.Q0 = Q0;
        diagnostics.Q1 = Q1;
        diagnostics.Q2 = Q2;        
      end
    end
    
  else
    
    Sys.d2psi = wignerd(2,phi(iOri),theta(iOri),0);
    if explicitFieldSweep
      EZ0_ = Sys.EZ0;
      EZ2_ = Sys.EZ2;
      if isfield(Sys,'NZ0')
        for iNuc = 1:numel(Sys.NZ0)
          NZ0_(iNuc) = Sys.NZ0(iNuc);  %#ok
        end
      end
    end
    
  end
    
  % Loop over all field values
  %-----------------------------------------------------------------------------
  for iB = 1:numel(B0)
    logmsg(1,'Field value %d/%d: %g mT',iB,numel(B0),B0(iB));
    
    if generalLiouvillian
      
      if explicitFieldSweep
        H = B0(iB)*HB + HG;
      end
      L = 2i*pi*H + Gamma;  % Hamiltonian: Hz -> rad s^-1
      nDim = size(L,1);
      
      if nDim~=BasisSize
        error('Matrix size (%d) inconsistent with basis size (%d). Please report.',nDim,BasisSize);
      end
      if any(isnan(L))
        error('Liouvillian matrix contains NaN entries! Please report.');
      end
      
    else
      
      if explicitFieldSweep
        Sys.EZ0 = EZ0_*B0(iB);
        Sys.EZ2 = EZ2_*B0(iB);
        if isfield(Sys,'NZ0')
          for iNuc = 1:numel(Sys.NZ0)
            Sys.NZ0(iNuc) = NZ0_(iNuc)*B0(iB);
          end
        end
      end
      Sys.DirTilt = Basis.DirTilt; % used in chili_lm
      
      % Build xlk array needed by chili_lm (rearranged from XLMK)
      maxL = numel(XLMK)-1; % maxmimum L in XLMK ( = 2*L from potential)
      xlk = [];
      for L_ = 0:maxL
        xlk(L_+1,1:2*L_+1) = XLMK{L_+1}(L_+1,:);  %#ok
      end
      Dynamics.xlk = xlk;
      Dynamics.maxL = maxL;      
      Dynamics.Diff = Dynamics.R;
      
      % Call mex function to get L matrix elements
      [r,c,Vals,nDim] = chili_lm(Sys,Basis.v,Dynamics,Opt.AllocationBlockSize);
      % (chili_lm constructs r/c/Vals for -1i*H + Gamma)
      Vals = conj(Vals); % make sure L = +1i*H + Gamma (assumes Gamma is real)
      L = sparse(r,c,Vals,BasisSize,BasisSize);
      
      if saveDiagnostics && iOri==1
        % extract H and Gamma from L = 1i*H + Gamma
        % (assumes only that H and Gamma are Hermitian)
        reL = real(L);
        imL = imag(L);
        Gamma = (reL+reL.')/2 + 1i*(imL-imL.')/2;
        H = (imL+imL.')/2 + 1i*(reL.'-reL)/2;
        H = H/(2*pi); % rad s^-1 -> Hz
      end
      
    end
    
    if saveDiagnostics && iOri==1
      diagnostics.iOri = iOri;
      diagnostics.L = L;
      diagnostics.H = H;
      diagnostics.Gamma = Gamma;
    end
    
    % Apply K-symmetrization if needed to obtain complex symmetric L for Lanczos
    % algorithm. L = i*H + Gamma is complex symmetric only if both H and Gamma
    % are real-valued.
    if generalLiouvillian && Opt.useLMKbasis && useLanczosSolver
      isComplexSymmetric = isreal(H);
      %maxerr = @(A)max(abs(A(:)));
      %imagerr = @(A)maxerr(imag(A))/maxerr(real(A));
      %Himag = imagerr(H)
      %Gimag = imagerr(Gamma)
      if ~isComplexSymmetric
        L = TT'*L*TT;
        %ksymmHimag = imagerr(TT'*H*TT);
        %ksymmHimag = imagerr(TT'*Gamma*TT);
        StartVector = TT'*StartVector;
      end
    end
    
    % Rescale for numerical stability
    if Opt.Rescale
      scale = max(abs(L(:)));
      L = L/scale;
      omega = omega0/scale;
    else
      scale = 1;
      omega = omega0;
    end
    
    maxDval = max(abs(real(L(:))));
    logmsg(2,'  size: %dx%d, maxabsreal: %g',length(L),length(L),full(maxDval));
    
    maxDvalLim = 2e3;
    if maxDval>maxDvalLim
      %  error(sprintf('Numerical instability, values in diffusion matrix are too large (%g)!',maxDval));
    end
    
    logmsg(2,'  non-zero elements: %d/%d (%0.2f%%)',nnz(L),length(L).^2,100*nnz(L)/length(L)^2);
    
    
    %===========================================================================
    % Computation of the spectral function
    %===========================================================================
    switch Opt.Solver
      
      case 'L' % Lanczos method
        if generalLiouvillian && usePotential
          maxabs = @(a)max(abs(a(:)));
          isComplexSymmetric = maxabs(L-L.')/maxabs(L) < 1e-10;
          if ~isComplexSymmetric
            error('L is not complex symmetric - cannot use Lanczos method.');
          end
        end
        [thisspec,converged,~,nSteps,breakdown] = chili_lanczos(L,StartVector,-1i*omega,Opt);
        if converged
          logmsg(2,'  Lanczos converged to within %g at iteration %d/%d',...
            Opt.Threshold,nSteps,BasisSize);
        else
          % Fall back to Arnoldi, which does not break down
          if breakdown
            logmsg(1,'  Lanczos broke down at iteration %d/%d, switching to Arnoldi',...
              nSteps,BasisSize);
          else
            logmsg(1,'  Lanczos did not converge to within %g after %d iterations, switching to Arnoldi',...
              Opt.Threshold,nSteps);
          end
          [thisspec,converged,~,nSteps] = chili_arnoldi(L,StartVector,-1i*omega,Opt);
          if converged
            logmsg(2,'  Arnoldi converged to within %g at iteration %d/%d',...
              Opt.Threshold,nSteps,BasisSize);
          else
            thisspec = ones(size(omega));
            logmsg(0,'  Arnoldi iteration did not converge to within %g after %d iterations!\n  Increase Options.LLMK (current settings [%d,%d,%d,%d])',...
              Opt.Threshold,nSteps,Opt.LLMK');
          end
        end

      case 'A' % Arnoldi method
        [thisspec,converged,~,nSteps] = chili_arnoldi(L,StartVector,-1i*omega,Opt);
        if converged
          logmsg(2,'  Arnoldi converged to within %g at iteration %d/%d',...
            Opt.Threshold,nSteps,BasisSize);
        else
          thisspec = ones(size(omega));
          logmsg(0,'  Arnoldi iteration did not converge to within %g after %d iterations!\n  Increase Options.LLMK (current settings [%d,%d,%d,%d])',...
            Opt.Threshold,nSteps,Opt.LLMK');
        end
        
      case '\' % MATLAB backslash solver for sparse linear system
        
        I = speye(size(L));
        rho0 = StartVector;
        if explicitFieldSweep
          Q = L - 1i*omega*I;
          thisspec(iB) = rho0'*(Q\rho0);
        else
          for iOmega = 1:numel(omega)
            Q = L - 1i*omega(iOmega)*I;
            thisspec(iOmega) = rho0'*(Q\rho0);
          end
        end
        
      case 'E' % eigenvalue method (sum of Lorentzians)
        % see e.g. G. Binsch, J. Am. Chem. Soc. 91, 1304 (1969)
        L = full(L);
        [U,Lam] = eig(L);
        Lam = diag(Lam);
        rho0 = StartVector;
        Amplitude = (rho0'*U).'.*(U\rho0);
        thisspec = 0;
        for iPeak = 1:numel(Amplitude)
          thisspec = thisspec + Amplitude(iPeak)./(Lam(iPeak)-1i*omega);
        end
        
      case 'B' % bi-conjugate gradients stabilized
        I = speye(size(L));
        rho0 = StartVector;
        for iOmega = 1:numel(omega)
          Q = L - 1i*omega(iOmega)*I;
          [u,~] = bicgstab(Q,rho0,Opt.Threshold,nDim);
          thisspec(iOmega) = rho0'*u;
        end
        
      case 'C' % conjugate gradients
        CGshift = 1e-6 + 1e-6i;
        [~,alpha,beta,err,stepsDone] = chili_conjgrad(L,StartVector,CGshift);        
        logmsg(1,'  step %d/%d: CG converged to within %g',...
          stepsDone,BasisSize,err);
        thisspec = chili_contfracspec(-1i*omega,alpha,beta);
        
    end
    
  end % field loop
  
  spec = spec + thisspec*Weights(iOri);
  
end % orientation loop

spec = spec/scale;

% Rescale to match rigid-limit chili intensities to pepper intensities
spec = spec*1e10;
spec = spec/2; % since chili uses normalized Sx and pepper uses unnormalized Sx
% Starting vector is normalized Sx, so scale by Tr(Sx^2)/N of the spin system
% (N = number of states), relative to S=1/2, to match pepper for any spin system
SxNormSq = sum(Sys.S.*(Sys.S+1))/3;  % Tr(Sx^2)/N
SxNormSq_S12 = (1/2)*(1/2+1)/3;  % Tr(Sx^2)/N for S=1/2
spec = spec*SxNormSq/SxNormSq_S12;

if Opt.highField
  spec = spec/2;
end

if Exp.FrequencySweep
  spec = spec*1e3;
end

% Save structure with internal data to workspace for diagnostics
if saveDiagnostics
  assignin('base',Opt.Diagnostics,diagnostics);
end

%===============================================================================



%===============================================================================
% Phasing
%===============================================================================
spec = real(exp(1i*Exp.mwPhase)*spec);


%===============================================================================
% Post-convolution
%===============================================================================
if doPostConvolution
  logmsg(1,'Postconvolution...');
  
  % Spin system with shf nuclei only
  pcidx = Opt.PostConvNucs;
  pcSys.g = mean(fullSys.g(:));
  pcSys.A = mean(fullSys.A(pcidx,:),2);
  if isfield(fullSys,'n')
    pcSys.n = fullSys.n(pcidx);
  end
  pcSys.Nucs = nuclist2string(fullSys.Nucs(pcidx));
  
  % Experimental parameters for isotropic shf spectrum
  if FieldSweep
    pcExp.Range = Exp.Range;
    pcExp.mwFreq = unitconvert(mean(Exp.Range),'mT->MHz',pcSys.g)/1e3; % GHz
    Range = Exp.Range;
  else
    pcExp.mwRange = Exp.mwRange;
    pcExp.Field = Exp.Field;
    Range = Exp.mwRange*1e3;
  end
  pcExp.Harmonic = 0;
  pcExp.nPoints = Exp.nPoints;
  
  % Linewidth for shf spectrum
  dx = diff(Range)/(Exp.nPoints-1); % field sweep: mT; freq sweep: MHz
  pcSys.lw = dx/12;
  
  % Simulate isotropic spectrum of shf nuclei
  spec_pc = garlic(pcSys,pcExp);
  spec_pc = spec_pc/sum(spec_pc);
  
  % Convolute SLE spectrum with isotropic spectrum
  spec = conv(spec,spec_pc,'same');
end


%===============================================================================
% Basis set analysis
%===============================================================================
if Opt.BasisAnalysis
  logmsg(1,'-------------------------------------------------------------------');
  logmsg(1,'Basis set analysis');
  omega_ = linspace(omega(1),omega(end),12);
  u_sum = 0;
  for iOmega = 1:numel(omega_)
    u = bicgstab(L+omega_(iOmega)*speye(size(L)),StartVector,1e-7,180);
    u_sum = u_sum + abs(u)/abs(StartVector'*u);
  end
  u_sum = u_sum/max(u_sum);

  Thr = [1e-3 1e-4 1e-5 1e-6 1e-8];
  for iThr = 1:numel(Thr)
    inc = (u_sum>Thr(iThr));
    LL = Indices(inc,1);
    Le = max(LL(mod(LL,2)==0));
    Lo = max(LL(mod(LL,2)~=0));
    jK = max(Indices(inc,2));
    K = max(Indices(inc,3));
    M = max(Indices(inc,4));
    fprintf('%0.2e: %3d %3d %3d %3d %3d, size %d\n',Thr(iThr),Le,Lo,jK,K,M,sum(inc));
  end
end


%===============================================================================
% Temperature: include Boltzmann equilibrium polarization
%===============================================================================
if FieldSweep
  DeltaE = Exp.mwFreq;  % GHz
else
  DeltaE = xAxis;  % GHz
end
if useTemperature
  Polarization = tanh(planck*DeltaE*1e9/(2*boltzm*Exp.Temperature));  % two-level
else
  % no temperature: high-temperature limit, with kT replaced by h*nuRef/2
  if FieldSweep, nuRef = Exp.mwFreq; else, nuRef = 1; end  % GHz
  Polarization = DeltaE/nuRef;
end
spec = spec.*Polarization;


%===============================================================================
% Convolutional broadening, derivatives, field modulation
%===============================================================================
% Convolution with Gaussian only. Lorentzian broadening is already included
% in the slow-motion simulation via T2, and the phase was applied above.
fwhmG = Sys.lw(1);
if FieldSweep
  unitstr = 'mT';
else
  unitstr = 'GHz';
  fwhmG = fwhmG/1e3; % MHz -> GHz
end
outspec = p_lineshape(xAxis,spec,[fwhmG 0],derivsLeft,0,Exp,unitstr,0,true);


%===============================================================================
%  Final processing
%===============================================================================
info = struct;  % not implemented yet
%===============================================================================


logmsg(1,'-------------------------------------------------------------------');

end
%===============================================================================
%===============================================================================
%===============================================================================



%===============================================================================
function Basis = processbasis(Basis,maxPotentialK,I,Symmetry)

nNuclei = numel(I);

nobetatilts = Symmetry.nobetatilts;
tensorsCollinear = Symmetry.tensorsCollinear;
axialSystem = Symmetry.axialSystem;

% Spatial basis parameters: evenLmax oddLmax Kmax Mmax jKmin evenK
%-------------------------------------------------------------------------------
Basis.evenLmax = Basis.LLMK(1);
Basis.oddLmax = Basis.LLMK(2);
Basis.Mmax = Basis.LLMK(3);
Basis.Kmax = Basis.LLMK(4);

% Set jKmin = +1 if tensorial coefficients are all real. This is the
% case when all tensors (g and A) are collinear and tilted relative to
% the diffusion tensor by the angles (0,beta,0).
if isempty(Basis.jKmin)
  Basis.jKmin = -1;
  if tensorsCollinear
    Basis.jKmin = +1;
  end
end

% Use only even K if there is no beta tilt (i.e. all +1 and -1 spherical
% tensor components in the Hamiltonian are zero).
if isempty(Basis.evenK)
  Basis.evenK = nobetatilts;
end

%{
% Use only even L values (oddLmax=0) and no K values (Kmax=0)
% in case of axial magnetic tensors, axial potential, 
% and no magnetic/diffusion tilt
if axialSystem && Basis.evenK && (isempty(maxPotentialK) || (maxPotentialK==0))
  Basis.oddLmax = 0;
  Basis.Kmax = 0;
end
%}

% Spin basis truncation parameters: pSmin, pImax, pImaxall
%-------------------------------------------------------------------------------

% pImax (maximum nuclear coherence order, for each nucleus)
if nNuclei==0
  Basis.pImax = 0;
  Basis.pImaxall = 0;
else
  if isempty(Basis.pImax)
    Basis.pImax = 2*I;
  end
  Basis.pImax = min(Basis.pImax,2*I);
  if isempty(Basis.pImaxall)
    Basis.pImaxall = sum(Basis.pImax);
  end
  Basis.pImaxall = min(Basis.pImaxall,sum(2*I));
end

% Set fields for fast two-nuclei code
if nNuclei==2
  Basis.pI1max = Basis.pImax(1);
  Basis.pI2max = Basis.pImax(2);
end

% Assemble output array of basis set parameters for chili_lm
%-------------------------------------------------------------------------------
deltaK = Basis.evenK+1;
Basis.v = [...
  Basis.evenLmax Basis.oddLmax Basis.Kmax Basis.Mmax, ...
  Basis.jKmin Basis.pSmin deltaK ...
  Basis.MpSymm ...
  Basis.pImax];

end


%===============================================================================
function [Dyn,err] = processdynamics(Dyn,FieldSweep)

err = '';

% Diffusion tensor, correlation time
%-------------------------------------------------------------------------------
% convert everything (tcorr, logcorr, logDiff, Diff) to R
if isfield(Dyn,'Diff')
  % Diff given
  if any(Dyn.Diff<0)
    error('Sys.Diff cannot be negative.');
  end
  Dyn.R = Dyn.Diff;
elseif isfield(Dyn,'logDiff')
  Dyn.R = 10.^Dyn.logDiff;
elseif isfield(Dyn,'tcorr')
  Dyn.R = 1/6./Dyn.tcorr;
elseif isfield(Dyn,'logtcorr')
  if Dyn.logtcorr>=0, error('Sys.logtcorr must be negative.'); end
  Dyn.R = 1/6./10.^Dyn.logtcorr;
else
  err = sprintf('You must specify a rotational correlation time or a diffusion tensor\n(Sys.tcorr, Sys.logtcorr, Sys.Diff or Sys.logDiff).');
  return
end

% check values
if any(Dyn.R<0)
  error('Negative diffusion rate or correlation times are not possible.');
elseif any(Dyn.R>1e12)
  fprintf('Diffusion rate very fast. Simulation might not converge.\n');
elseif any(Dyn.R<1e3)
  fprintf('Diffusion rate very slow. Simulation might not converge.\n');
end

% expand to rhombic tensor
switch numel(Dyn.R)
  case 1, Dyn.R = Dyn.R([1 1 1]);
  case 2, Dyn.R= Dyn.R([1 1 2]);
  case 3 % Diff already rhombic
  otherwise
    err = 'Sys.Diff must have 1, 2 or 3 elements (isotropic, axial, rhombic).';
    return
end

% Linewidth
%-------------------------------------------------------------------------------
if isfield(Dyn,'lw')
  if numel(Dyn.lw)>1
    if FieldSweep
      LorentzFWHM = Dyn.lw(2)*28 * 1e6; % mT -> MHz -> Hz (g = 2.0006)
    else
      LorentzFWHM = Dyn.lw(2)*1e6; % MHz -> Hz
    end
  else
    LorentzFWHM = 0;
  end
  if LorentzFWHM~=0
    % Lorentzian T2 from FWHM in freq domain 1/T2 = pi*FWHM
    Dyn.T2 = 1/LorentzFWHM/pi;
  else
    Dyn.T2 = Inf;
  end
end

end
