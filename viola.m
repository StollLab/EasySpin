function varargout = viola(Sys,Exp,Opt)

if nargin<3, Opt = struct('none',NaN); end

if ~isfield(Exp,'nPoints'), Exp.nPoints = 1024; end
if ~isfield(Opt,'nKnots'), Opt.nKnots = 1; end
if ~isfield(Exp,'MOMD'), Exp.MOMD = 0; end
if Exp.MOMD, Opt.nKnots = 5; else Opt.nKnots = 1; end

% Tensor forms: 1 cartesian 2 spherical
Pars.W.form = 1; % Lorentzian broadening tensor
Pars.g.form = 1; % g tensor
Pars.A.form = 1; % hyperfine tensor
Pars.R.form = 1; % diffusion tensor

% Line broadening
%--------------------------------------------------------------------------
if isfield(Sys,'lw'), error('Use Sys.lwpp instead of Sys.lw.'); end

if ~isfield(Sys,'lwpp'), Sys.lwpp = [0 0]; end
W = Sys.lwpp(2)*8.5;
Pars.W.x = W;
Pars.W.y = W;
Pars.W.z = W;
Pars.gib0 = Sys.lwpp(1)*10;
Pars.gib2 = 0;

% Spin Hamiltonian parameters
if isfield(Sys,'g')
  switch numel(Sys.g)
    case 3, g = Sys.g;
    case 2, g = Sys.g([1 1 2]);
    case 1, g = Sys.g([1 1 1]);
    otherwise
      error('Too many numbers in Sys.g.');
  end
else
  error('Sys.g is missing.');
end
Pars.g.x = g(1); % g tensor principal values
Pars.g.y = g(2);
Pars.g.z = g(3);

Pars.gamman = nucgval(Sys.Nucs); % Nuclear g value
Pars.A.in2 = 2*nucspin(Sys.Nucs); % Twice the nuclear spin I
A = Sys.A/2.8; % Hyperfine tensor principal values in G
Pars.A.x = A(1);
Pars.A.y = A(2);
Pars.A.z = A(3);
if ~isfield(Sys,'Apa'), pa = [0 0 0]; else pa = Sys.Apa; end
if ~isfield(Sys,'gpa'), pg = [0 0 0]; else pg = Sys.gpa; end
pag = eulang(erot(pa)*erot(pg)');
pag = pag*180/pi;
Pars.m.alpha = pag(1);
Pars.m.beta = pag(2);
Pars.m.gamma = pag(3);

% Brownian rotational diffusion tensor
%--------------------------------------------------------------------------
if ~isfield(Sys,'logDiff'), error('Need Sys.logDiff'); end
logDiff = Sys.logDiff;
if numel(logDiff)==1, logDiff = logDiff([1 1 1]); end
if numel(logDiff)==2, logDiff = logDiff([1 1 2]); end
Pars.R.x = logDiff(1);
Pars.R.y = logDiff(2);
Pars.R.z = logDiff(3);

if ~isfield(Sys,'Diffpa'), Sys.Diffpa = [0 0 0]; end
%Diffpa = Sys.Diffpa*180/pi;
Pars.d.alpha = Sys.Diffpa(1);
Pars.d.beta = Sys.Diffpa(2);
Pars.d.gamma = Sys.Diffpa(3);

% Orienting potential coefficients
%--------------------------------------------------------------------------
if ~isfield(Sys,'lambda'); Sys.lambda = zeros(1,5); end
if numel(Sys.lambda)<5, Sys.lambda(5) = 0; end
lambda = Sys.lambda;
Pars.c20 = lambda(1);
Pars.c22 = lambda(2);
Pars.c40 = lambda(3);
Pars.c42 = lambda(4);
Pars.c44 = lambda(5);

% Heisenberg exchange frequency
if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end
Pars.oss = Sys.Exchange*1e6;

% Director tilt
if ~isfield(Exp,'psi'), Exp.psi = 0; end
Pars.psi = Exp.psi;

% Non-Brownian diffusion models (not supported by EasySpin)
Pars.ipdf = 0; % Diffusion model parameter. 0 Brownian, 1 non-Brownian, 2 anisotropic viscosity
Pars.p.l = 0;
Pars.p.xy = 0;
Pars.p.zz = 0;
Pars.djf.pll = 0; % Discrete jump model: djf.prp= log10 of jump rate around x,y in sec^-1
Pars.djf.prp = 0; % djf.pll = log10 of jump rate around z in sec^-1
Pars.ist = 0; % Number of symmetry-related sites for discrete-jump motion about the z axis
Pars.model.l = 0;  % Model parameter for non-Brownian motion around z axis
Pars.model.xy = 0; % Model parameter for non-Brownian motion around x, y axes
Pars.model.zz = 0; % Model parameter for non-Brownian internal spin label motion around z axis

Pars.dc20 = 0;   % Director orienting potential coefficient (units of kT)
Pars.sigpsi = 0; % Director distribution parameter. Not currently implemented.


% Experimental parameters
if isfield(Exp,'CenterSweep')
  B0 = Exp.CenterSweep(1);
  range = Exp.CenterSweep(2);
elseif isfield(Exp,'Range')
  B0 = mean(Exp.Range);
  range = diff(Exp.Range);
else
  error('Need Exp.CenterSweep or Exp.Range.');
end
Pars.B0 = B0*10;              % Center field in G
Pars.range = range*10;        % Field range in G
Pars.freq = Exp.mwFreq;       % Spectrometer frequency in GHz
Pars.phase = 0;               % Phase in degrees, 0 absorption, 90 dispersion
Pars.nfld = Exp.nPoints;      % Number of points in spectrum
if ~isfield(Exp,'Harmonic'), Exp.Harmonic = 1; end
Pars.ideriv = Exp.Harmonic;   % Derivative order


% Basis set definitions
if ~isfield(Opt,'LLKM'), Opt.LLKM = [14 7 6 2]; end
Pars.lemx = Opt.LLKM(1);
Pars.lomx = Opt.LLKM(2);
Pars.kmn = 0;
Pars.kmx = Opt.LLKM(3);
Pars.mmn = 0;
Pars.mmx = Opt.LLKM(4);
Pars.ipnmx = 2*Pars.A.in2;
Pars.ipemn = 0;

Pars.nort = Opt.nKnots;  % Number of director orientations (>1 for MOMD)

Pars.nstep = 300;  % Number of CG steps allowed
if ~isfield(Opt,'Threshold'), Opt.Threshold = 1e-6; end
Pars.cgtol = Opt.Threshold;
Pars.shift = 0;

spc = meprf(Pars);

field = Exp.CenterSweep(1) + linspace(-1,1,Exp.nPoints)/2*Exp.CenterSweep(2);
switch nargout
  case 0
    plot(field,spc);
    axis tight;
    xlabel('magnetic field [mT]');
    ylabel('intensity [a.u.]');
  case 1
    varargout = {spc};
  case 2
    varargout = {field,spc};
end

function [spc,SPC] = meprf(PARM)
% Wrapper M-file to call eprl.dll using NLS parameter structure
%
% Usage:
%      [spc,SPC] = meprf(PARM)
%
%      PARM is a parameter structure created by eprfdef
%      spc is an array of spectrum intensity values corresponding to PARM.fld
%      SPC is matrix of spectra at each orientation for MOMD spectrum
%
fparm = zeros(42,1);
fparm(1)	=	PARM.phase;	
fparm(2)	=	PARM.gib0;
fparm(3)	=	PARM.gib2;
fparm(4)	=	PARM.W.x;	
fparm(5)	=	PARM.W.y;
fparm(6)	=	PARM.W.z;
fparm(7)	= PARM.g.x;
fparm(8)  = PARM.g.y;
fparm(9)  = PARM.g.z;
fparm(10)	=	PARM.A.x;
fparm(11)	=	PARM.A.y;
fparm(12)	=	PARM.A.z;
fparm(13)	=	PARM.R.x;	
fparm(14)	=	PARM.R.y;
fparm(15)	=	PARM.R.z;
fparm(16)	=	PARM.p.l;
fparm(17)	=	PARM.p.xy;
fparm(18)	=	PARM.p.zz;
fparm(19)	=	PARM.djf.pll;
fparm(20)	=	PARM.djf.prp;
fparm(21)	=	PARM.oss;
fparm(22)	=	PARM.psi;
fparm(23)	=	PARM.d.alpha;
fparm(24)	=	PARM.d.beta;
fparm(25)	=	PARM.d.gamma;
fparm(26)	=	PARM.m.alpha;
fparm(27)   =   PARM.m.beta;
fparm(28)   =   PARM.m.gamma;
fparm(29)	=	PARM.c20;
fparm(30)	=	PARM.c22;
fparm(31)	=	PARM.c40;
fparm(32)	=	PARM.c42;
fparm(33)	=	PARM.c44;
fparm(34)   =   PARM.dc20;
fparm(35)   =   PARM.sigpsi;
fparm(36)	=	PARM.B0;
fparm(37)   =	PARM.gamman;
fparm(38)	=	PARM.cgtol;
fparm(39)   =	real(PARM.shift);			
fparm(40)   =   imag(PARM.shift);
fparm(41)	=	PARM.range;
fparm(42)	=	PARM.freq;

%Integer parameters	
iparm = zeros(23,1);

iparm(1)	=	PARM.A.in2;
iparm(2)	=	PARM.ipdf;
iparm(3)	=	PARM.ist;
iparm(4)	=	PARM.model.l;
iparm(5)	=	PARM.model.xy;
iparm(6)	=	PARM.model.zz;
iparm(7)	=	PARM.lemx;
iparm(8)	=	PARM.lomx;
iparm(9)	=	PARM.kmn;
iparm(10)=	PARM.kmx;
iparm(11)=	PARM.mmn;
iparm(12)=	PARM.mmx;
iparm(13)=	PARM.ipnmx;
iparm(14)=  PARM.ipemn;
iparm(15)=	PARM.nort;
iparm(16)=	PARM.nstep;
iparm(17) = PARM.nfld;
iparm(18)=	PARM.ideriv;
iparm(19)=	PARM.W.form;
iparm(20)=	PARM.g.form;
iparm(21)=	PARM.A.form;
iparm(22)=	PARM.R.form;
iparm(23)=  0; % PARM.matflag; deleted
basis = [];
% if PARM.ptol > 0.0  basis=PARM.basis( PARM.weight>PARM.ptol, : ); else basis=[]; end
%
% Calculate MOMD spectrum if required
%
if PARM.nort > 1 && any( [ PARM.c20, PARM.c22, PARM.c40, PARM.c42, PARM.c44 ] )
   SPC = [];
   cospsi = 0:1/(PARM.nort-1):1;
   for cp = cospsi
      fparm(22) = acos( cp )*180/pi;
      SPC = [ SPC, eprf(fparm,iparm,basis) ];
   end 
   spc = trapz( cospsi', SPC' )';
%
% Calculate single-orientation spectrum
%
else
   spc = eprf(fparm,iparm,basis);
   %spc = eprl(fparm,iparm,basis);
   SPC = [];
end
