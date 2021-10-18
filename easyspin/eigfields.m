% eigfields  Resonance fields by the eigenfield method 
%
%   B = eigfields(Sys, Par,)
%   B = eigfields(Sys, Par, Opt)
%   [B, Int] = eigfields(...)
%
%   Calculates all resonance fields of spin system
%   Sys solving a generalized eigenvalue problem
%   in Liouville space together with transition pro-
%   babilitites.
%
%   Input:
%   - Sys: spin system specification structure
%   - Par: structure with fields
%        mwFreq - spectrometer frequency [GHz]
%        Mode - 'parallel' or 'perpendicular' (default)
%          direction of mirowave field relative to static field
%        Range - [Bmin Bmax] If set, compute only eigenfields
%           between Bmin and Bmax. [mT]
%   - Opt: options structure with fields
%        Threshold - if set, return only transitions with
%          relative intensity above Threshold.
%
%   Output:
%   - B:   cell array of all resonance fields [mT]
%   - Int: transition intensities [MHz^2/mT^2]

function varargout = eigfields(SpinSystem, Exp, Opt)

if nargin==0, help(mfilename); return; end

% Uses generalised Liouville space eigenvalue problem
% formulation by Belford et al.

% Belford, Belford, Burkhalter, J.Magn.Reson. 11, 251-265 (1973)
% ATTENTION: different definition of direct product, hence change
% in the definition of A and B

% Check Matlab version.
error(chkmlver);

% Add empty Options structure if not specified.
switch nargin
  case 3
  case 2
    Opt = [];
  otherwise
    error('Incorrect number of inputs!');
end

if nargout>2, error('Incorrect number of outputs.'); end

if isempty(Opt)
  Opt = struct;
end

if ~isstruct(SpinSystem) || ~isstruct(Exp) || ~isstruct(Opt)
  error('SpinSystem, Parameters and Options must be structures!');
end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;

% Mute warnings because of unavoidable division by zero.
OldWarningState = warning('off');

% Process SpinSystem structure.
%===================================================================
[SpinSystem,err] = validatespinsys(SpinSystem);
error(err);

% Process Parameter structure.
%===================================================================
if isfield(Exp,'Detection')
  error('Exp.Detection is obsolete. Use Exp.Mode instead.');
end

DefaultParameters.mwFreq = NaN;
DefaultParameters.Range = [0 realmax];
DefaultParameters.Mode = 'perpendicular';
DefaultParameters.Temperature = NaN; % not implemented!!

DefaultParameters.CrystalSymmetry = '';
DefaultParameters.CrystalOrientation = [];
DefaultParameters.MolFrame = [];

Exp = adddefaults(Exp,DefaultParameters);

if isnan(Exp.mwFreq), error('Parameters.mwFreq missing!'); end

if (diff(Exp.Range)<=0) || any(~isfinite(Exp.Range)) || ...
   ~isreal(Exp.Range) || any(Exp.Range<0) || (numel(Exp.Range)~=2)
  error('Parameters.Range is not valid!');
end

if isempty(Exp.Mode), Exp.Mode = 'perpendicular'; end

ParallelMode = (2==parseoption(Exp,'Mode',{'perpendicular','parallel'}));

if ~isnan(Exp.Temperature)
  warning('Thermal equilibrium populations not implemented. Parameters.Temperature is ignored!');
end

mwFreq = Exp.mwFreq*1e3;

% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,~,averageOverChi] = p_crystalorientations(Exp,Opt);


% Process options structure.
%===================================================================
DefaultOptions.Freq2Field = true;
DefaultOptions.Threshold = 0;
DefaultOptions.RejectionRatio = 1e-8; % UNDOCUMENTED!

Opt = adddefaults(Opt,DefaultOptions);

if (Opt.Freq2Field~=1) && (Opt.Freq2Field~=0)
  error('Options.Freq2Field incorrect!');
end

computeFreq2Field = Opt.Freq2Field;

if ~isnumeric(Opt.Threshold) || numel(Opt.Threshold)~=1 || Opt.Threshold<0
  error('Opt.Threshold must be a single nonnegative number.');
end

% Determine whether intensities should be calculated.
computeIntensities = nargout>1 || Opt.Threshold~=0;

msg = 'positions';
if computeIntensities
  msg = [msg, ', intensities'];
end
logmsg(1,'  computing %s',msg);


% Build Hamiltonian components.
%===================================================================
if iscell(SpinSystem)
  [F,Gx,Gy,Gz] = deal(SpinSystem);
else
  [F,Gx,Gy,Gz] = sham(SpinSystem);
end

% Build Liouville space operators.
A = eyekron(F) - kroneye(conj(F)) + mwFreq*eye(length(F)^2);

% Check if is positive-definite. If yes, a simple eigenvalue
% problem has to be solved, not the general one.
E = diag(eig(A));
SimpleEigenproblem = all(E>0);
if SimpleEigenproblem
  msg = 'reduced to simple eigenproblem';
else
  msg = 'general eigenproblem';
end
logmsg(1,'  %s, matrix size %dX%d',msg,length(A),length(A));
clear V E;

% Prepare vectors for intensity computation.
if computeIntensities
  vGx = reshape(Gx.',numel(Gx),1);
  vGy = reshape(Gy.',numel(Gy),1);
  vGz = reshape(Gz.',numel(Gz),1);
end

% Loop over all orientations of the spin system.
EigenFields = cell(1,nOrientations);
Intensities = cell(1,nOrientations);

for iOri = 1:nOrientations
  logmsg(3,'  orientation %d of %d',iOri,nOrientations);
  
  [xLab,yLab,zLab] = erot(Orientations(iOri,:),'rows');
  GzL = zLab(1)*Gx + zLab(2)*Gy + zLab(3)*Gz;

  B = kroneye(conj(GzL)) - eyekron(GzL);
  if SimpleEigenproblem, BB = A\B; end

  if computeIntensities
    % Eigenvectors correspond to reshape(|u><v|,[],1) = kron(conj(v),u)
    % For a given eigenfield, u and v are the states with Ev-Eu = mwFreq
    if SimpleEigenproblem
      [Vecs,Fields] = eig(BB);
      [Fields,idx] = sort(1./diag(Fields));
      Vecs = Vecs(:,idx);
    else
      [Vecs,Fields] = eig(A,B);
      [Fields,idx] = sort(diag(Fields));
      Vecs = Vecs(:,idx);
    end

    % Remove negative, nonfinite and complex eigenfields
    % and those above user limit Options.MaxField
    idx = (abs(imag(Fields))<Opt.RejectionRatio*abs(real(Fields))) & ...
      (Fields>0) & isfinite(Fields) & ...
      (Fields>Exp.Range(1)) & (Fields<Exp.Range(2));
    if ~any(idx)
      EigenFields{iOri} = [];
      Intensities{iOri} = [];
    else
      EigenFields{iOri} = real(Fields(idx));
      Vecs = Vecs(:,idx);
      
      % Normalize eigenvectors to unity
      Norms = sqrt(sum(abs(Vecs).^2));
      Vecs = Vecs./Norms(ones(size(Vecs,1),1),:);
      
      % Compute quantum-mechanical transition rate
      idx = ones(1,length(EigenFields{iOri}));
      if (ParallelMode)
        vGzL = zLab(1)*vGx + zLab(2)*vGy + zLab(3)*vGz;
        if averageOverChi
          TransitionRate = abs(sum(vGzL(:,idx).*Vecs)).^2;
        else
          TransitionRate = abs(sum(vGzL(:,idx).*Vecs)).^2;
        end
      else
        vGxL = xLab(1)*vGx + xLab(2)*vGy + xLab(3)*vGz;
        if averageOverChi
          vGyL = yLab(1)*vGx + yLab(2)*vGy + yLab(3)*vGz;
          % Calculate transition rate using <v|A|u> = trace(A|u><v|)
          TransitionRate = (abs(sum(vGxL(:,idx).*Vecs)).^2 + abs(sum(vGyL(:,idx).*Vecs)).^2)/2;
        else
          TransitionRate = abs(sum(vGxL(:,idx).*Vecs)).^2;
        end
      end
      
      % Compute polarization
      Polarization = 1;
      Polarization = Polarization/prod(2*SpinSystem.I+1);
      
      % Compute frequency-to-field domain conversion factor
      if computeFreq2Field
        % 1/(<v|G|v>-<u|G|u>) = 1/(trace(G|v><v|) - trace(G|u><u|)) =
        %   1/trace(A*(|v><v|-|u><u|)) = 1/trace(A*commute(|u><v|,|v><u|))
        % |u><u| = (|u><v|)(|v><u|)
        n = length(F);
        Vecs = reshape(Vecs,n,n,numel(Vecs)/n^2);
        dBdE = [];
        for iVec = 1:size(Vecs,3)
          V = Vecs(:,:,iVec);
          dBdE(iVec) = 1/abs(trace(GzL*commute(V,V')));
        end
      else
        dBdE = ones(size(TransitionRate));
      end
      
      % Combine factors
      Intensities{iOri} = Polarization*real(TransitionRate.*dBdE).';
      
      idx = Intensities{iOri}>=Opt.Threshold(1)*max(Intensities{iOri});
      EigenFields{iOri} = EigenFields{iOri}(idx);
      Intensities{iOri} = Intensities{iOri}(idx);
    end
  else
    
    if SimpleEigenproblem
      Fields = eig(BB);
      Fields = sort(1./Fields);
    else
      Fields = eig(A,B);
      Fields = sort(Fields);
    end

    inRange = (abs(imag(Fields))<Opt.RejectionRatio*abs(real(Fields))) & ...
      (Fields>0) & isfinite(Fields) & ...
      (Fields>Exp.Range(1)) & (Fields<Exp.Range(2));
    EigenFields{iOri} = real(Fields(inRange));
    Intensities = {[]};
    
  end
  
end

% One orientation: simple array instead of 1x1 cell array as output!
if nOrientations==1
  EigenFields = EigenFields{1};
  Intensities = Intensities{1};
end

% Prepare output.
switch nargout
  case 1, varargout = {EigenFields};
  case 2, varargout = {EigenFields,Intensities};
end

% Restore original warning state.
warning(OldWarningState);

return
%=======================================================================
%=======================================================================
