% sf_predefinedexperiments  Expand a predefined pulse experiment (macro)
%
%   Exp = sf_predefinedexperiments(Exp)
%
%   Expands a predefined experiment into a user-defined experiment. The
%   result is the same for all simulation engines, and saffron returns it
%   in info.Exp.
%
%   Input:
%     Exp  experiment structure, with fields
%       Sequence        name of the predefined experiment: '2pESEEM',
%                       '3pESEEM', '4pESEEM', 'HYSCORE', 'MimsENDOR',
%                       'DaviesENDOR'
%       tau, T, t1, t2  initial delays (µs)
%       dt              increment(s) (µs)
%       tp              mw pulse length(s) (µs), 0 or absent = ideal pulses;
%                       a scalar applies to all mw pulses, a vector gives
%                       one length per mw pulse
%       Range           RF range [min max] (MHz) for ENDOR
%       tprf            RF pulse length (µs) for ENDOR (default 20)
%       nPoints         number of points (default [256 256] for HYSCORE,
%                       1001 for ENDOR, 512 otherwise)
%       DetWindow       detection window (default 0)
%
%   Output:
%     Exp  experiment structure, with the macro parameters (tau, T, t1, t2,
%          dt, tp, Range, tprf) removed, and with the fields
%       Sequence        pulse sequence: cell array of pulses and delays (µs)
%       nPoints, Dim1, Dim2   incrementation scheme
%       Pathways        coherence pathway selection
%       DetWindow, DetPhase   detection window and phase
%       rfFreq          RF carrier frequency (MHz), for ENDOR
%       PhaseCycle      RF on/off cycle, for ENDOR
%       Macro           name of the predefined experiment

function Exp = sf_predefinedexperiments(Exp)

name = Exp.Sequence;
macros = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR','DaviesENDOR'};
if ~any(strcmp(name,macros))
  error('Exp.Sequence = ''%s'' is not a predefined experiment. Use one of: %s.',...
    name,strjoin(macros,', '));
end

if isfield(Exp,'PhaseCycle') && ~isempty(Exp.PhaseCycle)
  error('Exp.PhaseCycle cannot be used with predefined experiments. Use a user-defined experiment instead.');
end
if isfield(Exp,'DetPhase')
  error('Exp.DetPhase cannot be used with predefined experiments. Use a user-defined experiment instead.');
end
if isfield(Exp,'Pathways')
  error('Exp.Pathways cannot be used with predefined experiments. Use a user-defined experiment instead.');
end
for f = {'Dim1','Dim2','rfFreq'}
  if isfield(Exp,f{1})
    error('Exp.%s cannot be used with predefined experiments. Use a user-defined experiment instead.',f{1});
  end
end

isENDOR = any(strcmp(name,{'MimsENDOR','DaviesENDOR'}));

% Delays
if ~isfield(Exp,'tau')
  if any(strcmp(name,{'2pESEEM'}))
    Exp.tau = 0;
  else
    error('Exp.tau is missing.');
  end
end
if ~isscalar(Exp.tau) || Exp.tau<0
  error('Exp.tau must be a single nonnegative number (in µs).');
end
if ~strcmp(name,'2pESEEM') && Exp.tau==0
  error('Exp.tau must be larger than 0.');
end
if ~isfield(Exp,'T'), Exp.T = 0; end
if ~isfield(Exp,'t1'), Exp.t1 = 0; end
if ~isfield(Exp,'t2'), Exp.t2 = 0; end
if ~isscalar(Exp.T), error('Exp.T must be a single number (in µs).'); end
if ~isscalar(Exp.t1), error('Exp.t1 must be a single number (in µs).'); end
if ~isscalar(Exp.t2), error('Exp.t2 must be a single number (in µs).'); end

% Increments and number of points
if ~isENDOR
  if ~isfield(Exp,'dt')
    error('Exp.dt is missing.');
  end
end
if ~isfield(Exp,'nPoints')
  switch name
    case 'HYSCORE', Exp.nPoints = [256 256];
    case {'MimsENDOR','DaviesENDOR'}, Exp.nPoints = 1001;
    otherwise, Exp.nPoints = 512;
  end
end
if strcmp(name,'HYSCORE')
  Exp.nPoints = Exp.nPoints.*[1 1];
  if numel(Exp.nPoints)~=2
    error('Exp.nPoints needs 1 or 2 elements for HYSCORE.');
  end
  if numel(Exp.dt)>2
    error('Exp.dt needs either 1 or 2 elements, one per dimension. You gave %d.',numel(Exp.dt));
  end
elseif numel(Exp.nPoints)~=1
  error('Exp.nPoints must be a single number for %s.',name);
elseif isfield(Exp,'dt') && numel(Exp.dt)~=1
  error('Exp.dt must be a single value for %s. You gave %d.',name,numel(Exp.dt));
end

% mw pulses (along x)
nmwPulses = [2 3 4 4 3 3];
nmwPulses = nmwPulses(strcmp(name,macros));
if ~isfield(Exp,'tp') || isempty(Exp.tp)
  tp = zeros(1,nmwPulses);
else
  tp = Exp.tp;
  if isscalar(tp), tp = tp*ones(1,nmwPulses); end
  if numel(tp)~=nmwPulses
    error('Exp.tp must have 1 or %d elements (one per microwave pulse) for %s.',nmwPulses,name);
  end
  if any(tp<0)
    error('Exp.tp cannot be negative.');
  end
end
if strcmp(name,'DaviesENDOR') && any(tp==0)
  error(['Davies ENDOR needs finite, selective microwave pulses. Give their lengths ' ...
    'in Exp.tp (µs), e.g. Exp.tp = [0.2 0.1 0.2] for the inversion, pi/2 and pi pulses.']);
end
flips = {[pi/2 pi],[pi/2 pi/2 pi/2],[pi/2 pi/2 pi pi/2],[pi/2 pi/2 pi pi/2],...
  [pi/2 pi/2 pi/2],[pi pi/2 pi]};
flips = flips{strcmp(name,macros)};
P = cell(1,nmwPulses);
for k = 1:nmwPulses
  P{k}.Flip = flips(k);
  if tp(k)>0
    P{k}.tp = tp(k);
  end
end

% RF pulse for ENDOR
if isENDOR
  if ~isfield(Exp,'Range')
    error('Frequency range (Exp.Range, in MHz) must be given for an ENDOR experiment.');
  end
  if numel(Exp.Range)~=2 || Exp.Range(2)<=Exp.Range(1)
    error('Exp.Range must be [min max], in MHz.');
  end
  if ~isfield(Exp,'tprf'), Exp.tprf = 20; end
  RF.Channel = 'rf';
  RF.Flip = pi;
  RF.tp = Exp.tprf;
  Exp.rfFreq = Exp.Range(1);
  if Exp.nPoints>1
    rfStep = diff(Exp.Range)/(Exp.nPoints-1);
  else
    rfStep = 0;
  end
  % RF on/off difference: echo(off) - echo(on), the ENDOR effect
  onOff = [0 -1 1; 0 1 0];
end

% Sequence, indirect dimensions, pathways, detection phase. The detection
% phase makes the echo real and positive for ideal pulses: the pulse
% transfer prefactors of the echo pathways are -i/2 (2pESEEM), -i/8
% (3pESEEM, Mims ENDOR) and +i/8 (4pESEEM, HYSCORE).
switch name
  case '2pESEEM'
    Exp.Sequence = {P{1} Exp.tau P{2} Exp.tau};
    Exp.Dim1 = {'d1,d2' Exp.dt};
    Exp.Pathways = '+-';
    Exp.DetPhase = -pi/2;
  case '3pESEEM'
    Exp.Sequence = {P{1} Exp.tau P{2} Exp.T P{3} Exp.tau};
    Exp.Dim1 = {'d2' Exp.dt};
    Exp.Pathways = '+0-';
    Exp.DetPhase = -pi/2;
  case '4pESEEM'
    Exp.Sequence = {P{1} Exp.tau P{2} Exp.T P{3} Exp.T P{4} Exp.tau};
    Exp.Dim1 = {'d2,d3' Exp.dt};
    Exp.Pathways = '+00-';
    Exp.DetPhase = +pi/2;
  case 'HYSCORE'
    dt = Exp.dt.*[1 1];
    Exp.Sequence = {P{1} Exp.tau P{2} Exp.t1 P{3} Exp.t2 P{4} Exp.tau};
    Exp.Dim1 = {'d2' dt(1)};
    Exp.Dim2 = {'d3' dt(2)};
    Exp.Pathways = '+00-';
    Exp.DetPhase = +pi/2;
  case 'MimsENDOR'
    Exp.Sequence = {P{1} Exp.tau P{2} Exp.T RF 0 P{3} Exp.tau};
    Exp.Dim1 = {'rfFreq' rfStep};
    Exp.Pathways = '+00-';
    Exp.PhaseCycle = {[],[],onOff,[]};
    Exp.DetPhase = -pi/2;
  case 'DaviesENDOR'
    Exp.Sequence = {P{1} Exp.T RF 0 P{2} Exp.tau P{3} Exp.tau};
    Exp.Dim1 = {'rfFreq' rfStep};
    Exp.Pathways = '00+-';
    % the inverted spin packets give a negative echo; the RF pulse reduces
    % the inversion, so the ENDOR effect is echo(on) - echo(off)
    Exp.PhaseCycle = {[],[0 1 1; 0 -1 0],[],[]};
    Exp.DetPhase = -pi/2;
end
if ~isfield(Exp,'DetWindow') || isempty(Exp.DetWindow)
  Exp.DetWindow = 0;
end
Exp.Macro = name;

% Remove macro parameters
Exp = rmfield(Exp,intersect(fieldnames(Exp),{'tau','T','t1','t2','dt','tp','Range','tprf'}));

end
