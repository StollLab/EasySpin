% pulse      Pulse definition
%
% [t,IQ] = pulse(Par)
% [t,IQ] = pulse(Par,Opt)
%
% [t,IQ,modulation] = pulse(Par,Opt)
%
% Input:
%   Par = structure containing the following fields:
%     Par.tp          = pulse length, in µs
%     Par.TimeStep    = time step for waveform definition, in µs (default:
%                       determined automatically based on pulse parameters)
%     Par.Flip        = pulse flip angle, in radians (see Ref. 1)
%                       (default: pi)
%     Par.Amplitude   = pulse amplitude, in MHz; ignored if Par.Flip given
%                       (if Opt.Transmitter is used, the amplitude is 
%                        interpreted as the input amplitude (relative to 
%                        Par.InputAmplitude, e.g. 0.5 corresponds to a 
%                        pulse with one half of the maximum digital 
%                        amplitude of the AWG)
%     Par.Qcrit       = critical adiabaticity, used to calculate pulse
%                       amplitude for frequency-swept pulses [1]; if given
%                       takes precedence over Par.Amplitude and Par.Flip
%     Par.Frequency   = pulse frequency; center frequency for amplitude
%                       modulated pulses, [start-frequency end-frequency]
%                       for frequency swept pulses; (default: 0)
%     Par.Phase       = phase for the pulse in radians (default: 0 = +x)
%     Par.Type        = pulse shape name in a string with structure 'AM/FM'
%                       (or just 'AM'), where AM refers to the amplitude
%                       modulation function and FM to the frequency
%                       modulation function. The available options are
%                       listed below.
%                       If only a single keyword is given, the FM function
%                       is set to 'none'. For a pulse with constant
%                       amplitude, the AM needs to be specified as
%                       'rectangular'.
%                       Different AM functions can be multiplied by
%                       concatenating the keywords as 'AM1*AM2/FM'.
%                       (default: 'rectangular')
%     Par.*           = value for the pulse parameters defining the
%                       specified pulse shape. The pulse parameters
%                       required for each of the available modulation
%                       functions are listed below.
%     Par.I, Par.Q    = I and Q data describing an arbitrary pulse.
%                       The time axis is reconstructed based on Par.tp
%                       and the length of the I and Q vectors, all other
%                       input parameters (Amplitude, Flip, Frequency,
%                       Phase, etc.) are ignored.
%   If compensation for the resonator bandwidth or simulation of the effect
%   of the resonator is requested with Opt.Resonator, the following 
%   parameters also need to be defined:
%     Par.MagnitudeResponse  = magnitude response function (ideal or 
%                              experimental)
%     Par.faxis              = corresponding frequency axis in GHz
%     Par.mwFreq             = microwave frequency for the experiment in GHz
%   or
%     Par.ResonatorFrequency = resonator center frequency in GHz
%     Par.QL                 = loaded resonator Q-value
%     Par.mwFreq             = microwave frequency for the experiment in GHz
%   If compensation for transmitter nonlinearity or simulation of amplitude
%   compression is requested with Opt.Transmitter, the following parameters
%   also need to be defined:
%     Par.OutputAmplitude   = output amplitude as a function of input 
%                             amplitude
%     Par.InputAmplitude    = range of input amplitudes (normalized to 1)
%
%   Opt = optional structure with the following fields
%    Opt.IQ               = on/off; complex-valued pulse (on) or
%                           real-valued pulse (off) (default: on).
%    Opt.Resonator        = 'off' (default)
%                           'simulation'     - simulates the effect of the 
%                                              resonator on the pulse shape
%                           'BWcompensation' - resonator bandwidth compensation
%                                              for frequency-swept pulses (see 
%                                              Ref. 2)  (available for linear 
%                                              chirp and sech pulses).
%                           This option requires definition of the
%                           parameters Par.faxis, Par.MagnitudeResponse and
%                           Par.mwFreq.
%    Opt.Transmitter      = 'off' (default)
%                           'simulation'   - simulates the effect of the 
%                                            transmitter nonlinearity
%                           'compensation' - compensate for the transmitter 
%                                            nonlinearity
%                           This option requires definition of the
%                           parameters Par.InputAmplitude, Par.OutputAmplitude
%                           and Par.Amplitude (relative amplitude on the
%                           Par.InputAmplitude scale).
%    Opt.OverSampleFactor = oversampling factor for the determination of the 
%                           time step (default: 10)
%    Opt.CutoffFactor     = cutoff factor for truncation of the impulse 
%                           response function (used for resonator simulation) 
%                           (default:1/1000)
%
% Available pulse modulation functions:
%   - Amplitude modulation: rectangular, gaussian, sinc, halfsin, quartersin,
%                           sech, WURST, Gaussian pulse cascades (G3, Q3, 
%                           custom coefficients using 'GaussianCascade', see
%                           private\GaussianCascadeCoefficients.txt for 
%                           details), Fourier-series pulses (I-BURP 1/2,
%                           SNOB i2/i3, custom coefficients using 
%                           'FourierSeries',see 
%                           private\FourierSeriesCoefficients.txt for details)
%   - Frequency modulation: none, linear, tanh, uniformQ
%
% The parameters required for the different modulation functions are:
% Amplitude modulation:
% 'rectangular'         - none
% 'gaussian'            - tFWHM     = FWHM in µs
%                       Alternatively:
%                       - trunc     = truncation parameter (0 to 1)
% 'sinc'                - zerocross = width between the first zero-
%                                     crossing points in µs
% 'sech'                - beta      = dimensionless truncation parameter
%                       - n         = order of the secant function argument
%                                     (default: 1);
%                                     asymmetric sech pulses can be
%                                     obtained by specifying two values,
%                                     e.g. [6 1]
% 'WURST'               - nwurst    = WURST n parameter (determining the
%                                     steepness of the amplitude function)
% 'halfsin'             - none
% 'quartersin'          - trise     = rise time in µs for quarter sine
%                                     weighting at the pulse edges
% 'GaussianCascade'     - A0        = list of relative amplitudes
%                       - x0        = list of positions (in fractions of
%                                     tp)
%                       - FWHM     = list of FWHM (in fractions of tp)
% 'FourierSeries'       - A0        = initial amplitude coefficient
%                       - An        = list of Fourier coefficients for cos
%                       - Bn        = list of Fourier coefficients for sin
%
% Frequency modulation:
% For all frequency-modulated pulses, the frequency sweep range needs to be
% defined in Par.Frequency (e.g. Par.Frequency = [-50 50])
% 'linear'              no additional parameters
% 'tanh'                - beta      = dimensionless truncation parameter
% 'uniformQ'            The frequency modulation is calculated as the
%                       integral of the squared amplitude modulation
%                       function (for nth order sech/tanh pulses or in
%                       general to obtain offset-independent adiabaticity
%                       pulses, see Ref. 3.
%
% Output:   t          = time axis for defined waveform in µs
%           IQ         = real and imaginary part of the pulse function
%           modulation = structure with amplitude (modulation.A, in MHz),
%                        frequency (modulation.freq, in MHz) and phase
%                        (modulation.phase, in rad) modulation functions
%
% References:
% 1. The conversion from flip angles or critical adiabaticities to amplitudes 
%    is performed using the approximations described in:
%    Jeschke, G., Pribitzer, S., Doll, A. Coherence Transfer by Passage
%    Pulses in Electron Paramagnetic Resonance Spectroscopy.
%    J. Phys. Chem. B 119, 13570–13582 (2015). (DOI: 10.1021/acs.jpcb.5b02964)
% 2. Chirps with variable rate to compensate for the resonator bandwidth.
%    The bandwidth compensation is implemented as described in:
%    Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., Adiabatic and
%    fast passage ultra-wideband inversion in pulsed EPR
%    J. Magn. Reson. 230, 27-39 (2013). (DOI: 10.1016/j.jmr.2013.01.002)
%    and
%    Pribitzer, S., Doll, A. & Jeschke, G. SPIDYAN, a MATLAB library for
%    simulating pulse EPR experiments with arbitrary waveform excitation.
%    J. Magn. Reson. 263, 45–54 (2016). (DOI: 10.1016/j.jmr.2015.12.014)
% 3. Shaped pulses with offset-independent adiabaticity and determination
%    of their frequency modulation functions is described in:
%    Garwood, M., DelaBarre, L., The return of the frequency sweep: 
%    designing adiabatic pulses for contemporary NMR.
%    J. Magn. Reson. 153, 155-177 (2001). (DOI: 10.1006/jmre.2001.2340)
%

function varargout = pulse(varargin)

% ----------------------------------------------------------------------- %
% Input argument parsing
% ----------------------------------------------------------------------- %
switch nargin
  case 0
    help(mfilename);
    return
  case 1 % [t,IQ] = pulse(Par)
    Par = varargin{1};
    Opt = struct;
  case 2 % [t,IQ] = pulse(Par,Opt)
    Par = varargin{1};
    Opt = varargin{2};
  otherwise
    error('pulse() needs 1 or 2 input arguments.')
end

if ~isstruct(Par)
  error('First input argument must be a structure.');
end
if ~isstruct(Opt)
  error('Second input argument must be a structure.');
end

plotResults = (nargout==0);

if plotResults
  calculateExciteProfile = true;
else
  calculateExciteProfile = false;
end

% Set parameters to defaults
% ----------------------------------------------------------------------- %
if ~isfield(Par,'tp')
  error('Pulse length not defined in Par.tp.')
end
if ~isfield(Par,'Flip') || isempty(Par.Flip)
  if ~isfield(Par,'Amplitude') && ~isfield(Par,'Qcrit')
    Par.Flip = pi;
    Par.Amplitude = [];
  else
    Par.Flip = [];
  end
else
  Par.Amplitude = [];
end

if ~isfield(Par,'Frequency') || isempty(Par.Frequency)
  Par.Frequency = 0; % MHz
end
if ~isfield(Par,'Phase')
  Par.Phase = 0; % rad
end

% Options
% ----------------------------------------------------------------------- %
if ~isfield(Opt,'IQ')
  Opt.IQ = true;
end
if ~isfield(Opt,'Transmitter')
  Opt.Transmitter = 'off';
end
if ~isfield(Opt,'Resonator')
  Opt.Resonator = 'off';
end
if ~isfield(Opt,'OverSampleFactor')
  Opt.OverSampleFactor = 10;
end
if ~isfield(Opt,'CutoffFactor')
  Opt.CutoffFactor = 1/1000;
end

% ----------------------------------------------------------------------- %
% Calculate pulse function
% ----------------------------------------------------------------------- %
modulation = struct;

% Check if pulse I and Q data is given
if (isfield(Par,'I') && ~isempty(Par.I)) || ...
   (isfield(Par,'Q') && ~isempty(Par.Q))
  
  if ~isfield(Par,'Type') || isempty(Par.Type)
    Par.Type = 'user-IQ';
  end
  
  if ~isfield(Par,'I')
    Par.I = zeros(size(Par.Q));
  elseif ~isfield(Par,'Q')
    Par.Q = zeros(size(Par.I));
  elseif size(Par.I)~=size(Par.Q)
    error('I and Q input data have different lengths.');
  end
  if ~isvector(Par.I) || ~isvector(Par.Q)
    error('I and Q input data should be vectors.');
  end
  
  t = linspace(0,Par.tp,numel(Par.I));
  IQ = complex(Par.I,Par.Q);
  
  if isfield(Par,'TimeStep') && ~isempty(Par.TimeStep)
    t0 = t;
    t = t0(1):Par.TimeStep:t0(end);
    IQ = interp1(t0,IQ,t,'spline');
  else
    Par.TimeStep = t(2)-t(1);
  end
  
  % ------------------------------------------------------------------- %
  % Calculate amplitude compression simulation or compensation
  % ------------------------------------------------------------------- %
  if strcmp(Opt.Transmitter,'simulation') || ...
      strcmp(Opt.Transmitter,'compensation')
    IQ = AmplitudeNonlinearity(IQ,Par,Opt.Transmitter);
  end
  
  modulation.A = [];
  modulation.freq = [];
  modulation.phase = [];
   
else
  
  % Set pulse shape to rectangular if it is not specified
  if ~isfield(Par,'Type') || isempty(Par.Type)
    Par.Type = 'rectangular';
  end
  
  % Determine pulse shape from input string
  shape = regexp(Par.Type,'/','split');
  switch numel(shape)
    case 1
      AmplitudeModulation = lower(shape);
      FrequencyModulation = 'none';
    case 2
      AmplitudeModulation = lower(regexp(shape{1},'*','split'));
      if ~isempty(shape{2})
        FrequencyModulation = lower(shape{2});
      else
        FrequencyModulation = 'none';
      end
    otherwise
      error('Par.Type must contain the pulse shape in the form ''AM/FM''.');
  end
  
  % Check that all the required pulse parameters are given
  for na = 1:numel(AmplitudeModulation)
    switch AmplitudeModulation{na}
      case 'rectangular'
        % no additional parameters needed
        
      case 'gaussian'
        if (~isfield(Par,'tFWHM') || isempty(Par.tFWHM)) && ...
            (~isfield(Par,'trunc') || isempty(Par.trunc))
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.tFWHM or Par.trunc for the Gaussian envelope.']);
        elseif ~isfield(Par,'tFWHM') && isfield(Par,'trunc') && ...
            ~isempty(Par.trunc)
          % Convert truncation parameter to tFWHM
          Par.tFWHM = sqrt(-(Par.tp^2)/log2(Par.trunc));
          if Par.tFWHM==0
            Par.tFWHM = Par.TimeStep/2;
          end
        end
        
      case 'sinc'
        
        if ~isfield(Par,'zerocross') || isempty(Par.zerocross)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.zerocross for the sinc envelope.']);
        end
        
      case 'halfsin'
        % no additional parameters needed

      case 'quartersin'
        
        if ~isfield(Par,'trise') || isempty(Par.trise)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.trise for the quartersine envelope.']);
        end
        
      case 'sech'
        
        if ~isfield(Par,'beta') || isempty(Par.beta)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.beta parameter (in 1/µs) for the sech envelope.']);
        end
        if ~isfield(Par,'n') || isempty(Par.n)
          Par.n = 1;
        end
        
      case 'wurst'
        
        if ~isfield(Par,'nwurst') || isempty(Par.nwurst)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.nwurst parameter for the WURST envelope.']);
        end
        if numel(Par.nwurst)~=1 || mod(Par.nwurst,1) || Par.nwurst<1
          error('Pulseshape.nwurst must be a nonnegative integer (1,2,...).');
        end
        
      case {'gaussiancascade','g3','q3'}
        
        if ~(isfield(Par,'A0') && isfield(Par,'x0') && isfield(Par,'FWHM'))
          
          if strcmp(AmplitudeModulation{na},'gaussiancascade')
            error(['The amplitudes A0, positions x0 and FWHM of the ',...
                   'Gaussians are required as input.'])
          else
            % Load parameters from file
            fname = 'private\GaussianCascadeCoefficients.txt';
            fid = fopen(fname);
            while 1
              s = fgetl(fid);
              id = strfind(lower(s),AmplitudeModulation{na});
              if id
                for k = 1:6; fgetl(fid); end
                s = fgetl(fid);
                Par.x0 = sscanf(s,'%f');
                s = fgetl(fid);
                Par.A0 = sscanf(s,'%f');
                s = fgetl(fid);
                Par.FWHM = sscanf(s,'%f');
              end
              term = strfind(s,'end');
              if term
                break
              end
            end
            fclose(fid);
            clear fid s fname id
          end
          
        elseif ~(numel(Par.A0)==numel(Par.x0) && ...
                 numel(Par.A0)==numel(Par.FWHM))
          error(['The same number of parameters is required for the A0, ',...
                 'x0 and FWHM inputs.'])
        end
        
      case {'fourierseries','i-burp 1','i-burp 2','snob i2','snob i3'}
        
        if (~(isfield(Par,'A0') && isfield(Par,'An') && isfield(Par,'Bn')) || ...
            (isempty(Par.A0) && isempty(Par.An) && isempty(Par.Bn)))
          
          if strcmp(AmplitudeModulation{na},'fourierseries')
            error('The Fourier coefficients A0, An and Bn are required as input.')
          else
            % Load Fourier coefficients from file
            fname = 'private\FourierSeriesCoefficients.txt';
            fid = fopen(fname);
            while 1
              s = fgetl(fid);
              id = strfind(lower(s),AmplitudeModulation{na});
              if id
                for k = 1:5; fgetl(fid); end
                s = fgetl(fid);
                Par.A0 = sscanf(s,'%f');
                s = fgetl(fid);
                Par.An = sscanf(s,'%f');
                s = fgetl(fid);
                Par.Bn = sscanf(s,'%f');
              end
              term = strfind(s,'end');
              if term
                break
              end
            end
            fclose(fid);
            clear fid s fname id
          end
          
        elseif numel(Par.An)~=numel(Par.Bn)
          error('The same number of A and B Fourier coefficients is required.')
        end
               
      otherwise
        
        error('The amplitude modulation function ''%s'' is not defined.',...
              AmplitudeModulation{na});
        
    end
  end
  
  switch FrequencyModulation
    
    case 'none'
      
      if numel(Par.Frequency)>1
        error(['Frequency modulation is set to ''none'', but a frequency ',...
               'range is given in Par.Frequency. Please define a single ',...
               'pulse frequency.']);
      end
      
    case 'linear'
      
      if numel(Par.Frequency)~=2
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify frequency range for the linear chirp in Par.Frequency ',...
          '(in MHz).']);
      end
      
    case 'tanh'
      
      if numel(Par.Frequency)~=2
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify frequency range for tanh in Par.Frequency (in MHz).']);
      end
      if ~isfield(Par,'beta') || isempty(Par.beta)
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify dimensionless Par.beta parameter for tanh.']);
      end
           
    case 'uniformq'
      
      if numel(Par.Frequency)~=2
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify frequency range in Par.Frequency (in MHz).']);
      end
      
    otherwise
      
      error('The frequency modulation function ''%s'' is not defined.',...
            FrequencyModulation);
      
  end
  
  if any(ismember(AmplitudeModulation,'sech')) && ...
      strcmp(FrequencyModulation,'tanh') && ...
      (isfield(Par,'n') && ~isempty(Par.n) && any(Par.n~=1))
    warning(['For uniform adiabaticity pulses with nth order sech amplitude',...
             'modulation use Par.Type = ''sech/uniformQ''.']);
  end
  % For frequency-swept pulses Qcrit takes precedence over Flip and Amplitude
  if isfield(Par,'Qcrit') && ~isempty(Par.Qcrit)
    if strcmp(FrequencyModulation,'none')
      warning(['For a pulse without frequency modulation, the Qcrit input ',...
               'parameter is ignored.']);
      Par.Qcrit = [];
      if ~isfield(Par,'Amplitude') || isempty(Par.Amplitude)
        Par.Flip = pi;
      end
    else
      Par.Amplitude = [];
      Par.Flip = [];
    end
  end
    
  if ~strcmp(Opt.Resonator,'off')
    
    % Bandwidth compensation is implemented for these pulses
    if strcmp(Opt.Resonator,'BWcompensation') && ... 
       (~(strcmp(FrequencyModulation,'linear') && ...
          (strcmp(AmplitudeModulation,'rectangular') || ...
           strcmp(AmplitudeModulation,'quartersin')))) && ...
       (~(strcmp(FrequencyModulation,'tanh') && ...
          strcmp(AmplitudeModulation,'sech'))) && ...
       (~(strcmp(FrequencyModulation,'uniformq') && ...
          strcmp(AmplitudeModulation,'sech')))
      error(['Resonator bandwidth compensation is not implemented for ',...
             'the selected pulse. See documentation for more details.']);
    end
    
    if isfield(Par,'faxis') && ~isempty(Par.faxis)
      if ~isfield(Par,'MagnitudeResponse') || isempty(Par.MagnitudeResponse)
        error(['Use of the resonator option with an experimental transfer ',...
               'function requires the inputs Par.faxis and Par.MagnitudeResponse.']);
      end
    elseif isfield(Par,'ResonatorFrequency') && ~isempty(Par.ResonatorFrequency)
      if ~isfield(Par,'QL') || isempty(Par.QL)
        error(['Use of the resonator option with an ideal transfer function ',...
               'requires the inputs Par.ResonatorFrequency and Par.QL.']);
      end
    else
      error(['Use of the resonator option requires definition of the ',...
             'resonator magnitude response function (Par.faxis and ',...
             'Par.MagnitudeResponse or Par.ResonatorFrequency and Par.QL).']);
    end
    
    if ~isfield(Par,'mwFreq') || isempty(Par.mwFreq)
      error('Par.mwFreq is required for the resonator option.');
    end
    
  end
  
  if ~strcmp(Opt.Transmitter,'off')
    
      if (~isfield(Par,'InputAmplitude') || isempty(Par.InputAmplitude)) || ...
          (~isfield(Par,'OutputAmplitude') || isempty(Par.OutputAmplitude))
        error(['Insufficient input for the simulation of amplitude compression or '...
          'amplitude nonlinearity compensation. ',...
          'Specify the transmitter nonlinearity (Par.InputAmplitude, ',...
          'Par.OutputAmplitude).']);
      end
      if ~isfield(Par,'Amplitude') || isempty(Par.Amplitude) || ...
          Par.Amplitude>max(Par.InputAmplitude)
        error(['Insufficient input for the simulation of amplitude compression '...
          'or amplitude nonlinearity compensation. The pulse amplitude relative ',...
          'to the maximum input amplitude is required in Par.Amplitude.']);
      end
    
  end
   
  % Estimate pulse bandwidth (for timestep determination)
  % --------------------------------------------------------------------- %
  if ~isfield(Par,'TimeStep')
    
    % Determine bandwidth of frequency modulation
    switch FrequencyModulation
      case 'none'
        FM_BW = 0;
      otherwise
        FM_BW = abs(Par.Frequency(2) - Par.Frequency(1));
    end
    
    % Determine bandwidth of amplitude modulation (from Fourier transform)
    dt = 0.1e-3; % test time step
    t0 = 0:dt:Par.tp;
    ti0 = t0-0.5*Par.tp;
    A0 = ones(1,numel(t0));
    for na = 1:numel(AmplitudeModulation) % preliminary calculation of AM function
      switch AmplitudeModulation{na}
        case 'rectangular'
          A0 = A0.*ones(1,numel(t0));
        case 'gaussian'
          A0 = A0.*exp(-(4*log(2)*ti0.^2)/Par.tFWHM^2);
        case 'sinc'
          x_ = (2*pi*ti0)/Par.zerocross;
          A1 = sin(x_)./x_;
          A1(isnan(A1)) = 1;
          A0 = A0.*(A1/max(A1));
        case 'halfsin'
          A1 = cos(pi*ti0/Par.tp);
          A0 = A0.*A1;
        case 'quartersin'
          % Pulse edges weighted with a quarter period of a sine wave
          A1 = ones(1,numel(t0));
          if Par.trise~=0 && 2*Par.trise<Par.tp
            tpartial = 0:dt:Par.trise;
            npts = numel(tpartial);
            A1(1:npts) = sin(tpartial*(pi/(2*Par.trise)));
            A1(end-npts+1:end) = A1(npts:-1:1);
          end
          A0 = A0.*A1;
        case 'sech'
          n = min(Par.n); % Par.n contains two fields for asymmetric pulses
          A0 = A0.*sech(Par.beta*2^(n-1)*(ti0/Par.tp).^n);
        case 'wurst'
          A0 = A0.*(1 - abs(sin(pi*ti0/Par.tp)).^Par.nwurst);
        case {'gaussiancascade','g3','q3'}
          A0 = zeros(1,numel(t0));
          for j = 1:numel(Par.A0)
            A0 = A0 + Par.A0(j)*exp(-(4*log(2)/(Par.FWHM(j)*Par.tp)^2)*(t0-Par.x0(j)*Par.tp).^2);
          end
          A0 = A0/max(A0);
        case {'fourierseries','i-burp 1','i-burp 2','snob i2','snob i3'}
          A0 = zeros(1,numel(t0)) + Par.A0;
          for j = 1:numel(Par.An)
            A0 = A0 + Par.An(j)*cos(j*2*pi*t0/Par.tp) + Par.Bn(j)*sin(j*2*pi*t0/Par.tp);
          end
          A0 = A0/max(A0);
      end
    end
    % Fourier transform
    if nextpow2(numel(t0))<10
      zf = 2^10;
    else
      zf = 4*2^nextpow2(numel(t0));
    end
    A0ft = abs(fftshift(fft(A0,zf)));
    f = fdaxis(dt,zf);
    intg = cumtrapz(A0ft);
    [dummy,indmax] = min(abs(intg-0.5*max(intg)));
    indbw = find(A0ft(indmax:end)>0.1*max(A0ft));
    AM_BW = 2*(f(indmax+indbw(end))-f(indmax));
    
    BW = max([FM_BW AM_BW]);
    
  end
  
  % Set up time axis
  if ~isfield(Par,'TimeStep') || isempty(Par.TimeStep)
    
    % Automatically determine appropriate time step
    % ------------------------------------------------------------------- %
    % Calculate maximum frequency offset
    maxFreq = max(abs(mean(Par.Frequency)+[-1 1]*BW/2));
    % Use Nyquist theorem to calculate time step, incl. oversampling
    if maxFreq~=0
      Nyquist_dt = 1/(2*maxFreq);
      Par.TimeStep = Nyquist_dt/Opt.OverSampleFactor;
    else
      Par.TimeStep = 0.002; % µs
    end
    if Par.TimeStep>Par.tp
      Par.TimeStep = Par.tp;
    end
    Par.TimeStep = Par.tp/round(Par.tp/Par.TimeStep); % last time point = tp
    
  end
  t = 0:Par.TimeStep:Par.tp;
  ti = t - Par.tp/2;
  nPoints = numel(t);
  
  % ------------------------------------------------------------------- %
  % Amplitude modulation function (modulation.A)
  % ------------------------------------------------------------------- %
  modulation.A = ones(1,nPoints);
  for na = 1:numel(AmplitudeModulation)
    switch AmplitudeModulation{na}
      
      case 'rectangular'
        
        A = ones(1,nPoints);
        
      case 'gaussian'
        
        A = exp(-(4*log(2)*ti.^2)/Par.tFWHM^2);
        
      case 'sinc'

        x_ = 2*pi*ti/Par.zerocross;
        A = sin(x_)./x_;
        A(isnan(A)) = 1;
        A = A/max(A);
        
      case 'halfsin'
        
        % Half sine
        A = cos(pi*ti/Par.tp);
        
      case 'quartersin'
        
        % Pulse edges weighted with a quarter period of a sine wave
        A = ones(1,nPoints);
        if Par.trise~=0 && 2*Par.trise<Par.tp
          tpartial = 0:Par.TimeStep:Par.trise;
          npts = numel(tpartial);
          A(1:npts) = sin(tpartial*(pi/(2*Par.trise)));
          A(end-npts+1:end) = A(npts:-1:1);
        end
        
      case 'sech'
        
        if numel(Par.n)==1 % symmetric pulse
          if Par.n==1 % reduces numerical errors
            A = sech((Par.beta/Par.tp)*ti);
          else
            A = sech(Par.beta*2^(Par.n-1)*(ti/Par.tp).^Par.n);
          end
        else % asymmetric pulse
          A(1:round(nPoints/2)) = sech(Par.beta*2^(Par.n(1)-1)*(ti(1:round(nPoints/2))/Par.tp).^Par.n(1));
          A(round(nPoints/2)+1:nPoints) = sech(Par.beta*2^(Par.n(2)-1)*(ti(round(nPoints/2)+1:nPoints)/Par.tp).^Par.n(2));
        end
        
      case 'wurst'
        
        A = 1 - abs(sin(pi*ti/Par.tp)).^Par.nwurst;
        
      case {'gaussiancascade','g3','q3'}
        
        A = zeros(1,numel(t));
        for j = 1:numel(Par.A0)
          A = A + Par.A0(j)*exp(-(4*log(2)/(Par.FWHM(j)*Par.tp)^2)*(t-Par.x0(j)*Par.tp).^2);
        end
        A = A/max(A);
        
      case {'fourierseries','i-burp 1','i-burp 2','snob i2','snob i3'}
        
        A = zeros(1,numel(t)) + Par.A0;
        for j = 1:numel(Par.An)
          A = A + Par.An(j)*cos(j*2*pi*t/Par.tp) + Par.Bn(j)*sin(j*2*pi*t/Par.tp);
        end
        A = A/max(A);
        
    end
    modulation.A = modulation.A.*A;
  end
  
  % ------------------------------------------------------------------- %
  % Frequency (modulation.freq) and phase (modulation.phase) modulation 
  % functions
  % ------------------------------------------------------------------- %
  switch FrequencyModulation
    
    case 'none'
      
      modulation.freq = zeros(1,nPoints);
      modulation.phase = zeros(1,nPoints);
      
    case 'linear'
      
      k = (Par.Frequency(2)-Par.Frequency(1))/Par.tp; % frequency sweep (chirp) rate
      modulation.freq = k*ti;
      modulation.phase = 2*pi*((k/2)*ti.^2);
      
    case 'tanh'
      
      % Determine BWinf parameter from BW and beta parameters
      % (the frequency is swept from -BW/2 to +BW/2)
      Par.BWinf = (Par.Frequency(2)-Par.Frequency(1))/tanh(Par.beta/2);
      
      modulation.freq = (Par.BWinf/2)*tanh((Par.beta/Par.tp)*ti);
      modulation.phase = (Par.BWinf/2)*(Par.tp/Par.beta)*...
        log(cosh((Par.beta/Par.tp)*ti));
      modulation.phase = 2*pi*modulation.phase;
      
    case 'uniformq'
      % The frequency modulation is calculated as the integral of the
      % squared amplitude modulation function (for nth order sech/tanh or
      % in general to obtain offset-independent adiabaticity pulses, see
      %   Garwood, M., DelaBarre, L., J. Magn. Reson. 153, 155-177 (2001).
      %   http://dx.doi.org/10.1006/jmre.2001.2340
      
      modulation.freq = cumtrapz(ti,modulation.A.^2)/trapz(ti,modulation.A.^2); % F2
      modulation.freq = (Par.Frequency(2)-Par.Frequency(1))*(modulation.freq-1/2);
      modulation.phase = 2*pi*cumtrapz(ti,modulation.freq);
      modulation.phase = modulation.phase + abs(min(modulation.phase)); % zero phase offset at pulse center
      
  end
  
  % ------------------------------------------------------------------- %
  % Calculate bandwidth compensation
  % ------------------------------------------------------------------- %
  if strcmp(Opt.Resonator,'BWcompensation')
    
    % Variable-rate chirps with resonator bandwidth compensation, as
    % described in:
    %   Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
    %   Adiabatic and fast passage ultra-wideband inversion in
    %   pulsed EPR. J. Magn. Reson. 230, 27–39 (2013).
    %   http://dx.doi.org/10.1016/j.jmr.2013.01.002
    % and
    %   Doll, A., Frequency-swept microwave pulses for electron spin
    %   resonance, PhD Dissertation (2016), ETH Zürich, (for sech pulses
    %   see chapter 8, section 8.3.2, p. 133).
    % Implemented as in SPIDYAN, see:
    %   Pribitzer, S., Doll, A. & Jeschke, G. SPIDYAN, a MATLAB library
    %   for simulating pulse EPR experiments with arbitrary waveform
    %   excitation. J. Magn. Reson. 263, 45–54 (2016).
    %   http://dx.doi.org/10.1016/j.jmr.2015.12.014
    
    % Original amplitude and frequency modulation functions
    nu0 = modulation.freq;
    A0 = modulation.A;
    
    % Resonator profile in the frequency range of the pulse
    newaxis = nu0 + mean(Par.Frequency) + Par.mwFreq*1e3; % MHz
    if isfield(Par,'faxis')
      f = Par.faxis*1e3; % GHz to MHz
      H = Par.MagnitudeResponse;
      if min(newaxis)<min(f) || max(newaxis)>max(f)
        error(['The frequency sweep width of the pulse extends further than ',...
          'the given resonator profile. Please provide the resonator ',...
          'profile for the complete pulse frequency sweep width.'])
      end
    elseif isfield(Par,'ResonatorFrequency')
      [f,H] = transferfunction('ideal',Par.ResonatorFrequency,Par.QL);
      f = f*1e3;
      H = abs(H);
    end
    
    profile = interp1(f,H,newaxis);
    if strcmp(AmplitudeModulation,'sech')
      profile = A0.*profile;
    end
    
    % Frequency dependence of t and time-to-frequency mapping
    c_ = trapz(nu0,1./profile.^2)/t(end); % c_ = 2*pi/Qref
                                          % Qref = reference adiabaticity
    t_f = cumtrapz(nu0,(1/c_)*profile.^-2);
    nu_adapted = interp1(t_f,nu0,t,'pchip');
    
    % New frequency, phase and amplitude modulation functions
    modulation.freq = nu_adapted;
    modulation.phase = 2*pi*cumtrapz(t,modulation.freq);
    modulation.phase = modulation.phase + abs(min(modulation.phase));  % zero phase offset at pulse center
    
    modulation.A = interp1(nu0,modulation.A,nu_adapted,'pchip');
    
  end
  
  % ------------------------------------------------------------------- %
  % Calculate amplitude compression simulation or compensation
  % ------------------------------------------------------------------- %
  if ~strcmp(Opt.Transmitter,'off')
    modulation.A = modulation.A/max(modulation.A);
    [modulation.A,Par] = AmplitudeNonlinearity(modulation.A,Par,Opt.Transmitter);
  end
  
  % ------------------------------------------------------------------- %
  % Determine pulse amplitude from flip angle (if only Par.Flip is given)
  % or from the critical adiabaticity
  % ------------------------------------------------------------------- %
  if (~isfield(Par,'Amplitude') || isempty(Par.Amplitude))
    
    if strcmp(FrequencyModulation,'none')
      % Amplitude-modulated pulses: flip angle = integral
      
      if (isfield(Par,'Flip') && ~isempty(Par.Flip))
        Par.Amplitude = Par.Flip/(2*pi*trapz(t,modulation.A));
      end
      
    else
      % Frequency-modulated pulses
      
      % Q_crit = (2*pi*v1max)^2/k = minimum adiabaticity on resonance
      %   see Jeschke et al. (2015) J. Phys. Chem. B, 119, 13570–13582.
      %   http://dx.doi.org/10.1021/acs.jpcb.5b02964
      
      if ((~isfield(Par,'Qcrit') || isempty(Par.Qcrit)) && ...
         (isfield(Par,'Flip') && ~isempty(Par.Flip)))
        if Par.Flip>pi
          error('Pulse amplitude calculation from flip angle not applicable for angles larger than pi.');
        end
        Par.Qcrit = (2/pi)*log(2/(1+cos(Par.Flip)));
        Par.Qcrit = min(Par.Qcrit,5); % set Q_crit to finite value if it is infinite or large
      end
      
      if strcmp(Opt.Resonator,'off')
        switch FrequencyModulation
          case 'linear'
            sweeprate = abs(Par.Frequency(2)-Par.Frequency(1))/Par.tp;
          case 'tanh'
            sweeprate = Par.beta*abs(Par.BWinf)/(2*Par.tp);
          case 'uniformq'
            % Q = w1max^2*A(t)^2/(BW*dnu/dt) see eq. 17 in
            %   Garwood, M., DelaBarre, L., J. Magn. Reson. 153, 155-177
            %   (2001).
            %   http://dx.doi.org/10.1006/jmre.2001.2340
            [dummy,ind] = min(abs(ti));
            dnu = abs(diff(2*pi*modulation.freq/(t(2)-t(1))));
            sweeprate = dnu(ind)/(2*pi*(modulation.A(ind))^2);
        end
      else
        % Numerical computation
        % Q = w1max^2*A(t)^2/(BW*dnu/dt)
        [dummy,ind] = min(abs(ti));
        dnu = abs(diff(2*pi*modulation.freq/(t(2)-t(1))));
        sweeprate = dnu(ind)/(2*pi*(modulation.A(ind))^2);
      end
      
      Par.Amplitude = sqrt(2*pi*Par.Qcrit*sweeprate)/(2*pi);
      
    end
  end
    
  % ------------------------------------------------------------------- %
  % Calculate pulse IQ function
  % ------------------------------------------------------------------- %
  modulation.A = Par.Amplitude*modulation.A;
  totalphase = modulation.phase + 2*pi*mean(Par.Frequency)*t + Par.Phase;
  IQ = modulation.A.*exp(1i*totalphase);
  
  % Real-valued pulse
  if ~Opt.IQ
    IQ = real(IQ);
  end
  
end

% ------------------------------------------------------------------- %
% Calculate effect of the resonator on the pulse shape
% ------------------------------------------------------------------- %
if strcmp(Opt.Resonator,'simulation')
  
  % Calculate transfer function and impulse response
  if isfield(Par,'faxis')
    [f,H] = transferfunction('experimental',Par.faxis,Par.MagnitudeResponse);
  elseif isfield(Par,'ResonatorFrequency')
    [f,H] = transferfunction('ideal',Par.ResonatorFrequency,Par.QL);
  end
  [t,h] = impulseresponse(f,H,Par.mwFreq,Opt);
  
  % Convert impulse response to pulse time axis
  t_ = 0:Par.TimeStep:t(end);
  h_ = interp1(t,h,t_,'spline');
  
  % Convolution of impulse response function with the pulse shape
  IQ = conv(IQ,h_)/sum(abs(h_));
  
  % Redefine time axis
  tp = (numel(IQ)-1)*Par.TimeStep;
  t = 0:Par.TimeStep:tp(end);
  
end

% --------------------------------------------------------------------- %
% Excitation profile calculation
% --------------------------------------------------------------------- %

if calculateExciteProfile
  
  [exprof.offsets,exprof.M] = exciteprofile(t,IQ);
  
end

% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
if plotResults
  
  clf
  S.f = gcf;
  set(S.f,'WindowStyle','normal','Name','pulse output',...
    'numbertitle','off','Units','Normalized',...
    'Position',[0.1,0.3,0.8,0.40],'Color',[1 1 1]*0.8,...
    'Toolbar','figure');
  
  % Set colors
  colI = [0 0 1];
  colQ = [1 0 0];
  colBW = [1 1 1]*0.8;
  colnu = [0 0 1];
  colx = [0 0.5 0];
  coly = [1 0 0];
  colz = [0 0 1];
  
  % Set positions
  width = 0.25;
  height = 0.55;
  sep = (1-3*width)/4;
  btm = 0.20;
  boxpos = 0.77;
  
  S.htext = uicontrol('Style','edit','String',['Type = ' Par.Type],...
    'FontSize',10,'FontWeight','bold',...
    'Enable','Inactive','Units','Normalized',...
    'Position',[0.25*sep,0.9,1-0.5*sep,0.075]);
  
  % IQ plot
  S.label(1) = uicontrol('Style','text','String','Pulse amplitude:',...
    'FontSize',10,'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Background',[1 1 1]*0.8,'Units','Normalized',...
    'Position',[sep,0.75,width,0.1]);
  S.tick(1) = uicontrol('Style','checkbox',...
    'String','I','Value',1,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[3*sep boxpos 0.1 0.1]);
  S.tick(2) = uicontrol('Style','checkbox',...
    'String','Q','Value',1,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[3.75*sep boxpos 0.1 0.1]);
  S.ha(1) = axes('Units','Normalized','Position',[sep,btm,width,height]);
  hold on; box on;
  line([min(t) max(t)],[0 0],'Color',colBW);
  if ~isempty(modulation.A)
    if strcmp(FrequencyModulation,'none') && Par.Frequency==0
      S.hA = plot(t,modulation.A);
    else
      S.hA = plot(t,modulation.A,t,-modulation.A);
    end
    set(S.hA,'Color',[1 1 1]*0.9);
  end
  S.hI = plot(t,real(IQ),'Color',colI);
  S.hQ = plot(t,imag(IQ),'Color',colQ);
  Amax = max(abs(IQ));
  axis([t(1) t(end) -1*Amax*1.1 1*Amax*1.1]);
  xlabel(['{\itt} (',char(181),'s)'])
  ylabel('amplitude (MHz)')
  set(gca,'Layer','top')
  legend([S.hI S.hQ],'I','Q','Location','SouthEast')
  if ~Opt.IQ
    set(S.hQ,'Visible','off')
    legend(S.hI,'I','Location','SouthEast')
    set(S.tick(2),'Value',0,'Enable','off')
  end
  
  % Frequency modulation plot
  S.label(2) = uicontrol('Style','text','String','Frequency modulation:',...
    'FontSize',10,'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Background',[1 1 1]*0.8,'Units','Normalized',...
    'Position',[2*sep+width,0.75,width,0.1]);
  S.ha(2) = axes('Units','Normalized','Position',[2*sep+width,btm,width,height]);
  hold on; box on;
  if ~isempty(modulation.freq)
    line([min(t) max(t)],[0 0]+mean(Par.Frequency),'Color',colBW);
    if numel(Par.Frequency)==2
      line([min(t) max(t)],[0 0]+Par.Frequency(1),'Color',colBW);
      line([min(t) max(t)],[0 0]+Par.Frequency(2),'Color',colBW);
    end
    S.hnu = plot(t,modulation.freq+mean(Par.Frequency),'Color',colnu);
    freqmax = [min(modulation.freq) max(modulation.freq)]+mean(Par.Frequency);
    if freqmax(1)==freqmax(2)
      if freqmax(1)==0; sc = 1; else; sc = 0.1*freqmax(1); end
      freqmax = [freqmax(1) freqmax(2)]+sc*[-1 1];
      shift = 0;
    else
      shift = [-1 1]*0.1*(freqmax(2)-freqmax(1));
    end
    xlabel(['{\itt} (',char(181),'s)']);
    ylabel('frequency (MHz)');
    axis([t(1) t(end) freqmax+shift]);
  end
  
  % Excitation profile plot
  S.label(3) = uicontrol('Style','text','String','Excitation profiles:',...
    'FontSize',10,'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Background',[1 1 1]*0.8,'Units','Normalized',...
    'Position',[3.5*sep+2*width,0.75,width,0.1]);
  S.tick(3) = uicontrol('Style','checkbox',...
    'String','x','Value',0,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[3.5*sep+2.5*width boxpos 0.1 0.1]);
  S.tick(4) = uicontrol('Style','checkbox',...
    'String','y','Value',0,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[4.25*sep+2.5*width boxpos 0.1 0.1]);
  S.tick(5) = uicontrol('Style','checkbox',...
    'String','z','Value',1,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[5*sep+2.5*width boxpos 0.1 0.1]);
  S.ha(3) = axes('Units','Normalized','Position',[3.5*sep+2*width,btm,width,height]);
  hold on; box on;
  line([min(exprof.offsets) max(exprof.offsets)],[0 0],'Color',colBW);
  line([1 1]*mean(Par.Frequency),[-1 1],'Color',colBW);
  if numel(Par.Frequency)==2
    line([1 1]*Par.Frequency(1),[-1 1],'Color',colBW);
    line([1 1]*Par.Frequency(2),[-1 1],'Color',colBW);
  end
  if ~Opt.IQ
    line([0 0],[-1 1],'Color','k');
    line(-[1 1]*mean(Par.Frequency),[-1 1],'Color',colBW);
    if numel(Par.Frequency)==2
      line(-[1 1]*Par.Frequency(1),[-1 1],'Color',colBW);
      line(-[1 1]*Par.Frequency(2),[-1 1],'Color',colBW);
    end
  end
  S.h(1) = plot(exprof.offsets,exprof.M(1,:),'Color',colx);
  S.h(2) = plot(exprof.offsets,exprof.M(2,:),'Color',coly);
  S.h(3) = plot(exprof.offsets,exprof.M(3,:),'Color',colz);
  set(S.h(1),'Visible','off')
  set(S.h(2),'Visible','off')
  ylabel('{\itM}_i/{\itM}_0')
  legend(S.h,'x','y','z','Location','SouthEast')
  xlabel('frequency (MHz)')
  axis([exprof.offsets(1) exprof.offsets(end) -1 1])
  set(gca,'Layer','top')
  
  S.handles = [S.hI S.hQ S.h(1) S.h(2) S.h(3)];
  set(S.tick,'Callback',{@showhide,S});
  
end

% ----------------------------------------------------------------------- %
% Output
% ----------------------------------------------------------------------- %
switch nargout
  case 0
    % plotting
  case 2 % [t,IQ] = pulse(...)
    varargout = {t,IQ};
  case 3 % [t,IQ,modulation] = pulse(...)
    modulation.freq = modulation.freq + mean(Par.Frequency);
    modulation.phase = modulation.phase + 2*pi*mean(Par.Frequency)*t + Par.Phase;
    varargout = {t,IQ,modulation};
  otherwise
    error('The function pulse() needs 2 or 3 output arguments.')
end

end


% ------------------------------------------------------------------- %
% ------------------------------------------------------------------- %
% Subfunctions
% ------------------------------------------------------------------- %
% ------------------------------------------------------------------- %

function [f,H] = transferfunction(type,varargin)
% Calculate resonator transfer function for a given center frequency
% and loaded Q-value (type = 'ideal') or by extrapolation from an
% experimentally measured transfer function in a limited range
%
% See:
% 1. Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
%    Adiabatic and fast passage ultra-wideband inversion in
%    pulsed EPR. J. Magn. Reson. 230, 27–39 (2013).
%    http://dx.doi.org/10.1016/j.jmr.2013.01.002
% 2. Pribitzer, S., Doll, A. & Jeschke, G. SPIDYAN, a MATLAB library
%    for simulating pulse EPR experiments with arbitrary waveform
%    excitation. J. Magn. Reson. 263, 45–54 (2016). 
%    http://dx.doi.org/10.1016/j.jmr.2015.12.014

% Ideal transfer function (RLC series circuit)
Hideal = @(f,f0,Q,nu_max) nu_max./(1+1i*Q*(f/f0-f0./f));

switch type
  case 'ideal'
    
    % Calculate ideal resonator transfer function
    f0 = varargin{1}; % center frequency
    QL = varargin{2}; % loaded Q-value
    
    % Frequency axis
    fmax = 2*f0;
    N = 2^15; % number of frequency points
    df = 2*fmax/N; % GHz
    f = (-fmax:df:fmax-df); % GHz
    
    H = Hideal(f,f0,QL,1);
    
  case 'experimental'
    
    % Extrapolation of the experimental transfer function over
    % the full frequency range
    faxis = varargin{1}; % frequency axis in GHz
    MagnitudeResponse = varargin{2}; % experimental resonator profile
    
    % Frequency axis
    fmax = 2*ceil(max(faxis)); % Nyquist frequency, GHz
    N = 2^15; % number of frequency points
    df = 2*fmax/N; % GHz
    f = (-fmax:df:fmax-df); % GHz
    
    % Interpolation
    ind = find((f>faxis(1)) & (f<faxis(end)));
    MagnitudeResponse_ = interp1(faxis,MagnitudeResponse,f(ind),'spline');
    
    % Estimate center frequency and loaded Q
    [v1max,maxind] = max(MagnitudeResponse_);
    f0 = f(ind(1)+maxind); % GHz
    
    v1_3dB = (3/4)*v1max;
    ind_3dB = [find(MagnitudeResponse_>v1_3dB,1,'first') find(MagnitudeResponse_>v1_3dB,1,'last')];
    f_3dB = f(ind(1)-1+ind_3dB);
    BW_3dB = diff(f_3dB);
    QL = f0/BW_3dB;
    
    % Fit f0 and QL for best overlap
    fitfunc = @(x) sqrt(sum((abs(Hideal(f(ind),x(1),x(2),v1max))-MagnitudeResponse_).^2)/numel(MagnitudeResponse_));
    x0 = [f0 QL];
    x = fminsearch(fitfunc,x0);
    f0 = x(1);
    QL = x(2);
    
    Hid = Hideal(f,f0,QL,v1max);
    
    % Extrapolate the experimental transfer function by exponentially approaching |H0(f)|
    H0(1:ind(1)) = (MagnitudeResponse_(1)-abs(Hid(ind(1))))*exp(1/log(2)*(f(1:ind(1))-f(ind(1)))) + abs(Hid(1:ind(1)));
    H0(ind) = MagnitudeResponse_;
    H0(ind(end):numel(f)) = (MagnitudeResponse_(end)-abs(Hid(ind(end))))*exp(-1/log(2)*(f(ind(end):end)-f(ind(end)))) + abs(Hid(ind(end):end));
    
    % Estimate phase response
    betaid=atan(imag(Hid)./real(Hid));
    betaexp = -imag(hilberttrans(log(abs(H0))));
    % Combine phase response of the ideal function and the provided measured data
    fitfunc = @(x) sqrt(sum(((betaexp(ind)-x)-betaid(ind)).^2)/(ind(end)-ind(1)));
    x = fminsearch(fitfunc,0);
    beta(ind) = betaexp(ind)-x;
    beta(1:ind(1)) = ((betaexp(ind(1))-x)-betaid(ind(1)))*exp(1/log(2)*(f(1:ind(1))-f(ind(1)))) + betaid(1:ind(1));
    beta(ind(end):numel(f)) = ((betaexp(ind(end))-x)-betaid(ind(end)))*exp(-1/log(2)*(f(ind(end):end)-f(ind(end)))) + betaid(ind(end):end);
    
    % Transfer function
    H = abs(H0).*exp(1i*beta);
    ind0 = find(f==0);
    H(ind0) = H(ind0+1);
    
    % Define negative frequency as symmetric to the positive frequency
    H(1:numel(H)/2) = fliplr(H(numel(H)/2+1:end));
    
end

end

function [t,h] = impulseresponse(f,H,mwFreq,Opt)
% Calculate the impulse response function for the transfer function
% H (including downconversion by mwFreq)
%
% See:
% 1. Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
%    Adiabatic and fast passage ultra-wideband inversion in
%    pulsed EPR. J. Magn. Reson. 230, 27–39 (2013).
%    http://dx.doi.org/10.1016/j.jmr.2013.01.002
% 2. Pribitzer, S., Doll, A. & Jeschke, G. SPIDYAN, a MATLAB library
%    for simulating pulse EPR experiments with arbitrary waveform
%    excitation. J. Magn. Reson. 263, 45–54 (2016). 
%    http://dx.doi.org/10.1016/j.jmr.2015.12.014

% Normalize the transfer function and make it causal
H = H/max(H);
H_ = real(H) - 1i*imag(hilbert(real(H)));

% Impulse response function
df = unique(diff(f)); % GHz
dt = 1/(df(1)*numel(f)); % ns
t = -(numel(f)/2)*dt:dt:(numel(f)/2-1)*dt;
h0 = fftshift(ifft(ifftshift(H_)));

% Transformation to rotating frame
if mwFreq~=0
  [t,h_] = rfmixer(t*1e-3,real(h0),mwFreq,'IQdemod'); % time axis in µs
else
  h_ = h0;
end

% Windowing
h = h_.*chebwin(numel(h_)).';
h = h(round(numel(h_)/2)+1:end);
t = t(round(numel(h_)/2)+1:end); 

% Truncation
ind = find(fliplr(abs(h))>(Opt.CutoffFactor*max(abs(h))),1);
t = t(1:end-ind+1);
h = h(1:end-ind+1);

end

function [signal_out,varargout] = AmplitudeNonlinearity(A,Par,option)
% Compensation for/simulation of transmitter nonlinearity as
% described in:
%   Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
%   Adiabatic and fast passage ultra-wideband inversion in
%   pulsed EPR. J. Magn. Reson. 230, 27–39 (2013).
%   http://dx.doi.org/10.1016/j.jmr.2013.01.002

  % Normalize input amplitude
  maxInputAmplitude = max(Par.InputAmplitude);
  if maxInputAmplitude>1    
    Par.InputAmplitude = Par.InputAmplitude/maxInputAmplitude;
    Par.Amplitude = Par.Amplitude/maxInputAmplitude;
  end
  
  % Scale output amplitude
  nu1_max = max(Par.OutputAmplitude);
  Par.OutputAmplitude = Par.OutputAmplitude/nu1_max;
  
  switch option
    case 'simulation'

      F = griddedInterpolant(Par.InputAmplitude,Par.OutputAmplitude,'spline');
      
      % Compressed amplitude modulation function (renormalized to A0 = 1)
      if isreal(A)
        signal_out = sign(A).*F(abs(A)*Par.Amplitude);
      else
        A = A/max(abs(A));
        signal_out = sign(real(A)).*F(abs(real(A))*Par.Amplitude) + ...
                  1i*sign(imag(A)).*F(abs(imag(A))*Par.Amplitude);
      end
      
      if (~isfield(Par,'Flip') || isempty(Par.Flip)) && ...
          (~isfield(Par,'Qcrit') || isempty(Par.Qcrit))
        % Signal amplitude determined by provided transmitter 
        % nonlinearity characterization (for simulation under the
        % conditions of actual experiments)
        signal_out = nu1_max*signal_out;
        Par.Amplitude = 1;
      else
        % Re-normalize (if amplitude is defined through the flip
        % angle or Qcrit -> for simulation of the effects of
        % compression on the pulse shape)
        signal_out = signal_out/F(Par.Amplitude);
        Par.Amplitude = [];
      end

    case 'compensation'

      F = griddedInterpolant(Par.OutputAmplitude,Par.InputAmplitude,'spline');

      % Compensated amplitude function
      if isreal(A)
        signal_out = sign(A).*F(abs(A)*Par.Amplitude);
      else
        A = A/max(abs(A));
        signal_out = sign(real(A)).*F(abs(real(A))*Par.Amplitude) + ...
                  1i*sign(imag(A)).*F(abs(imag(A))*Par.Amplitude);
      end
      
      Par.Amplitude = 1;

  end
  
  if nargout==2
    varargout{1} = Par;
  end
  
end

function showhide(varargin)
% Callback for tick boxes

S = varargin{3}; % get calling handle structure

for i = 1:numel(S.tick)
  val = get(S.tick(i),'Value');
  if val==1
    set(S.handles(i),'Visible','on')
  else
    set(S.handles(i),'Visible','off');
  end
end

end