function varargout = pulse(varargin)

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% Pulse definition and excitation profile calculation
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%
% [t,y] = pulse(Exp)
% [t,y] = pulse(Exp,npulse)
% [t,y,mod] = pulse(Exp,npulse)
% [t,y] = pulse(Exp,npulse,Opt)
%
% [t,y,p] = pulse(Exp,npulse,Opt) if Opt.ExciteProfile = 1
% [t,y,mod,p] = pulse(Exp,npulse,Opt) if Opt.ExciteProfile = 1
%
% Input: Exp    = structure containing the following fields:
%          Exp.tp          = pulse length in us
%          Exp.timestep    = time step for waveform definition in us
%                            (default = determined based on input pulse
%                            parameters)
%          Exp.Amplitude   = pulse amplitude in MHz
%          Exp.Flip        = pulse flip angle in radians (*1)
%            (if Exp.Amplitude and Exp.Flip are both missing, a pulse with
%             amplitude 1 is returned; if both Exp.Amplitude and Exp.Flip
%             are given, Exp.Flip is ignored)
%          Exp.FreqOffset  = pulse center frequency offset with respect to
%                            Exp.mwFreq in MHz (default = 0)
%          Exp.Phase       = phase for the pulse in radians (default = 0,
%                            i.e. +x)
%          Exp.Pulse.Shape = pulse shape name in a string with structure
%                            'AM/FM' (or just 'AM'), where AM refers to the
%                            amplitude modulation function and FM to the
%                            frequency modulation function, the available
%                            options are listed below
%                            If only a single keyword is given, the FM
%                            function is set to 'none'. For a pulse with
%                            constant amplitude, the AM needs to be
%                            specified as 'rectangular'.
%                            Different AM functions can be multiplied by
%                            concatenating the keywords as 'AM1*AM2/FM'.
%                            (default = 'rectangular(/none)')
%          Exp.Pulse.par   = value for the pulse parameter 'par' defining
%                            the specified pulse shape, the pulse
%                            parameters required for each of the available
%                            modulation functions are listed below.
%          Exp.Pulse.I/Q   = I and Q data describing an arbitrary pulse,
%                            the time axis is reconstructed based on Exp.tp
%                            and the length of the I and Q vectors, all
%                            other input parameters (Amplitude, Flip,
%                            FreqOffset, Phase, etc.) are ignored
%
%        npulse = index/indices for the pulses in Exp input structure for
%                 which to return the pulse data
%
%        Opt    = optional structure with the following fields
%          Opt.OverSampleFactor = oversampling factor for the determination
%                                 of the time step (default = 10)
%          Opt.plot             = true/false, turn plotting of the results
%                                 on/off
%          Opt.ExciteProfile    = true/false, turn excitation profile
%                                 calculation on/off
%          Opt.offsets          = axis of frequency offsets in MHz for which
%                                 to compute the excitation profile
%                                 (default ±200 MHz, 201 pts or ±1.5*BW in
%                                 1 MHz steps if a bandwidth is defined for
%                                 one of the pulses)
%          Opt.nBCH             = number of steps combined in the excitation
%                                 profile computation using the Baker-
%                                 Campbell-Hausdorff series (*2)
%                                 (default = 3)
%
% Available pulse modulation functions:
%   - Amplitude modulation: rectangular, gaussian, sinc, quartersin, sech,
%                           nth order sech, WURST,
%   - Frequency modulation: linear, tanh, BW compensated, uniform
%                           adiabaticity
% ('nth order sech/tanh' is interpreted as 'nth order sech/uniform adiabaticity')
%
% The parameters required for the different modulation functions are:
% Amplitude modulation:
% 'rectangular'         - none
% 'gaussian'            - tFWHM     = FWHM in us
%                       Alternatively:
%                       - trunc     = truncation parameter (0 to 1)
% 'sinc'                - zerocross = width between the first zero-
%                                     crossing points in ns
% 'sech'                - beta      = dimensionless truncation parameter
% 'WURST'               - n         = WURST n parameter (determining the
%                                     steepness of the amplitude function)
% 'quartersin'          - trise     = rise time in us for quarter sine
%                                     weighting at the pulse edges
% 'nth order sech'      - beta      = dimensionless truncation parameter
%                       - n         = order of the secant function argument
%
% Frequency modulation:
% 'linear'              - BW        = frequency sweep width in MHz
%                                     (the frequency is swept over ±BW/2)
% 'tanh'                - BW        = frequency sweep width in MHz
%                       - beta      = dimensionless truncation parameter
%                       - SweepDirection = +1 or -1 (+1 default)
% 'BW compensated'(*3)  - BW        = frequency sweep width in MHz
%                                     (the frequency is swept over ±BW/2)
%                       Parameters required for resonator bandwidth
%                       compensation:
%                       - freqaxis  = frequency axis
%                       - v1        = magnitude response function (ideal or
%                                     experimental)
%                       - Exp.mwFreq needs to be defined
% 'uniform adiabaticity'- BW        = frequency sweep width in MHz
%                       the frequency modulation is calculated as the
%                       integral of the squared amplitude modulation
%                       function (for nth order sech/tanh pulses or in
%                       general to obtain offset-independent adiabaticity
%                       pulses, see Garwood, M., DelaBarre, L., J. Magn.
%                       Reson. 153, 155-177 (2001).
%
% Output:   t         = time axis for defined waveform in us
%           y         = real and imaginary part of the pulse function
%           mod       = structure with amplitude (mod.A in MHz), frequency
%                       (mod.nu in MHz) and phase (mod.phi in rad)
%                       modulation functions
%           Additionally, if Opt.ExciteProfile = 1:
%           p.dnu     = frequency offset axis for excitation profile in MHz
%           p.Mz      = excitation profile
%
% (*1) The conversion from flip angles to amplitudes is performed using the
%      approximations described in:
%      Jeschke, G., Pribitzer, S., Doll, A. Coherence Transfer by Passage
%      Pulses in Electron Paramagnetic Resonance Spectroscopy.
%      J. Phys. Chem. B 119, 13570–13582 (2015).
% (*2) Ernst, Bodenhausen, Wokaun, Principles of NMR in One and Two
%      Dimensions, Clarendon Press (1990), Chapter 3, p. 72-75.
% (*3) Chirps with variable rate to compensate for the resonator bandwidth.
%      The bandwidth compensation is implemented as described in:
%      Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., Adiabatic and
%      fast passage ultra-wideband inversion in pulsed EPR
%      J. Magn. Reson. 230, 27-39 (2013).
%

% ----------------------------------------------------------------------- %
% Input argument parsing
% ----------------------------------------------------------------------- %
Exp = varargin{1};
if nargin==2 % [t,y] = pulse(Exp,npulse) or [t,y] = pulse(Exp,Opt)
  if isstruct(varargin{2})
    Opt = varargin{2};
  else
    npulse = varargin{2};
  end
elseif nargin==3 % [t,y] = pulse(Exp,npulse,Opt)
  npulse = varargin{2};
  Opt = varargin{3};
elseif nargin>3
  error('The function pulse is not supported for more than 3 input arguments.')
end
if nargout==0
  Opt.plot = 1;
end

% Set experiment and option parameters to defaults
if ~isfield(Exp,'tp')
  error('Pulse length not defined in Exp.tp.')
end
if ~isfield(Exp,'Pulse')
  for i = 1:numel(Exp.tp)
    Exp.Pulse(i).Shape = 'rectangular';
  end
end
if ~isfield(Exp,'Amplitude') && ~isfield(Exp,'Flip')
  Exp.Amplitude(1:numel(Exp.tp)) = 1; % normalized amplitude
elseif isfield(Exp,'Amplitude') && numel(Exp.Amplitude)~=numel(Exp.tp)
  Exp.Amplitude(end:numel(Exp.tp)) = Exp.Amplitude(1);
elseif isfield(Exp,'Flip') && numel(Exp.Flip)~=numel(Exp.tp)
  Exp.Flip(end:numel(Exp.tp)) = Exp.Flip(end);
end
if ~isfield(Exp,'FreqOffset')
  Exp.FreqOffset(1:numel(Exp.tp)) = 0; % MHz
elseif numel(Exp.FreqOffset)~=numel(Exp.tp)
  Exp.FreqOffset(end:numel(Exp.tp)) = Exp.FreqOffset(end);
end
if ~isfield(Exp,'Phase')
  Exp.Phase(1:numel(Exp.tp)) = 0; % rad
elseif numel(Exp.Phase)~=numel(Exp.tp)
  Exp.Phase(end:numel(Exp.tp)) = Exp.Phase(end);
end
if isfield(Exp,'timestep') && numel(Exp.timestep)~=numel(Exp.tp)
  Exp.timestep(end:numel(npulse)) = Exp.timestep(end);
end
% If no pulse is selected, the functions returns outputs for each of the
% pulses defined in the input structure
if ~exist('npulse','var')
  npulse = 1:numel(Exp.tp);
end
if ~exist('Opt','var')
  Opt.ExciteProfile = 0;
  Opt.nBCH = 3;
  Opt.plot = 0;
end
if ~isfield(Opt,'plot')
  Opt.plot = 0;
end
if ~isfield(Opt,'OverSampleFactor')
  Opt.OverSampleFactor = 10;
end
if ~isfield(Opt,'ExciteProfile')
  Opt.ExciteProfile = 0;
end
if ~isfield(Opt,'nBCH')
  Opt.nBCH = 3;
end
if ~isfield(Opt,'offsets')
  if isfield(Exp.Pulse,'BW') % set range of offsets to ±1.5*BW in 1 MHz steps
    BW(1:numel(Exp.tp)) = 0;
    for i = 1:numel(Exp.tp)
      if ~isempty(Exp.Pulse(i).BW)
        BW(i) = Exp.Pulse(i).BW;
      end
    end
    Opt.offsets = -1.5*max(BW):1:1.5*max(BW);
  else
    limit = 200; % MHz, default
    npoints = 201;
    Opt.offsets = linspace(-limit,limit,npoints);
  end
end

% ----------------------------------------------------------------------- %
% Loop over pulses defined in input structure and calculate pulse function
% ----------------------------------------------------------------------- %
t = cell(1,numel(npulse));
y = cell(1,numel(npulse));
mod(1:numel(npulse)) = struct('A',[]);
if Opt.ExciteProfile==1
  p(npulse) = struct('offsets',Opt.offsets);
end
for np = 1:numel(npulse)
  
  n = npulse(np);
  
  % Check if pulse I and Q data is given
  if (isfield(Exp.Pulse(n),'I') && ~isempty(Exp.Pulse(n).I)) || ...
      (isfield(Exp.Pulse(n),'Q') && ~isempty(Exp.Pulse(n).Q))
    
    if ~isfield(Exp.Pulse(n),'Shape') || isempty(Exp.Pulse(n).Shape)
      Exp.Pulse(n).Shape = 'user-defined';
    end
    
    if ~isfield(Exp.Pulse(n),'I')
      Exp.Pulse(n).I = zeros(size(Exp.Pulse(n).Q));
    elseif ~isfield(Exp.Pulse(n),'Q')
      Exp.Pulse(n).Q = zeros(size(Exp.Pulse(n).I));
    elseif size(Exp.Pulse(n).I)~=size(Exp.Pulse(n).Q);
      error(['I and Q input data for pulse ',num2str(n),' have different lengths.'])
    end
    if ~isvector(Exp.Pulse(n).I) || ~isvector(Exp.Pulse(n).Q)
      error(['I and Q input data for pulse ',num2str(n),' should be vectors.'])
    end
    
    t{n} = linspace(0,Exp.tp(n),numel(Exp.Pulse(n).I));
    y{n} = Exp.Pulse(n).I + 1i*Exp.Pulse(n).Q;
    
    Exp.timestep(n) = t{n}(2)-t{n}(1);
    mod(n).A = [];
    mod(n).dnu = [];
    mod(n).phi = [];
    
  else
    
    % Set pulse shape to rectangular if it is not specified
    if ~isfield(Exp.Pulse(n),'Shape') || isempty(Exp.Pulse(n).Shape) || strcmp(Exp.Pulse(n).Shape,'')
      Exp.Pulse(n).Shape = 'rectangular';
    end
    
    % Determine pulse shape from input string
    shape = regexp(Exp.Pulse(n).Shape,'/','split');
    if numel(shape)==1
      AmplitudeModulation = shape;
      FrequencyModulation = 'none';
    else
      FrequencyModulation = shape{2};
      AmplitudeModulation = regexp(shape{1},'*','split');
    end
    
    % Check that all the required pulse parameters are given
    for na = 1:numel(AmplitudeModulation)
      switch AmplitudeModulation{na}
        
        case 'rectangular'
          
        case 'gaussian'
          
          if (~isfield(Exp.Pulse(n),'tFWHM') || isempty(Exp.Pulse(n).tFWHM)) && ...
              (~isfield(Exp.Pulse(n),'trunc') || isempty(Exp.Pulse(n).trunc))
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify tFWHM or trunc parameter for the Gaussian envelope.']);
          elseif ~isfield(Exp.Pulse(n),'tFWHM') && isfield(Exp.Pulse(n),'trunc') && ...
              ~isempty(Exp.Pulse(n).trunc)
            % Convert truncation parameter to tFWHM
            Exp.Pulse(n).tFWHM = sqrt(-(Exp.tp(n)^2)/log2(Exp.Pulse(n).trunc));
            if Exp.Pulse(n).tFWHM==0
              Exp.Pulse(n).tFWHM = Exp.timestep(n)/2;
            end
          end
          
        case 'sinc'
          
          if ~isfield(Exp.Pulse(n),'zerocross') || isempty(Exp.Pulse(n).zerocross)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify zero-crossing parameter in us for the sinc envelope.']);
          end
          
        case 'quartersin'
          
          if ~isfield(Exp.Pulse(n),'trise') || isempty(Exp.Pulse(n).trise)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify trise parameter in us for the quartersine envelope.']);
          end
          
        case 'sech'
          
          if ~isfield(Exp.Pulse(n),'beta') || isempty(Exp.Pulse(n).beta)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify beta parameter in 1/us for the sech envelope.']);
          end
          
        case 'nth order sech'
          
          if ~isfield(Exp.Pulse(n),'beta') || isempty(Exp.Pulse(n).beta)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify beta parameter in 1/us for the sech envelope.']);
          end
          if ~isfield(Exp.Pulse(n),'n') || isempty(Exp.Pulse(n).n)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify n for the nth order sech envelope.']);
          end
          
          if strcmp(FrequencyModulation,'tanh')
            FrequencyModulation = 'uniform adiabaticity';
          end
          
        case 'WURST'
          
          if ~isfield(Exp.Pulse(n),'n') || isempty(Exp.Pulse(n).n)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify n parameter for the WURST envelope.']);
          end
          
        otherwise
          
          error(['The amplitude modulation function selected for pulse ',num2str(n),' is not defined.']);
          
      end
    end
    
    switch FrequencyModulation
      
      case 'none'
        
      case 'linear'
        
        if ~isfield(Exp.Pulse(n),'BW') || isempty(Exp.Pulse(n).BW)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify BW parameter in MHz for the linear chirp.']);
        end
        
      case 'tanh'
        
        if ~isfield(Exp.Pulse(n),'BW') || isempty(Exp.Pulse(n).BW)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify BW parameter in MHz for tanh.']);
        end
        if ~isfield(Exp.Pulse(n),'beta') || isempty(Exp.Pulse(n).beta)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify dimensionless beta parameter for tanh.']);
        end
        if ~isfield(Exp.Pulse(n),'SweepDirection') || isempty(Exp.Pulse(n).SweepDirection)
          Exp.Pulse(n).SweepDirection = +1;
        end
        
      case 'BW compensated'
        
        if (~isfield(Exp.Pulse(n),'freqaxis') || isempty(Exp.Pulse(n).freqaxis)) || ...
            (~isfield(Exp.Pulse(n),'v1') || isempty(Exp.Pulse(n).v1))
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify the resonator magnitude response function.']);
        end
        if ~isfield(Exp,'mwFreq') || isempty(Exp.mwFreq)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Exp.mwFreq is required to compute resonator bandwidth compensation.']);
        end
        
      case 'uniform adiabaticity'
        
        if ~isfield(Exp.Pulse(n),'BW') || isempty(Exp.Pulse(n).BW)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify BW parameter in MHz for tanh.']);
        end
        
      otherwise
        
        error(['The frequency modulation function selected for pulse ',num2str(n),' is not defined.']);
        
    end
    
    % Set up time axis
    if ~isfield(Exp,'timestep') || (isfield(Exp,'timestep') && numel(Exp.timestep)<n)
      if ~strcmp(FrequencyModulation,'none')
        % Compute time step for frequency modulated pulses based on the
        % Nyquist sampling theorem considering the maximum frequency
        % of the pulse and the frequency offset
        maxFreq = max(abs([Exp.FreqOffset(n)-Exp.Pulse(n).BW/2 ...
          Exp.FreqOffset(n)+Exp.Pulse(n).BW/2]));
        Exp.timestep(n) = 1/(2*Opt.OverSampleFactor*maxFreq);
      elseif numel(AmplitudeModulation)==1 && strcmp(AmplitudeModulation,'rectangular')
        maxFreq = Exp.FreqOffset(n);
        if maxFreq~=0
          Exp.timestep(n) = 1/(2*Opt.OverSampleFactor*maxFreq);
        else
          Exp.timestep(n) = 0.002; % ns, default time step
        end
      else
        % Determine "equivalent maximum frequencies" for AM functions
        equivmaxFreq = zeros(1,numel(AmplitudeModulation));
        for na = 1:numel(AmplitudeModulation)
          if strcmp(AmplitudeModulation,'gaussian')
            equivmaxFreq(na) = 1/(2*Exp.Pulse(n).tFWHM);
          elseif strcmp(AmplitudeModulation,'sinc')
            equivmaxFreq(na) = 1/Exp.Pulse(n).zerocross;
          elseif strcmp(AmplitudeModulation,'sech')
            % from approximate FWHM of sech function
            equivmaxFreq(na) = Exp.Pulse(n).beta/(4*Exp.tp(n)*asech(0.5));
          elseif strcmp(AmplitudeModulation,'nth order sech')
            % from rise time to 1/2 of maximum value
            thalf = (Exp.tp(n)/2)*(1-(2*asech(0.5)/Exp.Pulse(n).beta)^(1/Exp.Pulse(n).n));
            equivmaxFreq(na) = 1/(4*thalf);
          elseif strcmp(AmplitudeModulation,'WURST')
            % from rise time to 1/2 of maximum value
            thalf = Exp.tp(n)*(1/2 - (1/pi)*asin(2^(-1/Exp.Pulse(n).n)));
            equivmaxFreq(na) = 1/(4*thalf);
          elseif strcmp(AmplitudeModulation,'quartersin')
            equivmaxFreq(na) = 1/(4*Exp.Pulse(n).trise);
          end
        end
        maxFreq = max(abs([equivmaxFreq Exp.FreqOffset]));
        Exp.timestep(n) = 1/(2*Opt.OverSampleFactor*maxFreq);
      end
      Exp.timestep(n) = Exp.tp(n)/round(Exp.tp(n)/Exp.timestep(n)); % last time point = tp
    end
    t{n} = 0:Exp.timestep(n):Exp.tp(n);
    ti = t{n} - Exp.tp(n)/2;
    
    % ------------------------------------------------------------------- %
    % Amplitude modulation function
    % ------------------------------------------------------------------- %
    mod(n).A = ones(1,numel(t{n}));
    A = zeros(numel(AmplitudeModulation),numel(t{n}));
    for na = 1:numel(AmplitudeModulation)
      switch AmplitudeModulation{na}
        
        case 'rectangular'
          
          A(na,1:numel(t{n})) = 1;
          
        case 'gaussian'
          
          A(na,:) = exp(-(4*log(2)*ti.^2)/Exp.Pulse(n).tFWHM^2);
          
        case 'sinc'
          
          A(na,:) = sin((2*pi*ti)/Exp.Pulse(n).zerocross)./((2*pi*ti)/Exp.Pulse(n).zerocross);
          A(na,Exp.tp(n)/(2*Exp.timestep(n))+1) = 1;
          A(na,:) = A(na,:)/max(A(na,:));
          
        case 'quartersin'
          
          % Pulse edges weighted with a quarter period of a sine wave
          A(na,1:numel(t{n})) = 1;
          if Exp.Pulse(n).trise~=0 && 2*Exp.Pulse(n).trise<Exp.tp(n)
            npts = Exp.Pulse(n).trise/Exp.timestep(n);
            tpartial = 0:Exp.timestep(n):Exp.Pulse(n).trise;
            A(na,1:npts+1) = sin(tpartial*(pi/(2*Exp.Pulse(n).trise)));
            A(na,end-npts:end) = sin((pi/2)+tpartial*(pi/(2*Exp.Pulse(n).trise)));
          end
          
        case 'sech'
          
          A(na,:) = sech((Exp.Pulse(n).beta/Exp.tp(n))*ti);
          
        case 'nth order sech'
          
          A(na,:) = sech(Exp.Pulse(n).beta*2^(Exp.Pulse(n).n-1)*(ti/Exp.tp(n)).^Exp.Pulse(n).n);
          
        case 'WURST'
          
          A(na,:) = (1-(abs(sin(pi*ti/Exp.tp(n)))).^Exp.Pulse(n).n);
          
      end
      
      mod(n).A = mod(n).A.*A(na,:);
    end
    
    % ------------------------------------------------------------------- %
    % Frequency (mod.nu) and phase (mod.phi) modulation functions
    % ------------------------------------------------------------------- %
    switch FrequencyModulation
      
      case 'none'
        
        mod(n).nu(1:numel(t{n})) = 0;
        mod(n).phi(1:numel(t{n})) = 0;
        
      case 'linear'
        
        k = Exp.Pulse(n).BW/Exp.tp(n); % rate of change
        mod(n).nu = -Exp.Pulse(n).BW/2+k*t{n};
        mod(n).phi = 2*pi*(-Exp.Pulse(n).BW/2*t{n}+(k/2)*t{n}.^2);
        
      case 'BW-compensated'
        
        % Variable-rate chirps with resonator bandwidth
        % compensation as described in:
        % Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
        % Adiabatic and fast passage ultra-wideband inversion in
        % pulsed EPR. J. Magn. Reson. 230, 27–39 (2013).
        
        % Constant rate chirp
        k = Exp.Pulse(n).BW/Exp.tp(n); % rate of change
        nu_linear = -Exp.Pulse(n).BW/2+k*t{n};
        
        v1_range = interp1(Exp.Pulse(n).freqaxis,Exp.Pulse(n).v1,...
          nu_linear+Exp.FreqOffset+Exp.mwFreq);
        
        % Frequency dependence of t and time to frequency mapping
        const = trapz(1./v1_range.^2)/t{n}(end); % const = 2pi/Qref,
        % Qref = reference adiabaticity
        t_f = cumtrapz((1/const)*(1./v1_range.^2));
        nu_adapted = interp1(t_f,nu_linear,t{n},'pchip');
        
        mod(n).nu = nu_adapted;
        mod(n).phi = 2*pi*cumtrapz(t{n},mod(n).nu);
        
      case 'tanh'
        
        % Determine BWinf parameter from BW and beta parameters
        % (the frequency is swept from -BW/2 to +BW/2)
        Exp.Pulse(n).BWinf = Exp.Pulse(n).BW/tanh(Exp.Pulse(n).beta/2);
        
        mod(n).nu = (Exp.Pulse(n).BWinf/2)*tanh((Exp.Pulse(n).beta/Exp.tp(n))*ti);
        mod(n).phi = (Exp.Pulse(n).BWinf/2)*(Exp.tp(n)/Exp.Pulse(n).beta)*...
          log(cosh((Exp.Pulse(n).beta/Exp.tp(n))*ti));
        mod(n).phi = 2*pi*Exp.Pulse(n).SweepDirection*mod(n).phi;
        
      case 'uniform adiabaticity'
        % The frequency modulation is calculated as the integral of the
        % squared amplitude modulaton function (for nth order sech/tanh or
        % in general to obtain offset-independent adiabaticity pulses, see
        % Garwood, M., DelaBarre, L., J. Magn. Reson. 153, 155-177 (2001).
        
        mod(n).nu = cumtrapz(ti,mod(n).A.^2)/trapz(ti,mod(n).A.^2); % F2
        mod(n).nu = Exp.Pulse(n).BW*(mod(n).nu-1/2);
        mod(n).phi = 2*pi*cumtrapz(ti,mod(n).nu);
        mod(n).phi = mod(n).phi+abs(min(mod(n).phi)); % zero phase offset at pulse center
        
    end
    
    % ------------------------------------------------------------------- %
    % Determine pulse amplitude from flip angle (if only Exp.Flip is given)
    % ------------------------------------------------------------------- %
    if (isfield(Exp,'Flip') && ~isempty(Exp.Flip(n))) && ...
        (~isfield(Exp,'Amplitude') || (isfield(Exp,'Amplitude') && numel(Exp.Amplitude)<n))
      switch FrequencyModulation
        
        case 'none' % amplitude modulated pulses: beta = integral
          
          Exp.Amplitude(n) = Exp.Flip(n)/(2*pi*trapz(t{n},mod(n).A));
          
        case {'linear','BW compensated','tanh','uniform adiabaticity'}
          % see Jeschke et al. (2015) J. Phys. Chem. B, 119, 13570–13582.
          % Q_crit = (2*pi*v1max)^2/k = minimum adiabaticity on resonance
          
          if Exp.Flip(n)>pi
            error('Pulse amplitude calculation from flip angle not applicable for angles larger than pi.')
          end
          Q_crit = (2/pi)*(log(2/(1+cos(Exp.Flip(n)))));
          if Q_crit>5 % set Q_crit to finite value if it is infinite or very large
            Q_crit = 5;
          end
          
          if strcmp(FrequencyModulation,'linear') || strcmp(FrequencyModulation,'BW compensated')
            sweeprate = Exp.Pulse(n).BW/Exp.tp(n);
          elseif strcmp(FrequencyModulation,'tanh')
            sweeprate = Exp.Pulse(n).beta*Exp.Pulse(n).BWinf;
          elseif strcmp(FrequencyModulation,'uniform adiabaticity')
            sweeprate = Exp.Pulse(n).beta*Exp.Pulse(n).BW;
          end
          
          Exp.Amplitude(n) = sqrt(2*pi*Q_crit*sweeprate)/(2*pi);
          
      end
    end
    
    % ------------------------------------------------------------------- %
    % Pulse function
    % ------------------------------------------------------------------- %
    mod(n).A = Exp.Amplitude(n)*mod(n).A;
    y{n} = mod(n).A.*cos(mod(n).phi + 2*pi*Exp.FreqOffset(n)*t{n} + Exp.Phase(n))...
      +1i*mod(n).A.*cos(mod(n).phi + 2*pi*Exp.FreqOffset(n)*t{n} + Exp.Phase(n) - (pi/2));
    
  end
  
  % --------------------------------------------------------------------- %
  % Excitation profile calculation
  % --------------------------------------------------------------------- %
  if Opt.ExciteProfile==1
    
    % Spin operators
    Sx = sop(1/2,'x');
    Sy = sop(1/2,'y');
    Sz = sop(1/2,'z');
    
    p0 = -Sz;
    
    p(n).Mz = zeros(1,numel(Opt.offsets));
    for k = 1:length(Opt.offsets)
      
      Ham0 = Opt.offsets(k)*Sz;
      
      p1 = p0;
      
      if min(y{n})==max(y{n}) % used for rectangular pulses
        Ham = real(y{n}(1))*Sx+imag(y{n}(1))*Sy+Ham0;
        U = expm(-2i*pi*Ham*Exp.timestep(n));
        for j = 1:numel(t{n})
          p1 = U*p1*U';
        end
      else
        if Opt.nBCH==1 % exact solution
          for j = 1:numel(t{n})
            
            Ham = real(y{n}(j))*Sx+imag(y{n}(j))*Sy+Ham0;
            U = expm(-2i*pi*Ham*Exp.timestep(n));
            p1 = U*p1*U';
          end
        else % average Hamiltonian approximation (Baker-Campbell-Hausdorff series)
        % see Ernst, Bodenhausen, Wokaun, Principles of NMR in One and Two
        % Dimensions, Clarendon Press (1990), Chapter 3, p. 72-75.
          ll = 0;
          for j = 1:numel(t{n})/Opt.nBCH
            
            Ham_BCH = zeros(size(Ham0));
            
            for l1 = 1:Opt.nBCH
              Ham_BCH = Ham_BCH -2i*pi*Exp.timestep(n)*(real(y{n}(ll+l1))*Sx+imag(y{n}(ll+l1))*Sy+Ham0);
              for l2 = l1+1:Opt.nBCH
                Ham_BCH = Ham_BCH -...
                  2*(pi*Exp.timestep(n))^2*(...
                  (real(y{n}(ll+l2))*Sx+imag(y{n}(ll+l2))*Sy+Ham0)*(real(y{n}(ll+l1))*Sx+imag(y{n}(ll+l1))*Sy+Ham0)-...
                  (real(y{n}(ll+l1))*Sx+imag(y{n}(ll+l1))*Sy+Ham0)*(real(y{n}(ll+l2))*Sx+imag(y{n}(ll+l2))*Sy+Ham0));
              end
            end
            
            ll = ll+Opt.nBCH;
            U = expm(Ham_BCH);
            p1 = U*p1*U';
          end
          
        end
      end
      
      p(n).Mz(k) = -sum(sum(Sz.*p1.'));
    end
    
    if isfield(Opt,'plot') && Opt.plot==1
      clf
      subplot(2,1,1)
      hold on; box on;
      plot(t{n},real(y{n}),'b')
      plot(t{n},imag(y{n}),'r')
      xlabel('t [\mu s]')
      ylabel('\nu_1 [MHz]')
      axis tight
      subplot(2,1,2)
      hold on; box on;
      plot(Opt.offsets,real(p(n).Mz),'b');
      xlabel('Frequency offset [MHz]')
      ylabel('M_z')
      axis tight
      ylim([-0.5 0.5])
    end
    
  else
    
    if Opt.plot==1
      clf
      hold on; box on;
      plot(t{n},real(y{n}),'b')
      plot(t{n},imag(y{n}),'r')
      xlabel('t [\mu s]')
      ylabel('\nu_1 [MHz]')
      axis tight
    end
    
  end
  
end

% ----------------------------------------------------------------------- %
% Output
% ----------------------------------------------------------------------- %
if numel(npulse)==1
  t = t{npulse};
  y = y{npulse};
  mod = mod(npulse);
end
if nargout==1
  error('The function pulse needs to be called with at least two output arguments.')
elseif nargout==2 % [t,y] = pulse(...)
  varargout{1} = t;
  varargout{2} = y;
elseif nargout==3 % [t,y,mod] = pulse(...) or [t,y,p] = pulse(...)
  varargout{1} = t;
  varargout{2} = y;
  if Opt.ExciteProfile==1;
    varargout{3} = p;
  else
    varargout{3} = mod;
  end
elseif nargout==4 % [t,y,mod,p] = pulse(...)
  if Opt.ExciteProfile==0;
    error('The function pulse returns only up to three output arguments for the selected options.')
  end
  varargout{1} = t;
  varargout{2} = y;
  varargout{3} = mod;
  varargout{4} = p;
end
