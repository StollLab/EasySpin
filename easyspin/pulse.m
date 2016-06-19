function varargout = pulse(varargin)

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% Pulse definition and excitation profile calculation
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%
% [t,y] = pulse(Exp)
% [t,y] = pulse(Exp,iPulse)
% [t,y] = pulse(Exp,Opt)
% [t,y] = pulse(Exp,iPulse,Opt)
%
% [t,y,p] = pulse(Exp,iPulse,Opt)
% [t,y,p,modulation] = pulse(Exp,iPulse,Opt)
%
% Input: 
% Exp = structure containing the following fields:
%       Exp.tp              = pulse length in us
%       Exp.TimeStep        = time step for waveform definition in us
%                         (default = determined based on input pulse
%                          parameters)
%       Exp.Amplitude       = pulse amplitude in MHz
%       Exp.Flip            = pulse flip angle in radians (*1)
%        (if Exp.Amplitude and Exp.Flip are both missing, a pulse with
%         amplitude 1 is returned; if both Exp.Amplitude and Exp.Flip are
%         given, Exp.Flip is ignored)
%       Exp.CenterFreq      = pulse center frequency (default = 0)
%       Exp.Phase           = phase for the pulse in radians (default = 0,
%                             i.e. +x)
%       Exp.PulseShape.Type = pulse shape name in a string with structure
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
%                            (default = 'rectangular')
%       Exp.PulseShape.*   = value for the pulse parameters defining
%                            the specified pulse shape, the pulse
%                            parameters required for each of the available
%                            modulation functions are listed below.
%       Exp.PulseShape.I/Q = I and Q data describing an arbitrary pulse,
%                            the time axis is reconstructed based on Exp.tp
%                            and the length of the I and Q vectors, all
%                            other input parameters (Amplitude, Flip,
%                            CenterFreq, Phase, etc.) are ignored
%
% iPulse = index/indices for the pulses in Exp input structure for
%                 which to return the pulse data
%
% Opt = optional structure with the following fields
%       Opt.OverSampleFactor = oversampling factor for the determination
%                              of the time step (default = 10)
%       Opt.Detect           = 'Sz' (default), 'Sy', 'Sx', 'all', define
%                               detection operator for excitation profile
%                               calculation
%       Opt.Offsets          = axis of frequency offsets in MHz for which
%                              to compute the excitation profile
%                              (default approximately ±1.5*BW centered at
%                               Exp.CenterFreq, 201 points)
%
% Available pulse modulation functions:
%   - Amplitude modulation: rectangular, gaussian, sinc, quartersin, sech, WURST
%   - Frequency modulation: none, linear, tanh, BWcompensated, uniformQ
%
% The parameters required for the different modulation functions are:
% Amplitude modulation:
% 'rectangular'         - none
% 'gaussian'            - tFWHM     = FWHM in us
%                       Alternatively:
%                       - trunc     = truncation parameter (0 to 1)
% 'sinc'                - zerocross = width between the first zero-
%                                     crossing points in us
% 'sech'                - beta      = dimensionless truncation parameter
%                       - n         = order of the secant function argument
%                                     (default = 1)
% 'WURST'               - nwurst    = WURST n parameter (determining the
%                                     steepness of the amplitude function)
% 'quartersin'          - trise     = rise time in us for quarter sine
%                                     weighting at the pulse edges
%
% Frequency modulation:
% 'linear'              - BW        = frequency sweep width in MHz
%                                     (the frequency is swept over ±BW/2)
%                       - SweepDirection = +1 or -1 (+1 default)
% 'tanh'                - BW        = frequency sweep width in MHz
%                       - beta      = dimensionless truncation parameter
%                       - SweepDirection = +1 or -1 (+1 default)
% 'BWcompensated'(*2)   - BW        = frequency sweep width in MHz
%                                     (the frequency is swept over ±BW/2)
%                       - SweepDirection = +1 or -1 (+1 default)
%                       Parameters required for resonator bandwidth
%                       compensation:
%                       - freqaxis  = frequency axis
%                       - v1        = magnitude response function (ideal or
%                                     experimental)
%                       - Exp.mwFreq needs to be defined
% 'uniformQ'            - BW        = frequency sweep width in MHz
%                       the frequency modulation is calculated as the
%                       integral of the squared amplitude modulation
%                       function (for nth order sech/tanh pulses or in
%                       general to obtain offset-independent adiabaticity
%                       pulses, see Garwood, M., DelaBarre, L., J. Magn.
%                       Reson. 153, 155-177 (2001).
%
% Output:   t          = time axis for defined waveform in us
%           y          = real and imaginary part of the pulse function
%           If three or four output arguments are requested, the excitation
%           profile is calculated. Additional output fields are:
%           p          = structure with frequency offset axis (p.offsets) 
%                        for excitation profile in MHz and excitation
%                        profile (Mi/M0, i = x,y,z, p.Mz/p.My/p.Mx)
%           modulation = structure with amplitude (modulation.A in MHz),
%                        frequency (modulation.nu in MHz) and phase
%                        (modulation.phi in rad) modulation functions
%
% (*1) The conversion from flip angles to amplitudes is performed using the
%      approximations described in:
%      Jeschke, G., Pribitzer, S., Doll, A. Coherence Transfer by Passage
%      Pulses in Electron Paramagnetic Resonance Spectroscopy.
%      J. Phys. Chem. B 119, 13570–13582 (2015).
% (*2) Chirps with variable rate to compensate for the resonator bandwidth.
%      The bandwidth compensation is implemented as described in:
%      Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., Adiabatic and
%      fast passage ultra-wideband inversion in pulsed EPR
%      J. Magn. Reson. 230, 27-39 (2013).
%

% ----------------------------------------------------------------------- %
% Input argument parsing
% ----------------------------------------------------------------------- %
switch nargin
  case 1 % [t,y] = pulse(Exp)
    Exp = varargin{1};
    Opt = struct;
  case 2 % [t,y] = pulse(Exp,iPulse) or [t,y] = pulse(Exp,Opt)
    Exp = varargin{1};
    if isstruct(varargin{2})
      Opt = varargin{2};
    else
      Opt = struct;
      iPulse = varargin{2};
    end
  case 3 % [t,y] = pulse(Exp,iPulse,Opt)
    Exp = varargin{1};
    iPulse = varargin{2};
    Opt = varargin{3};
  otherwise
    error('The function pulse is not supported for more than 3 input arguments.')
end

plotResults = (nargout==0);

if nargout==0
  Opt.ExciteProfile = true;
end
if (nargout==3 || nargout==4) && ~isfield(Opt,'ExciteProfile')
  Opt.ExciteProfile = true;
end
if ~isfield(Opt,'ExciteProfile')
  Opt.ExciteProfile = false;
end

% Set experiment and option parameters to defaults
if ~isfield(Exp,'tp')
  error('Pulse length not defined in Exp.tp.')
end
% If no pulse is selected, the functions returns outputs for each of the
% pulses defined in the input structure
if ~exist('iPulse','var')
  iPulse = 1:numel(Exp.tp);
end
if ~isfield(Exp,'Amplitude') && ~isfield(Exp,'Flip')
  Exp.Amplitude(1:numel(Exp.tp)) = 1; % normalized amplitude
elseif isfield(Exp,'Amplitude') && numel(Exp.Amplitude)~=numel(Exp.tp)
  Exp.Amplitude(end:numel(Exp.tp)) = Exp.Amplitude(1);
end
if isfield(Exp,'Flip') && numel(Exp.Flip)~=numel(Exp.tp)
  Exp.Flip(end:numel(Exp.tp)) = Exp.Flip(end);
end
if isfield(Exp,'Flip') && numel(Exp.tp)~=numel(Exp.Flip)
  error('The number of pulses in Exp.tp and Exp.Flip/Exp.Amplitude do not agree.')
end
if ~isfield(Exp,'CenterFreq')
  Exp.CenterFreq(1:numel(Exp.tp)) = 0; % MHz
elseif numel(Exp.CenterFreq)~=numel(Exp.tp)
  Exp.CenterFreq(end:numel(Exp.tp)) = Exp.CenterFreq(end);
end
if ~isfield(Exp,'Phase')
  Exp.Phase(1:numel(Exp.tp)) = 0; % rad
elseif numel(Exp.Phase)~=numel(Exp.tp)
  Exp.Phase(end:numel(Exp.tp)) = Exp.Phase(end);
end
if isfield(Exp,'TimeStep') && numel(Exp.TimeStep)~=numel(Exp.tp)
  Exp.TimeStep(end:numel(iPulse)) = Exp.TimeStep(end);
end
if ~isfield(Exp,'PulseShape')
  for i = 1:numel(Exp.tp)
    Exp.PulseShape(i).Type = 'rectangular';
  end
end
if numel(Exp.PulseShape)<numel(Exp.tp)
  for i = numel(Exp.PulseShape)+1:numel(Exp.tp)
    Exp.PulseShape(i).Type = 'rectangular';
  end
end

% Options
% ----------------------------------------------------------------------- %
if ~isfield(Opt,'OverSampleFactor')
  Opt.OverSampleFactor = 10;
end
if ~isfield(Opt,'Detect')
  Opt.Detect = 'Sz';
end
if ~isfield(Opt,'nOffsets')
  Opt.nOffsets = 201;
end

% ----------------------------------------------------------------------- %
% Loop over pulses defined in input structure and calculate pulse function
% ----------------------------------------------------------------------- %
t = cell(1,numel(iPulse));
y = cell(1,numel(iPulse));
modulation(1:numel(iPulse)) = struct('A',[]);
p(1:numel(iPulse)) = struct('offsets',[]);

for n = iPulse
  thisPulse = Exp.PulseShape(n);
  
  % Check if pulse I and Q data is given
  if (isfield(thisPulse,'I') && ~isempty(thisPulse.I)) || ...
     (isfield(thisPulse,'Q') && ~isempty(thisPulse.Q))
    
    if ~isfield(thisPulse,'Type') || isempty(thisPulse.Type)
      thisPulse.Type = 'user-IQ';
    end
    
    if ~isfield(thisPulse,'I')
      thisPulse.I = zeros(size(thisPulse.Q));
    elseif ~isfield(thisPulse,'Q')
      thisPulse.Q = zeros(size(thisPulse.I));
    elseif size(thisPulse.I)~=size(thisPulse.Q);
      error(['I and Q input data for pulse ',num2str(n),' have different lengths.'])
    end
    if ~isvector(thisPulse.I) || ~isvector(thisPulse.Q)
      error(['I and Q input data for pulse ',num2str(n),' should be vectors.'])
    end
    
    t{n} = linspace(0,Exp.tp(n),numel(thisPulse.I));
    nPoints = numel(t{n});
    y{n} = thisPulse.I + 1i*thisPulse.Q;
    
    Exp.TimeStep(n) = t{n}(2)-t{n}(1);
    modulation(n).A = [];
    modulation(n).dnu = [];
    modulation(n).phi = [];
    
    if Opt.ExciteProfile && ~isfield(Opt,'Offsets')
      error('The excitation profile of a user-defined pulse can only be calculated if the desired range of frequency offsets is specified in Opt.Offsets.');
    end
    
  else
    
    % Set pulse shape to rectangular if it is not specified
    if ~isfield(thisPulse,'Type') || isempty(thisPulse.Type) || strcmp(thisPulse.Type,'')
      thisPulse.Type = 'rectangular';
    end
    
    % Determine pulse shape from input string
    shape = regexp(thisPulse.Type,'/','split');
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
          
          if (~isfield(thisPulse,'tFWHM') || isempty(thisPulse.tFWHM)) && ...
             (~isfield(thisPulse,'trunc') || isempty(thisPulse.trunc))
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify tFWHM or trunc parameter for the Gaussian envelope.']);
          elseif ~isfield(thisPulse,'tFWHM') && isfield(thisPulse,'trunc') && ...
              ~isempty(thisPulse.trunc)
            % Convert truncation parameter to tFWHM
            thisPulse.tFWHM = sqrt(-(Exp.tp(n)^2)/log2(thisPulse.trunc));
            if thisPulse.tFWHM==0
              thisPulse.tFWHM = Exp.TimeStep(n)/2;
            end
          end
          
        case 'sinc'
          
          if ~isfield(thisPulse,'zerocross') || isempty(thisPulse.zerocross)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify zero-crossing parameter in us for the sinc envelope.']);
          end
          
        case 'quartersin'
          
          if ~isfield(thisPulse,'trise') || isempty(thisPulse.trise)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify trise parameter in us for the quartersine envelope.']);
          end
          
        case 'sech'
          
          if ~isfield(thisPulse,'beta') || isempty(thisPulse.beta)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify beta parameter in 1/us for the sech envelope.']);
          end
          if ~isfield(thisPulse,'n') || isempty(thisPulse.n)
            thisPulse.n = 1;
          end
          
        case 'WURST'
          
          if ~isfield(thisPulse,'nwurst') || isempty(thisPulse.nwurst)
            error(['Pulse AM function of pulse ',num2str(n),' not sufficiently defined. ',...
              'Specify nwurst parameter for the WURST envelope.']);
          end
          
        otherwise
          
          error(['The amplitude modulation function selected for pulse ',num2str(n),' is not defined.']);
          
      end
    end
    
    switch FrequencyModulation
      
      case 'none'
        
      case 'linear'
        
        if ~isfield(thisPulse,'BW') || isempty(thisPulse.BW)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify BW parameter in MHz for the linear chirp.']);
        end
        
      case 'tanh'
        
        if ~isfield(thisPulse,'BW') || isempty(thisPulse.BW)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify BW parameter in MHz for tanh.']);
        end
        if ~isfield(thisPulse,'beta') || isempty(thisPulse.beta)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify dimensionless beta parameter for tanh.']);
        end
        
      case 'BWcompensated'
        
        if (~isfield(thisPulse,'freqaxis') || isempty(thisPulse.freqaxis)) || ...
            (~isfield(thisPulse,'v1') || isempty(thisPulse.v1))
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify the resonator magnitude response function.']);
        end
        if ~isfield(Exp,'mwFreq') || isempty(Exp.mwFreq)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Exp.mwFreq is required to compute resonator bandwidth compensation.']);
        end
        
      case 'uniformQ'
        
        if ~isfield(thisPulse,'BW') || isempty(thisPulse.BW)
          error(['Pulse FM function of pulse ',num2str(n),' not sufficiently defined. ',...
            'Specify BW parameter in MHz.']);
        end
        
      otherwise
        
        error(['The frequency modulation function selected for pulse ',num2str(n),' is not defined.']);
        
    end
    if ~isfield(thisPulse,'SweepDirection') || isempty(thisPulse.SweepDirection)
      thisPulse.SweepDirection = +1;
    end
    
    % Set up time axis
    if ~isfield(Exp,'TimeStep') || (isfield(Exp,'TimeStep') && numel(Exp.TimeStep)<n)
      
      % Automatically determine appropriate time step
      %------------------------------------------------------------------ %
      % Determine bandwidth of frequency modulation
      switch FrequencyModulation
        case 'none'
          FM_BW = 0;
        otherwise
          FM_BW = thisPulse.BW;
      end
      
      % Determine bandwidth of amplitude modulation
      AM_BW = 0;
      for na = 1:numel(AmplitudeModulation)
        switch AmplitudeModulation{na}
          case 'rectangular'
            AM_BW = AM_BW + 4/Exp.tp(n);
          case 'gaussian'
            AM_BW = AM_BW + 1/(2*thisPulse.tFWHM);
          case 'sinc'
            AM_BW = AM_BW + 1/thisPulse.zerocross;
          case 'sech'
            % from approximate FWHM of sech function
            AM_BW = AM_BW + thisPulse.beta/(4*Exp.tp(n)*asech(0.5));
            %case 'nth order sech'
            % from rise time to 1/2 of maximum value
            %thalf = (Exp.tp(n)/2)*(1-(2*asech(0.5)/thisPulse.beta)^(1/thisPulse.n));
            %AM_BW = AM_BW + 1/(4*thalf);
          case 'WURST'
            % from rise time to 1/2 of maximum value
            thalf = Exp.tp(n)*(1/2 - (1/pi)*asin(2^(-1/thisPulse.nwurst)));
            AM_BW = AM_BW + 1/(4*thalf);
          case 'quartersin'
            AM_BW = AM_BW + 1/(4*thisPulse.trise);
        end
      end
      
      % Calculate maximum frequency offset
      BW = max([FM_BW AM_BW]);
      maxFreq = max(abs(Exp.CenterFreq(n)+[-1 1]*BW/2));
      % Use Nyquist theorem to calculate time step, incl. oversampling
      if maxFreq~=0
        Nyquist_dt = 1/(2*maxFreq);
        Exp.TimeStep(n) = Nyquist_dt/Opt.OverSampleFactor;
      else
        Exp.TimeStep(n) = 0.002; % us
      end
      if Exp.TimeStep(n)>Exp.tp(n)
        Exp.TimeStep(n) = Exp.tp(n);
      end
      Exp.TimeStep(n) = Exp.tp(n)/round(Exp.tp(n)/Exp.TimeStep(n)); % last time point = tp
    end
    t{n} = 0:Exp.TimeStep(n):Exp.tp(n);
    ti = t{n} - Exp.tp(n)/2;
    nPoints = numel(t{n});
    
    % ------------------------------------------------------------------- %
    % Amplitude modulation function (modulation.A)
    % ------------------------------------------------------------------- %
    modulation(n).A = ones(1,nPoints);
    for na = 1:numel(AmplitudeModulation)
      switch AmplitudeModulation{na}
        
        case 'rectangular'
          
          A = ones(1,nPoints);
          
        case 'gaussian'
          
          A = exp(-(4*log(2)*ti.^2)/thisPulse.tFWHM^2);
          
        case 'sinc'
          
          A = sin((2*pi*ti)/thisPulse.zerocross)./((2*pi*ti)/thisPulse.zerocross);
          A(round(Exp.tp(n)/(2*Exp.TimeStep(n)))+1) = 1;
          A = A/max(A);
          
        case 'quartersin'
          
          % Pulse edges weighted with a quarter period of a sine wave
          A = ones(1,nPoints);
          if thisPulse.trise~=0 && 2*thisPulse.trise<Exp.tp(n)
            tpartial = 0:Exp.TimeStep(n):thisPulse.trise;
            npts = numel(tpartial);
            A(1:npts) = sin(tpartial*(pi/(2*thisPulse.trise)));
            A(end-npts+1:end) = fliplr(A(na,1:npts));
          end
          
        case 'sech'
          
          if thisPulse.n==1 % reduces numerical errors
            A = sech((thisPulse.beta/Exp.tp(n))*ti);
          else
            A = sech(thisPulse.beta*2^(thisPulse.n-1)*(ti/Exp.tp(n)).^thisPulse.n);
          end
          
        case 'WURST'
          
          A = (1-(abs(sin(pi*ti/Exp.tp(n)))).^thisPulse.nwurst);
          
      end
      modulation(n).A = modulation(n).A.*A;
    end
    
    % ------------------------------------------------------------------- %
    % Frequency (modulation.nu) and phase (modulation.phi) modulation functions
    % ------------------------------------------------------------------- %
    switch FrequencyModulation
      
      case 'none'
        
        modulation(n).nu(1:nPoints) = 0;
        modulation(n).phi(1:nPoints) = 0;
        
      case 'linear'
        
        k = thisPulse.BW/Exp.tp(n); % frequency sweep (chirp) rate
        modulation(n).nu = k*ti;
        modulation(n).phi = 2*pi*thisPulse.SweepDirection*((k/2)*ti.^2);
        
      case 'BWcompensated'
        
        % Variable-rate chirps with resonator bandwidth
        % compensation, as described in:
        % Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
        % Adiabatic and fast passage ultra-wideband inversion in
        % pulsed EPR. J. Magn. Reson. 230, 27–39 (2013).
        % http://dx.doi.org/10.1016/j.jmr.2013.01.002
        
        % Constant rate chirp
        k = thisPulse.BW/Exp.tp(n); % rate of change
        nu_linear = k*ti;
        
        v1_range = interp1(thisPulse.freqaxis,thisPulse.v1,...
          nu_linear+Exp.CenterFreq(n)+Exp.mwFreq*10^3);
        
        % Frequency dependence of t and time to frequency mapping
        const = trapz(1./v1_range.^2)/t{n}(end); % const = 2pi/Qref,
        % Qref = reference adiabaticity
        t_f = cumtrapz((1/const)*(1./v1_range.^2));
        nu_adapted = interp1(t_f,nu_linear,t{n},'pchip');
        
        modulation(n).nu = nu_adapted;
        modulation(n).phi = 2*pi*thisPulse.SweepDirection*cumtrapz(t{n},modulation(n).nu);
        modulation(n).phi = modulation(n).phi+abs(min(modulation(n).phi));
        
      case 'tanh'
        
        % Determine BWinf parameter from BW and beta parameters
        % (the frequency is swept from -BW/2 to +BW/2)
        thisPulse.BWinf = thisPulse.BW/tanh(thisPulse.beta/2);
        
        modulation(n).nu = (thisPulse.BWinf/2)*tanh((thisPulse.beta/Exp.tp(n))*ti);
        modulation(n).phi = (thisPulse.BWinf/2)*(Exp.tp(n)/thisPulse.beta)*...
          log(cosh((thisPulse.beta/Exp.tp(n))*ti));
        modulation(n).phi = 2*pi*thisPulse.SweepDirection*modulation(n).phi;
        
      case 'uniformQ'
        % The frequency modulation is calculated as the integral of the
        % squared amplitude modulation function (for nth order sech/tanh or
        % in general to obtain offset-independent adiabaticity pulses, see
        % Garwood, M., DelaBarre, L., J. Magn. Reson. 153, 155-177 (2001)
        % http://dx.doi.org/10.1006/jmre.2001.2340
        
        modulation(n).nu = cumtrapz(ti,modulation(n).A.^2)/trapz(ti,modulation(n).A.^2); % F2
        modulation(n).nu = thisPulse.BW*(modulation(n).nu-1/2);
        modulation(n).phi = 2*pi*cumtrapz(ti,modulation(n).nu);
        modulation(n).phi = modulation(n).phi+abs(min(modulation(n).phi)); % zero phase offset at pulse center
        modulation(n).phi = thisPulse.SweepDirection*modulation(n).phi;
        
    end
    
    % ------------------------------------------------------------------- %
    % Determine pulse amplitude from flip angle (if only Exp.Flip is given)
    % ------------------------------------------------------------------- %
    if (isfield(Exp,'Flip') && ~isempty(Exp.Flip(n))) && ...
       (~isfield(Exp,'Amplitude') || (isfield(Exp,'Amplitude') && numel(Exp.Amplitude)<n))
      switch FrequencyModulation
        
        case 'none' % amplitude modulated pulses: beta = integral
          
          Exp.Amplitude(n) = Exp.Flip(n)/(2*pi*trapz(t{n},modulation(n).A));
          
        case {'linear','BWcompensated','tanh','uniformQ'}
          % see Jeschke et al. (2015) J. Phys. Chem. B, 119, 13570–13582.
          % Q_crit = (2*pi*v1max)^2/k = minimum adiabaticity on resonance
          
          if Exp.Flip(n)>pi
            error('Pulse amplitude calculation from flip angle not applicable for angles larger than pi.')
          end
          Q_crit = (2/pi)*log(2/(1+cos(Exp.Flip(n))));
          if Q_crit>5 % set Q_crit to finite value if it is infinite or very large
            Q_crit = 5;
          end
          
          if strcmp(FrequencyModulation,'linear') || strcmp(FrequencyModulation,'BWcompensated')
            sweeprate = thisPulse.BW/Exp.tp(n);
          elseif strcmp(FrequencyModulation,'tanh')
            sweeprate = thisPulse.beta*thisPulse.BWinf/(2*Exp.tp(n));
          elseif strcmp(FrequencyModulation,'uniformQ')
            % Q = w1max^2*A(t)^2/(BW*dnu/dt) see eq. 17 in Garwood, M., DelaBarre, L.,
            % J. Magn. Reson. 153, 155-177 (2001)
            [~,ind] = min(abs(ti));
            dnu = diff(2*pi*modulation(n).nu/(t{n}(2)-t{n}(1)));
            sweeprate = dnu(ind)/(2*pi*(modulation(n).A(ind))^2);
          end
          
          Exp.Amplitude(n) = sqrt(2*pi*Q_crit*sweeprate)/(2*pi);
          
      end
    end
    
    % ------------------------------------------------------------------- %
    % Calculate pulse IQ function
    % ------------------------------------------------------------------- %
    modulation(n).A = Exp.Amplitude(n)*modulation(n).A;
    totalphase = modulation(n).phi + 2*pi*Exp.CenterFreq(n)*t{n} + Exp.Phase(n);
    y{n} = modulation(n).A.*exp(1i*totalphase);
    
  end
  
  % --------------------------------------------------------------------- %
  % Excitation profile calculation
  % --------------------------------------------------------------------- %
  
  if Opt.ExciteProfile
    
    % Set up offset axis for excitation profile calculation
    if ~isfield(Opt,'Offsets')
      
      % Estimate bandwidth of the pulse
      if ~strcmp(FrequencyModulation,'none')
        
        BW = thisPulse.BW;
      
      else
        
        % Estimate bandwidth from Fourier transform of amplitude modulation
        % function
        if nextpow2(numel(y{n}))<10
          zf = 2^10;
        else
          zf = 4*2^nextpow2(numel(y{n}));
        end
        yft = abs(fftshift(fft(y{n},zf)));
        f = fdaxis(Exp.TimeStep(n),zf);
        intg = cumtrapz(yft);
        [~,indmax] = min(abs(intg-0.5*max(intg)));
        indbw = find(yft(indmax:end)<0.1*max(yft),1);
        BW = 2*(f(indmax+indbw)-f(indmax));
        
      end
      p(n).offsets = linspace(-0.75*BW,0.75*BW,Opt.nOffsets) + Exp.CenterFreq(n);
      
    else
      p(n).offsets = Opt.Offsets;
    end
    nOffsets = numel(p(n).offsets);
    
    % Spin operators
    Sx = sop(1/2,'x');
    Sy = sop(1/2,'y');
    Sz = sop(1/2,'z');
    
    % Equilibrium density matrix
    Density0 = -Sz;
    
    % Detection operator
    switch Opt.Detect
      case 'Sz', Det = -Sz; varname = 'Mz';
      case 'Sy', Det = -Sy; varname = 'My';
      case 'Sx', Det = -Sx; varname = 'Mx';
      case 'all', Det{1} = -Sx; Det{2} = -Sy; Det{3} = -Sz;
      otherwise
        error('Opt.Detect has unknown value. Use ''Sx'', ''Sy'', ''Sz'', or ''all''.');
    end
    
    % Pre-allocate result array
    if ~iscell(Det)
      p(n).(varname) = zeros(1,nOffsets);
    else
      p(n).Mx = zeros(1,nOffsets);
      p(n).My = zeros(1,nOffsets);
      p(n).Mz = zeros(1,nOffsets);
    end
    
    for iOffset = 1:nOffsets
      
      Ham0 = p(n).offsets(iOffset)*Sz;
      
      % Compute pulse propagator
      if min(y{n})==max(y{n}) % rectangular pulses
        
        Ham = real(y{n}(1))*Sx+imag(y{n}(1))*Sy+Ham0;
        tp = Exp.TimeStep(n)*(nPoints-1);
        UPulse = expm(-2i*pi*Ham*tp);
        
      else % general pulses
        
        UPulse = eye(2);
        Isignal = real(y{n});
        Qsignal = imag(y{n});
        for it = 1:nPoints-1
          
          Ham = Isignal(it)*Sx + Qsignal(it)*Sy + Ham0;
          
          % dU = expm(-2i*pi*Ham*Exp.TimeStep(n));
          % Fast matrix exponential for a traceless, antihermitian 2x2 matrix
          M = -2i*pi*Exp.TimeStep(n)*Ham; % M = [a b; -b' -a]
          q = sqrt(M(1,1)^2-abs(M(1,2))^2);
          if abs(q)<1e-10
            dU = eye(2) + M;
          else
            dU = cosh(q)*eye(2) + (sinh(q)/q)*M;
          end
          UPulse = dU*UPulse;
        end
        
      end
      
      % Propagate density matrix
      Density = UPulse*Density0*UPulse';
      
      % Calculate observables
      % (using trace(A*B) = sum(sum(A.*B.')))
      if ~iscell(Det)
        p(n).(varname)(iOffset) = 2*real(sum(sum(Det.*Density.')));
      else
        p(n).Mx(iOffset) = 2*real(sum(sum(Det{1}.*Density.')));
        p(n).My(iOffset) = 2*real(sum(sum(Det{2}.*Density.')));
        p(n).Mz(iOffset) = 2*real(sum(sum(Det{3}.*Density.')));
      end
      
    end
    
    if plotResults
      if n==iPulse(1)
        clf
        cc = winter(numel(iPulse));
        l = 1;
      end
      subplot(2,1,1)
      hold on; box on;
      plot(t{n},real(y{n}),'Color',cc(l,:))
      plot(t{n},imag(y{n}),':','Color',cc(l,:))
      xlabel('t (\mus)')
      ylabel('\nu_1 (MHz)')
      legend('I','Q')
      axis tight
      subplot(2,1,2)
      hold on; box on;
      switch Opt.Detect
        case 'Sz'
          plot(p(n).offsets,real(p(n).Mz),'Color',cc(l,:));
          ylabel('M_z/M_0')
        case 'Sy'
          plot(p(n).offsets,real(p(n).My),'Color',cc(l,:));
          ylabel('M_y/M_0')
        case 'Sx'
          plot(p(n).offsets,real(p(n).Mx),'Color',cc(l,:));
          ylabel('M_x/M_0')
        case 'all'
          plot(p(n).offsets,real(p(n).Mx),...
               p(n).offsets,real(p(n).My),...
               p(n).offsets,real(p(n).Mz));
          ylabel('M_i/M_0')
          legend('x','y','z')
      end
      lgd(l) = strcat(num2str(n),{' '},thisPulse.Type);
      l = l+1;
      legend(lgd,'Location','Best');
      xlabel('frequency offset (MHz)')
      axis tight
      ylim([-1 1])
    end
    
  end
  
end % loop over pulses

% ----------------------------------------------------------------------- %
% Output
% ----------------------------------------------------------------------- %
if numel(iPulse)==1
  t = t{iPulse};
  y = y{iPulse};
  modulation = modulation(iPulse);
end

switch nargout
  case 0 % plotting
  case 2 % [t,y] = pulse(...)
    varargout = {t,y};
  case 3 % [t,y,p] = pulse(...)
    varargout = {t,y,p};
  case 4 % [t,y,p,modulation] = pulse(...)
    varargout = {t,y,p,modulation};
  otherwise
    error('The function pulse() needs 2, 3, or 4 output arguments.')
end
