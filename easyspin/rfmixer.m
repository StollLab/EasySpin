% rfmixer  Digital up- or downconversion
%
%   [tOut,signalOut] = rfmixer(tIn,signalIn,mwFreq,type) 
%   [tOut,signalOut] = rfmixer(tIn,signalIn,mwFreq,type,Opt) 
%
%   Mixes the input signal with the LO frequency mwFreq. Depending on the 
%   selected mixer type, the function acts as a double-sideband (DSB) mixer,
%   single-sideband (SSB) mixer or performs IQ modulation, demodulation or 
%   frequency shifting.
%
%   Resampling of the input signal before up/downconversion can be
%   performed if the new sampling time step is given in Opt.dt.
%
%   If no outputs are requested, the results are plotted.
%
%   Input:
%   - tIn:        time axis, in ?s
%   - signalIn:   input signal vector, in-phase part only for DSB,
%                 SSB or IQ demodulation mode, in-phase and quadrature
%                 part for IQ mixer and IQ frequency shift operation
%   - mwFreq:     LO frequency, in GHz
%                 (for IQ frequency shifts, the direction of the shift
%                 has to be included)
%   - type :      mode of operation, the options are:
%                 'DSB' = double sideband mixer
%                 'USB' = single sideband mixer, upper sideband selected
%                 'LSB' = single sideband mixer, lower sideband selected
%                 'IQmod' = IQ modulator
%                 'IQdemod' = IQ demodulator
%                 'IQshift' = IQ frequency shifter (up- or downconversion)
%   - Options:
%     Opt.dt =     time step for the output signal, in ?s.
%                  If no time step for resampling is given and the input
%                  time axis is too large, a new time step is computed as
%                  1/(2*Opt.OverSampleFactor*maxfreq).
%     Opt.OverSampleFactor     = factor for determining new time step in
%                                signal resampling (default = 1.25).
%     Opt.InterpolationMethod  = interpolation method for signal resampling
%                                (default = 'spline').
%     Opt.BandwidthThreshold   = threshold for input bandwidth determination
%                                (default = 0.1).
%     Opt.HilbertThreshold     = threshold for amplitude and phase FT overlap,
%                                the Hilbert transform corresponds to the
%                                quadrature signal only if amplitude and phase
%                                FT do not overlap (Ref. 1). If the maximum of
%                                the product of the normalized amplitude and
%                                phase FT (extracted from real input signal
%                                by Hilbert transform) is larger than the
%                                specified threshold, recovery of the I and
%                                Q data from the real input signal is not
%                                possible (default = 0.05).
%     Opt.NoiseCutoffThreshold = threshold for selection of part of the input
%                                signal used to compute amplitude and cos(phase)
%                                FT overlap (the evaluation of the validity
%                                of the Hilbert transform is affected by
%                                the noise level in the signal and baseline).
%
%   Output:
%   - tOut:          time axis for the output signal, in ?s.
%   - signalOut:     output signal vector, in-phase component only for DSB,
%                    SSB or IQmod, in-phase and quadrature component for
%                    IQdemod and IQshift
%
% References:
% 1. Bedrosian's product theorem is explained in:
%    Boashash, B., Estimating and interpreting the instantaneous frequency 
%    of a signal. I. Fundamentals
%    Proc. IEEE 80, 520 (1992). (DOI: 10.1109/5.135376)
%

function varargout = rfmixer(varargin)

if nargin==0
  help(mfilename);
  return
end

% Argument parsing
%-----------------------------------------------------------
switch nargin
  case {1,2,3}
    error('rfmixer() requires at least four input arguments.');
  case 4
    t = varargin{1};
    signal = varargin{2};
    mwFreq = varargin{3};
    type = varargin{4};
    Opt = struct;
  case 5
    t = varargin{1};
    signal = varargin{2};
    mwFreq = varargin{3};
    type = varargin{4};
    Opt = varargin{5};
  otherwise
    error('rfmixer() cannot take %d input arguments.',nargin);
end

% Input checks
%-----------------------------------------------------------

% Check input signal is consistent with mode
if any(strcmpi(type,{'DSB','USB','LSB','IQdemod'}));
  if ~isreal(signal)
    error('A real input signal is required for the selected mixer type.')
  end
elseif any(strcmpi(type,{'IQmod','IQshift'}));
  if isreal(signal)
%     warning(['Both input and quadrature components of the input signal are ',...
%            'required for the selected mixer type. The quadrature component ',...
%            'is assumed to be zero.'])
    signal = signal + 1i*ones(size(signal))*1e-300;
  end
else
  error('The selected mixer type is not available.')
end
realInput = isreal(signal);

dim = size(signal);
if min(dim)>1 || numel(dim)>2
  error('Input must be a row or column vector.');
end
if numel(t)~=numel(signal)
  error('The signal and its time axis must have equal lengths.');
end
if size(t)~=size(signal)
  t = reshape(t,size(signal));
end

% Default options
%-----------------------------------------------------------
% Threshold for input bandwidth determination
if ~isfield(Opt,'BandwidthThreshold'), Opt.BandwidthThreshold = 0.1; end 
% Interpolation method for resampling
if ~isfield(Opt,'InterpolationMethod'), Opt.InterpolationMethod = 'spline'; end 
% Factor for oversampling in the determination of the new time step in resampling
if ~isfield(Opt,'OverSampleFactor'), Opt.OverSampleFactor = 1.25; end
% Threshold for part of the input signal used to compute amplitude and
% cos(phase) FT overlap (the evaluation of the validity of the Hilbert
% transform is affected by the noise level in the signal and baseline)
if ~isfield(Opt,'NoiseCutoffThreshold'), Opt.NoiseCutoffThreshold = 0.1; end
% Threshold for amplitude and cos(phase) FT overlap
if ~isfield(Opt,'HilbertThreshold'), Opt.HilbertThreshold = 0.05; end
% Hilbert transform corresponds to quadrature signal only if amplitude
% and cos(phase) FT do not overlap.
% (based on Bedrosian's product theorem; see
%  Boashash, Proc. IEEE 80, 1992, 520, http://dx.doi.org/10.1109/5.135376)

% Signal resampling
%-----------------------------------------------------------
% Determine max frequency of input and output signals
f = fdaxis(t); % GHz or MHz
signalFT = fftshift(fft(signal));
inputband = f(abs(signalFT)>Opt.BandwidthThreshold*max(abs(signalFT)));
maxFreqIn = max(inputband);

if abs(mwFreq*1e3) > maxFreqIn % for upconversion
  maxFreqOut = abs(mwFreq*1e3) + maxFreqIn;
else % for downconversion
  maxFreqOut = max([abs(mwFreq*1e3) maxFreqIn]);
end
nyqdt = 1/(2*maxFreqOut); % ?s, dt for Nyquist criterion

tIn = t;
dtIn = tIn(2) - tIn(1); % input signal time step, in ?s
if isfield(Opt,'dt')
  if Opt.dt>nyqdt
    warning('Maximum frequency exceeds Nyquist frequency for specified resampling time step.')
  end
else
  if dtIn > nyqdt
    Opt.dt = 1/(2*Opt.OverSampleFactor*maxFreqOut); % default resampling time step, in ?s
  else
    Opt.dt = dtIn;
  end
end
tOut = tIn(1):Opt.dt:tIn(end);
signal_rs = interp1(tIn,signal,tOut,Opt.InterpolationMethod);

% Define I and Q data for real input signals
%-----------------------------------------------------------
if realInput
  if any(strcmpi(type,{'USB','LSB','IQdemod'}))
    
    % Use Hilbert transform to construct quadrature signal
    signal_rs = signal_rs + 1i*imag(hilberttrans(signal_rs));
    
    % Neglect part of the signal amplitude below Opt.NoiseCutoffThreshold
    ind = find(abs(signal_rs)>Opt.NoiseCutoffThreshold*max(abs(signal_rs)));
    signal_nonzero = signal_rs(ind(1):ind(end));

    % Check that input satisfies conditions for Hilbert transform to
    % correspond to the quadrature signal (based on Bedrosian's theorem)
    amplitude_FT = fft(abs(signal_nonzero));
    cosphase_FT = fft(cos(angle(signal_nonzero)));
    amplitude_FT = amplitude_FT/max(amplitude_FT);
    cosphase_FT = cosphase_FT/max(cosphase_FT);
    
    if max(abs(amplitude_FT.*cosphase_FT))>Opt.HilbertThreshold
      error(['SSB conversion and digital downconversion require the quadrature signal',...
             ' of the real input signal, which is calculated internally using the Hilbert',...
             ' transform. In this case, it cannot be recovered. (Use the quadrature signal', ...
             ' as input or increase Opt.HilbertThreshold to still perform the calculation.)']);
    end
    
  end
end

% Up/downconversion
%-----------------------------------------------------------
if strcmpi(type,'USB')
  LO = abs(mwFreq);
elseif strcmpi(type,'LSB') || strcmpi(type,'IQdemod')
  LO = -abs(mwFreq);
else
  LO = mwFreq;
end
signalOut = signal_rs.*exp(2i*pi*LO*1e3*tOut);

if any(strcmpi(type,{'DSB','USB','LSB','IQmod'}))
  signalOut = real(signalOut);
end

% Output, plotting
%-----------------------------------------------------------
switch nargout
  case 0
    subplot(2,1,1);
    if ~realInput
      plot(tOut,real(signal_rs),tOut,imag(signal_rs));
      legend('I','Q');
      legend boxoff
    else
      plot(tOut,real(signal_rs));
    end
    title('Input signal');
    xlabel('{\itt} (\mus)');
    ylim([-1.1 1.1]*max(abs(signal_rs)))
    subplot(2,1,2);
    if ~isreal(signalOut)
      plot(tOut,real(signalOut),tOut,imag(signalOut));
      legend('I','Q');
      legend boxoff
    else
      plot(tOut,real(signalOut));
    end
    title('Output signal');
    xlabel('{\itt} (\mus)');
    ylim(1.1*[-1 1]*max(abs(signalOut)));
    
  case 1
    error('rfmixer() needs to be called with at least two output arguments.\n')
    
  case 2
    varargout{1} = tOut;
    varargout{2} = signalOut;
   
  otherwise
    error('rfmixer() cannot return %d output arguments.',nargout);
end

return
