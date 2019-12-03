% resonator      Simulation of/compensation for the effect of the resonator
%                on a pulse
%
%  [tOut,signalOut] = resonator(tIn,signalIn,mwFreq,nu,TransferFunction,'simulate')
%  [tIn,signalIn] = resonator(tOut,signalOut,mwFreq,nu,TransferFunction,'compensate')
%  [tOut,signalOut] = resonator(tIn,signalIn,mwFreq,nu0,QL,'simulate')
%  [tIn,signalIn] = resonator(tOut,signalOut,mwFreq,nu0,QL,'compensate')
%  ... = resonator(...,Opt)
%
%  If the option 'simulate' is selected, resonator() simulates the effect
%  of the resonator on the input signal signal0.
%  If the option 'compensate' is selected, resonator() adapts the pulse 
%  shape to compensate for the resonator transfer function.
%
%  The resonator transfer function can be given directly by providing
%  a frequency axis and the corresponding transfer function as the second
%  and third input argument. If the input is real, it is assumed to be the
%  magnitude response and the phase response is estimated as described in
%  reference 1. The transfer function is expanded over the required frequency
%  range by exponentially approaching an ideal transfer function with
%  parameters estimated by fit to the provided transfer function.
%
%  Alternatively, the resonator center frequency and the loaded Q-value
%  can be provided as the second and third input argument and the transfer
%  function is calculated based on the ideal transfer function of an RLC
%  series circuit (for details see reference 1).
%
% References:
% 1. Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., Adiabatic and
%    fast passage ultra-wideband inversion in pulsed EPR
%    J. Magn. Reson. 230, 27-39 (2013). (DOI: 10.1016/j.jmr.2013.01.002)
%
%  Input:
%   - t0                   = time axis for the input signal (in µs)
%   - signal0              = input signal vector
%   - mwFreq               = microwave frequency for the input signal
%                            in GHz
%   - nu/nu0               = frequency axis for the resonator transfer 
%                            function (in GHz) or resonator center 
%                            frequency (in GHz)
%   - TransferFunction/QL = resonator transfer function or magnitude
%                            response or loaded Q-value
%   - 'simulate'/'compensate'
%   - Options structure with the following fields:
%        Opt.CutoffFactor     = cutoff factor for truncation of the pulse after
%                               convolution/deconvolution with the resonator
%                               transfer function
%                               (default:1/1000)
%        Opt.TimeStep         = time step in µs (if it is not provided the
%                               ideal time step is estimated based on the 
%                               Nyquist condition with an oversampling factor)
%        Opt.OverSampleFactor = oversampling factor for the determination of the 
%                               time step (default: 10)
%        Opt.N                = multiplication factor used in the determination
%                               of the width of the frequency domain window considered
%                               for the Fourier convolution/deconvolution, the center
%                               frequency and bandwidth of the input signal are 
%                               estimated from the magnitude FT and the width 
%                               is set to Opt.N times the estimated bandwidth
%        Opt.Window           = type of apodization window used on the frequency
%                               domain output (over the width determined by Opt.N)
%                               (see apowin() for available options) (default: 'gau')
%        Opt.alpha            = alpha parameter the apodization function defined in
%                               Opt.Window (see apowin() for details) (default: 0.6)
%
%  Output:
%   - t         = time axis for the output signal (in µs)
%   - signal    = signal modified by the resonator transfer function or
%                 compensated for the resonator transfer function
%

function [t,signal] = resonator(t0,signal0,mwFreq,varargin)

% Input argument check
% -----------------------------------------------------------------------
if nargin==0
  help(mfilename);
  return
end

if nargin<6
  error('resonator() requires at least 6 input arguments.')
end

dim = size(signal0);
if min(dim)>1
  error('Input must be a row or column vector.');
end
if numel(dim)>2
  error('Input must be a row or column vector.');
end
if numel(t0)~=numel(signal0)
  error('The time axis and input signal need to be vectors of the same length.')
end

% Check input arguments and get transfer function
if numel(varargin{1})==1 && numel(varargin{2})==1
  [f,H] = transferfunction('ideal',varargin{1},varargin{2});
else
  if numel(varargin{1})~=numel(varargin{2})
    error('The frequency axis and resonator transfer function need to have the same length.')
  end
  [f,H] = transferfunction('experimental',varargin{1},varargin{2});
end

option = varargin{3};

% Default options
if nargin==7
  Opt = varargin{4};
else
  Opt = struct();
end
if ~isfield(Opt,'CutoffFactor') || isempty(Opt.CutoffFactor)
  Opt.CutoffFactor = 1/1000;
end
if ~isfield(Opt,'OverSampleFactor') || isempty(Opt.OverSampleFactor)
  Opt.OverSampleFactor = 10;
end
if ~isfield(Opt,'TimeStep') || isempty(Opt.TimeStep)
  estimateTimeStep = true;
else
  estimateTimeStep = false;
end
if ~isfield(Opt,'N') || isempty(Opt.N)
  Opt.N = 20;
end
if ~isfield(Opt,'Window') || isempty(Opt.Window)
  Opt.Window = 'gau';
end
if (~isfield(Opt,'alpha') || isempty(Opt.alpha))
  Opt.alpha = 0.6;
end

% ---------------------------------------------------------------------
% Simulation or compensation for the resonator
% ---------------------------------------------------------------------
  
% Pulse Fourier transform
dt = (t0(2)-t0(1))/4;
signal_ = interp1(t0,signal0,0:dt:t0(end),'spline');  
N = round((numel(f)-numel(signal_))/2);
signal_ = [zeros(1,N-1) signal_ zeros(1,N)];

FT = ifftshift(fft(fftshift(signal_)));
f_ = fdaxis(dt,numel(FT));

% Extract pulse frequency response
if isreal(signal0) && mwFreq==0
  [~,ind0] = min(abs(f_));
  intg = cumtrapz(abs(FT(ind0:end)));
  [~,indmax] = min(abs(intg-0.5*max(intg)));
  indmax = ind0+indmax;
  indbw = find(abs(FT(indmax:end))>0.01*max(abs(FT)),1,'last');
else
  intg = cumtrapz(abs(FT));
  [~,indmax] = min(abs(intg-0.5*max(intg)));
  indbw = find(abs(FT(indmax:end))>0.01*max(abs(FT)),1,'last');
end
startind = indmax-Opt.N*indbw;
endind = indmax+Opt.N*indbw;
delta = max(max(endind-numel(f_),1-startind),0);
indpulse = startind+delta:endind-delta;
f_pulse = f_(indpulse);
FT_pulse = FT(indpulse);

% Interpolation of the transfer function onto the same axis
H_ = interp1((f-mwFreq*1e3),H,f_pulse,'pchip');

switch option
  case 'simulate'
  % Fourier convolution
  FTc_pulse = FT_pulse.*H_;
  case 'compensate'
  % Fourier deconvolution
  FTc_pulse = FT_pulse./H_;
end

% Windowing
if any(strcmp(Opt.Window,{'gau','exp','kai'}))
  windowfunction = apowin(Opt.Window,numel(FTc_pulse),Opt.alpha);
else
  windowfunction = apowin(Opt.Window,numel(FTc_pulse));
end
FTc_pulse = FTc_pulse.*windowfunction.';

% Inverse Fourier transform
FTc = zeros(size(FT));
FTc(indpulse) = FTc_pulse;
t_ = (-0.5*numel(FTc):1:0.5*numel(FTc)-1)*dt;
signal_ = fftshift(ifft(ifftshift(FTc)));

% Extract pulse
switch option
  case 'simulate'
    % Start of the pulse determined by input signal
    endind = find(abs(signal_)>(Opt.CutoffFactor*max(abs(signal_))),1,'last');
    [~,startind] = min(abs(t_+t0(end)/2));
    t_ = t_(startind:endind)-t_(startind);
    signal_ = signal_(startind:endind);
  case 'compensate'
    ind = find(abs(signal_)>(Opt.CutoffFactor*max(abs(signal_))));
    ind = ind(1):1:ind(end);
    t_ = t_(ind)-t_(ind(1));
    signal_ = signal_(ind);
end
    
% Update time step
if estimateTimeStep
  
  % Get maximum frequency
  [maxvalue,indmax] = max(abs(FTc));
  indbw = find(abs(FTc(indmax:end))>0.05*maxvalue,1,'last');
  maxFreq = 2*(f_(indmax+indbw-1)-f_(indmax));
  
  % Define new time step
  Nyquist_dt = 1/(2*maxFreq);
  newTimeStep = Nyquist_dt/Opt.OverSampleFactor;
  inputTimeStep = t0(2)-t0(1);
  if newTimeStep<inputTimeStep
    Opt.TimeStep = newTimeStep;
  else
    Opt.TimeStep = inputTimeStep;
  end
end
  
t = 0:Opt.TimeStep:t_(end);
signal = interp1(t_,signal_,t,'spline');

if isreal(signal0) && mwFreq==0
  signal = 2*real(signal);
end
  
end

% Transfer function calculation
% -----------------------------------------------------------------------
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
    N = 2^18; % number of frequency points
    df = fmax/N; % GHz
    f = (0:1:(N-1))*df; % GHz
    
    H = Hideal(f,f0,QL,1);
    H(1) = H(2);
    
  case 'experimental'
    
    % Extrapolation of the experimental transfer function over
    % the full frequency range
    FrequencyAxis = varargin{1}; % frequency axis in GHz
    FrequencyResponse = varargin{2}; % experimental resonator profile
    
    % Frequency axis
    fmax = 2*ceil(max(FrequencyAxis)); % Nyquist frequency, GHz
    N = 2^14; % number of frequency points
    df = fmax/N; % GHz
    f = (0:1:(N-1))*df; % GHz
    
    % Interpolation
    ind = find((f>FrequencyAxis(1)) & (f<FrequencyAxis(end)));
    FrequencyResponse_ = interp1(FrequencyAxis,FrequencyResponse,f(ind),'spline');
    
    % Separate real and imaginary contributions, if available
    if ~isreal(FrequencyResponse_)
      FrequencyResponse_r = real(FrequencyResponse_);
      FrequencyResponse_i = imag(FrequencyResponse_);
      FrequencyResponse_ = abs(FrequencyResponse_);
    end
    
    % Estimate center frequency and loaded Q
    [v1max,maxind] = max(FrequencyResponse_);
    f0 = f(ind(1)+maxind); % GHz
    
    v1_3dB = (3/4)*v1max;
    ind_3dB = [find(FrequencyResponse_>v1_3dB,1,'first') find(FrequencyResponse_>v1_3dB,1,'last')];
    f_3dB = f(ind(1)-1+ind_3dB);
    BW_3dB = diff(f_3dB);
    QL = f0/BW_3dB;
    
    % Fit f0 and QL for best overlap
    fitfunc = @(x) sqrt(sum((abs(Hideal(f(ind),x(1),x(2),v1max))-FrequencyResponse_).^2)/numel(FrequencyResponse_));
    x0 = [f0 QL];
    x = fminsearch(fitfunc,x0);
    f0 = x(1);
    QL = x(2);
    
    Hid = Hideal(f,f0,QL,v1max);
    Hid(1) = Hid(2);
    
    % Extrapolate the experimental transfer function by exponentially approaching |H0(f)|
    while FrequencyResponse_(1)-abs(Hid(ind(1)))<-0.05*v1max
      FrequencyResponse_ = FrequencyResponse_(2:end);
      ind = ind(2:end);
    end
    while FrequencyResponse_(end)-abs(Hid(ind(end)))<-0.05*v1max
      FrequencyResponse_ = FrequencyResponse_(1:end-1);
      ind = ind(1:end-1);
    end
    H0(1:ind(1)) = (FrequencyResponse_(1)-abs(Hid(ind(1))))*exp(1/log(2)*(f(1:ind(1))-f(ind(1)))) + abs(Hid(1:ind(1)));
    H0(ind) = FrequencyResponse_;
    H0(ind(end):numel(f)) = (FrequencyResponse_(end)-abs(Hid(ind(end))))*exp(-1/log(2)*(f(ind(end):end)-f(ind(end)))) + abs(Hid(ind(end):end));
    
    % Phase response
    betaid = atan(imag(Hid)./real(Hid));
    if exist('FrequencyResponse_i','var')
      betaexp(ind) = atan(FrequencyResponse_i./FrequencyResponse_r);
    else
      % Estimate phase response
      betaexp = -imag(hilberttrans(log(abs(H0))));
    end
    % Combine phase response of the ideal function and the provided measured data
    fitfunc = @(x) sqrt(sum(((betaexp(ind)-x)-betaid(ind)).^2)/(ind(end)-ind(1)));
    x = fminsearch(fitfunc,0);
    beta(ind) = betaexp(ind)-x;
    beta(1:ind(1)) = ((betaexp(ind(1))-x)-betaid(ind(1)))*exp(1/log(2)*(f(1:ind(1))-f(ind(1)))) + betaid(1:ind(1));
    beta(ind(end):numel(f)) = ((betaexp(ind(end))-x)-betaid(ind(end)))*exp(-1/log(2)*(f(ind(end):end)-f(ind(end)))) + betaid(ind(end):end);
    
    % Transfer function
    H = abs(H0).*exp(1i*beta);    
    H = real(H) - 1i*imag(hilberttrans(real(H)));
    
end

H = H/max(real(H));
f = f*1e3; % GHz to MHz

end
