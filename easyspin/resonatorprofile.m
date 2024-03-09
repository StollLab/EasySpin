function varargout = resonatorprofile(nu,nu0,Qu,beta,type)
% resonatorprofile  Resonator profile
%
%  out = resonatorprofile(nu,nu0,Q,beta,type)
%  [nu,out] = resonatorprofile([],nu0,Q,beta,type)
%
%  H = resonatorprofile(nu,nu0,Q,'transferfunction')
%  H = resonatorprofile(nu,nu0,Qu,beta,'transferfunction')
%  Gamma = resonatorprofile(nu,nu0,Qu,beta,'voltagereflection')
%  P = resonatorprofile(nu,nu0,Qu,beta,'powerreflection')
%
%  Calculates the reflection coefficient for a reflection resonator, i.e. the
%  fraction of voltage from an incoming wave that is reflected by the resonator.
%
%  Inputs:
%    nu      frequency range array, GHz
%    nu0     resonator frequency, GHz
%    Qu      unloaded Q-factor of the resonator
%    beta    coupling coefficient
%              beta<1: undercoupled
%              beta=1: critically coupled (matched)
%              beta>1: overcoupled
%    type    type of output to return, 'transferfunction',
%            'voltagereflection' or 'powerreflection'
%
%    If nu is [], then an appropriate frequency range is chosen automatically.
%
%  Outputs:
%    if type = 'transferfunction', the transfer function H is returned
%    if type = 'voltagereflecton', the voltage reflection coefficient Gamma 
%                (ratio of E-field amplitudes of reflected and incoming wave) is returned
%    if type = 'powerreflecton', the power reflection coefficient abs(Gamma)^2 is returned
%
%    If no output is requested, the results are plotted.
%
%  Example:
%    nu = linspace(9.2,9.8,1001);
%    nu0 = 9.5;
%    Qu = 1000;
%    beta = 1;
%    resonatorprofile(nu,nu0,Qu,beta,'voltagereflection');

if nargin==0
  help(mfilename);
  return
end

if nargin==4
  if ischar(beta)
    type = beta;
    beta = [];
  else
    error('Please specify the type of requested output type, ''transferfunction'', ''voltagereflection'', or ''powerreflection''.')
  end
end

if isempty(beta) && contains(type,'reflection')
  error('The coupling coefficient beta is required as a fourth input argument to resonatorprofile().')
end

if isempty(nu)
  dnu = nu0/Qu*10;
  nu = linspace(max(nu0-dnu,0),nu0+dnu,1001);
end

% Transfer function for an ideal RLC series circuit
H = @(nu,nu0,Q) 1./(1+1i*Q*(nu/nu0-nu0./nu));

switch type
  case 'transferfunction'
    if ~isempty(beta)
      Q = Qu/(1+beta);
    else
      Q = Qu;
    end
    out = H(nu,nu0,Q);
  case 'voltagereflection'
    out = (1 - beta*H(nu,nu0,Qu))./(1 + beta*H(nu,nu0,Qu)); % voltage reflection coefficient
  case 'powerreflection'
    out = abs((1 - beta*H(nu,nu0,Qu))./(1 + beta*H(nu,nu0,Qu))).^2; % power reflection coefficient
end

switch nargout
  case 0
    switch type
      case 'transferfunction'
        plot(nu,real(out))       % transfer function
        ylabel('H, real part of the transfer function');
      case 'voltagereflection'
        plot(nu,abs(out));       % reflected voltage
        ylabel('|\Gamma|, norm of voltage reflection coefficient');
      case 'powerreflection'
        plot(nu,out);            % reflected power
        ylabel('P_{reflected}, power reflection coefficient');
    end
    xlabel('frequency (GHz)');
    if ~isempty(beta)
      title(sprintf('Qu = %d, beta = %g',Qu,beta));
    else
      title(sprintf('Q = %d',Q));
    end
    ylim([0 1]);
    line([1 1]*nu0,ylim,'Color','r');
  case 1
    varargout = {out};
  case 2
    varargout = {nu,out};
end
