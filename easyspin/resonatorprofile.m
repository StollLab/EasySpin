function varargout = resonatorprofile(nu,nu0,Qu,beta,type)
% resonatorprofile  Various quantities characterizing a microwave resonator
%
%  resonatorprofile(nu,nu0,Qu,beta,type)
%  result = resonatorprofile(nu,nu0,Qu,beta,type)
%  [nu,result] = resonatorprofile([],nu0,Qu,beta,type)
%
%  H = resonatorprofile(nu,nu0,Qu,beta,'transferfunction')
%  Gamma = resonatorprofile(nu,nu0,Qu,beta,'voltagereflection')
%  Gamma2 = resonatorprofile(nu,nu0,Qu,beta,'powerreflection')
%  Gamma2 = resonatorprofile(nu,nu0,Qu,beta,'powertransmission')
%  S11 = resonatorprofile(nu,nu0,Qu,beta,'S11')
%  RL = resonatorprofile(nu,nu0,Qu,beta,'returnloss')
%  VSWR = resonatorprofile(nu,nu0,Qu,beta,'VSWR')
%
%  Calculates characteristics of a microwave reflection resonator such as the
%  reflected voltage or power and the transmitted power.
%
%  Inputs:
%    nu      frequency range array, GHz
%    nu0     resonator frequency, GHz
%    Qu      unloaded Q-factor of the resonator
%    beta    coupling coefficient between resonator and transmission line
%              0<beta<1: undercoupled
%              beta=1: critically coupled (matched)
%              beta>1: overcoupled
%    type    type of output to return, possible values:
%              'transferfunction'
%              'voltagereflection'
%              'powerreflection',
%              'powertransmission',
%              'S11',
%              'returnloss'
%              'VSWR'
%
%    If nu is [], then an appropriate frequency range is chosen automatically.
%
%  Outputs:
%    type=='transferfunction':  transfer function H (= admittance)
%    type=='voltagereflection': the voltage reflection coefficient Gamma 
%                (ratio of E-field amplitudes of reflected and incoming wave)
%    type=='powerreflection': the power reflection coefficient abs(Gamma)^2
%    type=='S11': the S11 scattering parameter
%    type=='VSWR': the voltage standing wave ratio
%    type=='returnloss': the return loss, in dB
%
%    If no output is requested, the results are plotted.
%
%  Example:
%    nu = linspace(9.2,9.8,1001);  % GHz
%    nu0 = 9.5;  % GHz
%    Qu = 1000;
%    beta = 1;
%    resonatorprofile(nu,nu0,Qu,beta,'powerreflection');

% References:
% - Doll et al, J.Magn.Reson. 2013, https://doi.org/10.1016/j.jmr.2013.01.002
% - Krymov et al, J.Magn.Reson. 2003, https://doi.org/10.1016/S1090-7807(03)00109-5
% - Pozar, Microwave Engineering, 4th ed. (2012), Ch.6

if nargin==0
  help(mfilename);
  return
end

% Input parameter checks
%-------------------------------------------------------------------------------
if nargin~=5
  error('Five inputs (nu, nu0, Qu, beta, type) are required.');
end

if ~isempty(nu)
  if ~isnumeric(nu) || ~isreal(nu) || any(nu<0)
    error('The frequency axis (1st input) must be a 1D array of positive numbers.');
  end
end

if ~isnumeric(nu0) || ~isscalar(nu0) || ~isreal(nu0) || nu0<=0
  error('The resonator frequency (2nd input) must be a positive number.');
end

if ~isnumeric(Qu) || ~isscalar(Qu) || ~isreal(Qu) || Qu<=0
  error('The unloaded Q factor (3th input) must be a positive number.');
end

if ~isnumeric(beta) || ~isscalar(beta) || ~isreal(beta) || beta<=0
  error(sprintf(['The coupling coefficient beta (4th input) must be a positive number.\n'...
    '(1 for critical coupling, 0<beta<1 for undercoupling, >1 for overcoupling.']));
end

validStrings = {'powerreflection','voltagereflection','transferfunction','powertransmission','S11','returnloss','VSWR'};
type = validatestring(type,validStrings,mfilename,'type (5th input)');

% Calculations
%-------------------------------------------------------------------------------
% Calculate loaded Q factor
QL = Qu/(1+beta);  % loaded Q

% Set up frequency axis if not provided
if isempty(nu)
  bw = nu0/QL;  % bandwidth
  dnu = bw*10;
  nu = linspace(max(nu0-dnu,0),nu0+dnu,1001);
end

% Voltage reflection coefficient (=S11) for loaded resonator
xi = nu/nu0 - nu0./nu;  % == approx 2*(nu-nu0)/nu0
Gamma = (1-beta+1i*Qu*xi)./(1+beta+1i*Qu*xi);  % Krymov Eqs.(6) and (7)

% Calculate desired quantity
switch type
  case 'transferfunction'
    H = 1./(1+1i*QL*xi);  % Doll 2013 Eq.(5)
    result = H;

  case 'voltagereflection'
    result = Gamma;

  case 'powerreflection'
    Gamma2 = 1 - 4*beta./(Qu^2*xi.^2+(1+beta).^2);
    %Gamma2 = abs(Gamma).^2;  % equivalent
    result = Gamma2;

  case 'powertransmission'
    Ptrans = 4*beta./(Qu^2*xi.^2+(1+beta).^2);  % Lorentzian
    %Ptrans = 1 - abs(Gamma).^2;  % equivalent
    result = Ptrans;

  case 'S11'
    S11 = Gamma;
    result = S11;

  case 'returnloss'
    RL = -20*log10(abs(Gamma));  % Pozar, p. 58, eq. (2.38)
    result = RL;

  case 'VSWR'
    VSWR = (1+abs(Gamma))./(1-abs(Gamma));  % Pozar, p. 58, eq. (2.41)
    result = VSWR;

  otherwise
    error(sprintf(['Type (5th input) must be one of the following:\n  ''transferfunction'', ''voltagereflection'','...
      '''powerreflection'', ''powertransmission'', ''S11'', ''VSWR'', ''returnloss''.']));
end

% Plotting, output
%-------------------------------------------------------------------------------
ylims = [];
switch nargout
  case 0
    switch type
      case 'transferfunction'
        plot(nu,real(result),nu,imag(result))
        ylabel('Re({\itH}), Im({\itH}), transfer function');
        legend('Re({\itH})', 'Im({\itH})','AutoUpdate','off');
        legend boxoff
        ylims = [-1 1];

      case 'voltagereflection'
        plot(nu,real(result),nu,imag(result),nu,abs(result));
        ylabel('Re({\it\Gamma}), Im({\it\Gamma}), |{\it\Gamma}|, voltage reflection coefficient');
        ylims = [-1 1];
        legend('Re({\it\Gamma})','Im({\it\Gamma})','|{\it\Gamma}|','AutoUpdate','off');
        legend boxoff
        yline(0)

      case 'powerreflection'
        plot(nu,result);
        ylabel('|{\it\Gamma}|^2, power reflection coefficient');
        ylims = [0 1];

      case 'powertransmission'
        plot(nu,result);
        ylabel('{\itP}_{trans}, power transmission coefficient');
        ylims = [0 1];

      case 'S11'
        plot(nu,real(S11),nu,imag(S11),nu,abs(S11));
        ylabel('Re({\itS}_{11}),  Im({\itS}_{11}),  |{\itS}_{11}|');
        ylims = [-1 1];
        legend('Re({\itS}_{11})','Im({\itS}_{11})','|{\itS}_{11}|','AutoUpdate','off');
        legend boxoff
        yline(0)

      case 'returnloss'
        plot(nu,RL)
        ylabel('return loss, dB');

      case 'VSWR'
        semilogy(nu,VSWR);
        ylabel('VSWR, voltage standing wave ratio');
        ylims = [1 max(VSWR)];
        grid on
    end
    axis tight
    if ~isempty(ylims)
      ylim(ylims)
    end
    xlabel('frequency (GHz)');
    xline(nu0,'--','');
    title(sprintf('{\\nu}_0 = %g GHz, {\\it{Q}}_u = %g, \\beta = %g, {\\it{Q}}_L = %g',nu0,Qu,beta,Qu/(1+beta)));

  case 1
    varargout = {result};

  case 2
    varargout = {nu,result};
end

end
