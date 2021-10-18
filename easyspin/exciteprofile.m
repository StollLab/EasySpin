% exciteprofile      Excitation profile calculation for arbitrary pulses
%
% [offsets,Mag] = exciteprofile(t,IQ)
% [offsets,Mag] = exciteprofile(t,IQ,offsets)
%
% Input:
%   t           = time axis for defined waveform in µs
%   IQ          = in-phase and quadrature part of the pulse function
%   offsets     = axis of frequency offsets in MHz for which to compute the
%                 excitation profile (default: approximately ±2*BW centered
%                 at the pulse center frequency, 201 points)
% Output:
%   offsets     = axis of frequency offsets in MHz
%   Mag         = excitation profile as Mi/M0, i = x,y,z in array
%

function varargout = exciteprofile(varargin)

% ----------------------------------------------------------------------- %
% Input argument parsing
% ----------------------------------------------------------------------- %
switch nargin
  case 0
    help(mfilename);
    return
  case 2 % ... = exciteprofile(t,IQ)
    t = varargin{1};
    IQ = varargin{2};
    Opt = struct;
  case 3 % ... = exciteprofile(t,IQ)
    t = varargin{1};
    IQ = varargin{2};
    Opt.Offsets = varargin{3};
  otherwise
    error('exciteprofile() needs 2 or 3 input arguments.')
end

if numel(t)~=numel(IQ)
  error('The number of points of the time axis (%d) and the number of points in the pulse shape function IQ (%d) do not match.',...
    numel(t),numel(IQ));
end

plotResults = (nargout==0);

% Options
if ~isfield(Opt,'nOffsets') % undocumented
  Opt.nOffsets = 201;
end

% ----------------------------------------------------------------------- %
% Estimate pulse bandwidth and set up offset axis
% --------------------------------------------------------------------- %
if ~isfield(Opt,'Offsets')
  
  % Fourier transform
  if nextpow2(numel(t))<10
    zf = 2^10;
  else
    zf = 4*2^nextpow2(numel(t));
  end
  IQft = abs(fftshift(fft(IQ,zf)));
  f = fdaxis(t(2)-t(1),zf);
  indbw = find(IQft>0.5*max(IQft));
  BW = abs(f(indbw(end))-f(indbw(1))); % ?? real pulses with contributions at +- freq
  CenterFreq = mean([f(indbw(end)) f(indbw(1))]);
  
  offsets = linspace(-BW,BW,Opt.nOffsets) + CenterFreq;
  
else
  
  offsets = Opt.Offsets;
  
end

nPoints = numel(t);
nOffsets = numel(offsets);
 
% --------------------------------------------------------------------- %
% Excitation profile calculation
% --------------------------------------------------------------------- %
 
% Spin operators
Sx = sop(1/2,'x');
Sy = sop(1/2,'y');
Sz = sop(1/2,'z');

% Equilibrium density matrix
Density0 = -Sz;

% Pre-allocate result array
Mag = zeros(3,nOffsets);

Isignal = real(IQ);
Qsignal = imag(IQ);
for iOffset = 1:nOffsets
  
  Ham0 = offsets(iOffset)*Sz;
  
  % Compute pulse propagator
  if min(IQ)==max(IQ) % rectangular pulses
    
    Ham = Isignal(1)*Sx + Qsignal(1)*Sy + Ham0;
    tp = (t(2)-t(1))*(nPoints-1);
    
    % UPulse = expm(-2i*pi*Ham*tp);
    % Fast matrix exponential for a traceless, antihermitian 2x2 matrix
    M = -2i*pi*tp*Ham; % M = [a b; -b' -a]
    q = sqrt(M(1,1)^2-abs(M(1,2))^2);
    if abs(q)<1e-10
      UPulse = eye(2) + M;
    else
      UPulse = cosh(q)*eye(2) + (sinh(q)/q)*M;
    end
    
  else % general pulses
    
    eye2 = eye(2);
    UPulse = eye2;
    for it = 1:nPoints-1
      
      Ham = Isignal(it)*Sx + Qsignal(it)*Sy + Ham0;
      
      % dU = expm(-2i*pi*Ham*Par.TimeStep);
      % Fast matrix exponential for a traceless, antihermitian 2x2 matrix
      M = -2i*pi*(t(2)-t(1))*Ham; % M = [a b; -b' -a]
      q = sqrt(M(1)^2-abs(M(3))^2);
      if abs(q)<1e-10
        dU = eye2 + M;
      else
        dU = cosh(q)*eye2 + (sinh(q)/q)*M;
      end
      UPulse = dU*UPulse;
    end
    
  end
  
  % Propagate density matrix
  Density = UPulse*Density0*UPulse';
  
  % Calculate observables
  % (using trace(A*B) = sum(sum(A.*B.')))
  Mag(1,iOffset) = -2*real(sum(sum(Sx.*Density.')));
  Mag(2,iOffset) = -2*real(sum(sum(Sy.*Density.')));
  Mag(3,iOffset) = -2*real(sum(sum(Sz.*Density.')));
  
end


% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
if plotResults
  
  clf
  S.f = gcf;
  set(S.f,'WindowStyle','normal','Name','exciteprofile output',...
    'numbertitle','off','Units','Normalized',...
    'Position',[0.25,0.30,0.50,0.40],'Color',[1 1 1]*0.8,...
    'Toolbar','figure');
  
  % Set colors
  colI = [0 0 1];
  colQ = [1 0 0];
  colBW = [1 1 1]*0.8;
  colx = [0 0.5 0];
  coly = [1 0 0];
  colz = [0 0 1];
  
  % Set positions
  width = 0.4;
  height = 0.65;
  sep = (1-2*width)/3;
  btm = 0.20;
  boxpos = 0.87;
  
  % IQ plot
  S.label(1) = uicontrol('Style','text','String','Pulse shape:',...
    'FontSize',10,'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Background',[1 1 1]*0.8,'Units','Normalized',...
    'Position',[sep,0.85,width,0.1]);
  S.tick(1) = uicontrol('Style','checkbox',...
    'String','I','Value',1,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[3*sep boxpos 0.1 0.1]);
  S.tick(2) = uicontrol('Style','checkbox',...
    'String','Q','Value',1,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[3.75*sep boxpos 0.1 0.1]);
  S.ha(1) = axes('Units','Normalized','Position',[sep,btm,width,height]);
  hold on; box on;
  line([min(t) max(t)],[0 0],'Color',colBW);
  S.hI = plot(t,real(IQ),'Color',colI);
  S.hQ = plot(t,imag(IQ),'Color',colQ);
  Amax = max(abs(IQ));
  axis([t(1) t(end) -1*Amax*1.1 1*Amax*1.1]);
  xlabel('{\itt} (\mus)')
  ylabel('amplitude (MHz)')
  set(gca,'Layer','top')
  legend([S.hI S.hQ],'I','Q','Location','SouthEast')
    
  % Excitation profile plot
  S.label(2) = uicontrol('Style','text','String','Excitation profiles:',...
    'FontSize',10,'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Background',[1 1 1]*0.8,'Units','Normalized',...
    'Position',[2.5*sep+width,0.85,width,0.1]);
  S.tick(3) = uicontrol('Style','checkbox',...
    'String','x','Value',0,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[2.5*sep+1.5*width boxpos 0.1 0.1]);
  S.tick(4) = uicontrol('Style','checkbox',...
    'String','y','Value',0,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[3.25*sep+1.5*width boxpos 0.1 0.1]);
  S.tick(5) = uicontrol('Style','checkbox',...
    'String','z','Value',1,'Background',[1 1 1]*0.8,...
    'Units','Normalized','Position',[4*sep+1.5*width boxpos 0.1 0.1]);
  S.ha(2) = axes('Units','Normalized','Position',[2.5*sep+width,btm,width,height]);
  hold on; box on;
  line([min(offsets) max(offsets)],[0 0],'Color',colBW);
  S.h(1) = plot(offsets,Mag(1,:),'Color',colx);
  S.h(2) = plot(offsets,Mag(2,:),'Color',coly);
  S.h(3) = plot(offsets,Mag(3,:),'Color',colz);
  set(S.h(1),'Visible','off')
  set(S.h(2),'Visible','off')
  ylabel('{\itM}_i/{\itM}_0')
  legend(S.h,'x','y','z','Location','SouthEast')
  xlabel('frequency (MHz)')
  axis([offsets(1) offsets(end) -1 1])
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
  case 2 % [offsets,Mag] = exciteprofile(...)
    varargout = {offsets,Mag};
  otherwise
    error('Two output arguments are required.')
end

end

% Callback for tick boxes
function showhide(varargin)

S = varargin{3}; % get calling handle structure

for i = 1:numel(S.tick)
  val = get(S.tick(i),'Value');
  if val==1;
    set(S.handles(i),'Visible','on')
  else
    set(S.handles(i),'Visible','off');
  end
end

end