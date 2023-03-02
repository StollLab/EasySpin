% rotview    Display rotation between two frame
%
%    rotview()
%    rotview([a b c])
%    rotview(a,b,c)
%    rotview(R)
%
% Plots two frames and allows adjustment of their relative orientation.
% Euler angles and angles between axis pairs are shown. If inputs are
% given, they are used as the starting Euler angles.
%
% Input:
%    a, b, c:    Euler angles alpha, beta, and gamma (in radians)
%    R:          rotation matrix, 3x3
%

function rotview(varargin)

% Parameters
%-------------------------------------------------------------------------------
switch nargin
  case 1
    R = varargin{1};
    if isequal(size(R),[3 3])
      [alpha,beta,gamma] = eulang(R);
    elseif numel(varargin{1})==3
      angles_in = varargin{1};
      alpha = angles_in(1);
      beta = angles_in(2);
      gamma = angles_in(3);
    else
      error('For one input argument, it must be a 3x3 rotation matrix or a 3-element vector of Euler angles.');
    end
  case 3
    alpha = varargin{1};
    beta = varargin{2};
    gamma = varargin{3};
  otherwise
    alpha = 30*pi/180;
    beta = 20*pi/180;
    gamma = 40*pi/180;
end

% convert to degrees
alpha = alpha*180/pi;
beta = beta*180/pi;
gamma = gamma*180/pi;

showTensor = true;

% Initialize quantities derived from the Euler angles
xA = 0; yA = 0; zA = 0; xB = 0; yB = 0; zB = 0;
x1 = 0; x2 = 0; y1 = 0; zzplane = 0; xyAplane = 0; xyBplane = 0;
R = [];
minangle = [-360 -360 -360]; % [0 0 0];
maxangle =  [360 360 360]; % [360 180 360];

if alpha<minangle(1)
  alpha = minangle(1);
end
if alpha>maxangle(1)
  alpha = maxangle(1);
end
if beta<minangle(2)
  beta = minangle(2);
end
if beta>maxangle(2)
  beta = maxangle(2);
end
if gamma<minangle(3)
  gamma = minangle(3);
end
if gamma>maxangle(3)
  gamma = maxangle(3);
end

calculateframes();


% UI setup
%------------------------------------------------------------------
sz = [1000 600]; % figure size
screensize = get(0,'ScreenSize');
scalefact = min(0.9*(screensize(3:4)/sz));
if scalefact>1, scalefact = 1; end
sz = sz*scalefact;
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
hFig = figure(443123);
set(hFig,'Position',[xpos, ypos, sz(1), sz(2)],'Units','pixels');
set(hFig,'Resize','off');
set(hFig,'Name','Rotation viewer',...
  'NumberTitle','off','WindowStyle','normal','MenuBar','none','Color','w');
clf(hFig)

x0 = 10;
y0 = 10;
uicontrol('Style','slider','Tag','alphaslider',...
  'Position',[x0+50 y0+40 120 15],...
  'Min',minangle(1),'Max',maxangle(1),'Value',alpha,...
  'SliderStep',[1 10]/(maxangle(1)-minangle(1)),...
  'callback',@(src,evt) rotview_updateangle(1,src));
uicontrol('Style','text','Position',[x0,y0+40,50,15],...
  'String','alpha (deg)','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','edit','Tag','alphatext',...
  'Position',[x0+180,y0+40,50,15],...
  'String','','HorizontalA','left',...
  'callback',@(src,evt) rotview_updateangle(1,src));

uicontrol('Style','slider','Tag','betaslider',...
  'Position',[x0+50 y0+20 120 15],...
  'Min',minangle(2),'Max',maxangle(2),'Value',beta,...
  'SliderStep',[1 10]/(maxangle(2)-minangle(2)),...
  'callback',@(src,evt) rotview_updateangle(2,src));
uicontrol('Style','text','Position',[x0,y0+20,50,15],...
  'String','beta (deg)','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','edit','Tag','betatext',...
  'Position',[x0+180,y0+20,50,15],...
  'String','','HorizontalA','left',...
  'callback',@(src,evt) rotview_updateangle(2,src));

uicontrol('Style','slider','Tag','gammaslider',...
  'Position',[x0+50 y0 120 15],...
  'Min',minangle(3),'Max',maxangle(3),'Value',gamma,...
  'SliderStep',[1 10]/(maxangle(3)-minangle(3)),...
  'callback',@(src,evt) rotview_updateangle(3,src));
uicontrol('Style','text','Position',[x0,y0,50,15],...
  'String','gamma (deg)','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','edit','Tag','gammatext',...
  'Position',[x0+180,y0,50,15],...
  'String','','HorizontalA','left',...
  'callback',@(src,evt) rotview_updateangle(3,src));

labelStrings = {'general (A->B)','Sys.gFrame (M->g)','Sys.AFrame (M->A)',...
  'Sys.DFrame (M->D)','Sys.QFrame (M->Q)','Sys.eeFrame (M->ee)',...
  'Sys.DiffFrame (M->Diff)','Exp.MolFrame (C->M)','Exp.SampleFrame (L->C)'};
labelID = 1;
uicontrol('Style','text','Position',[x0,y0+65,50,15],...
  'String','labeling','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','popupmenu','Tag','labelpopupmenu',...
  'Position',[x0+50 y0+70 180 15],...
  'String',labelStrings,...
  'callback',@rotview_labelpopupmenu);

uicontrol('Style','text','Position',[x0,y0+90,50,15],...
  'String','ellipsoid','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','checkbox','Tag','tensorcheck',...
  'Position',[x0+50 y0+90 120 15],...
  'Value',showTensor,...
  'Background','w',...
  'callback',@rotview_toggletensor);

% Rotation matrix display
%-------------------------------------------------------------------------------
xR0 = x0+50;
yR0 = y0+110;
dx = 50;
dy = 20;
xyzStr = {'x','y','z'};
for r = 1:3
  hRow(r) = uicontrol('Style','text','String',[xyzStr{r} '_B'],...
    'Position',[xR0-dx,yR0+(3-r)*dy,dx,dy]);
end
set(hRow,'Background','w','ForegroundColor','k','HorizontalAlignment','right');

for c = 1:3
  hCol(c) = uicontrol('Style','text','String',[xyzStr{c} '_A'],...
    'Position',[xR0+(c-1)*dx,yR0+(3-0)*dy,dx,dy]);
end
set(hCol,'Background','w','ForegroundColor','k','HorizontalAlignment','right');

for r = 1:3
  for c = 1:3
    hMatrixElement(r,c) = uicontrol(...
      'Style','text',...
      'HorizontalAl','right',...
      'Position',[xR0+(c-1)*dx,yR0+(3-r)*dy,dx,dy]);
  end
end
set(hMatrixElement,'Background','w','ForegroundColor','k');

% Plotting
%------------------------------------------------------------------

% Axes
hline(1) = line('Color','r');
hline(2) = line('Color',[0 0.6 0]);
hline(3) = line('Color','b');
set(hline(1:3),'Linewidth',2,'LineStyle','--');
hline(4) = line('Color','r');
hline(5) = line('Color',[0 0.6 0]);
hline(6) = line('Color','b');
set(hline(4:6),'Linewidth',2,'LineStyle','-');
hline(7) = line('Color','k');
hline(8) = line('Color','k');
hline(9) = line('Color','k');
hline(10) = line('Color','k');

% Text labels
for k = 1:6
  t(k) = text();
end
setaxislabels();

% Tensor ellipsoid
[xTensor,yTensor,zTensor] = sphere(100);
cTensor = [0.2 0.3 0.4];
ht = surface(cTensor(1)*xTensor,cTensor(2)*yTensor,cTensor(3)*zTensor,'Tag','tensor');
shading interp
if ~showTensor
  ht.Visible = 'off';
end

% Enclosing sphere
h = surface(xTensor,yTensor,zTensor,'Tag','enclose');
h.EdgeColor = 'none';
h.FaceColor = [1 1 1]*0.95;
h.FaceAlpha = 0.05;
light

% Planes
pz = patch('FaceColor','none');
pA = patch('FaceColor','r');
pB = patch('FaceColor','g');
set([pz pA pB],'FaceAlpha',0.1);

% Adjust view
axis equal
axis([-1 1 -1 1 -1 1]);
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
box off
axis off

set(gca,'Position',[0.2 0 0.8 1]);
set(gca,'CameraViewAngleMode','manual');
set(gca,'CameraViewAngle',8);
rotate3d(gca,'on');

view([167 16]);

calculateframes;
updateplot;

  function updateplot
    
    set(findobj('Tag','alphatext'),'String',sprintf('%2.3f',alpha));
    set(findobj('Tag','betatext'),'String',sprintf('%2.3f',beta));
    set(findobj('Tag','gammatext'),'String',sprintf('%2.3f',gamma));
    set(findobj('Tag','alphaslider'),'Value',alpha);
    set(findobj('Tag','betaslider'),'Value',beta);
    set(findobj('Tag','gammaslider'),'Value',gamma);

    % Names of the start and end frame
    setaxislabels();
    
    % Angles between all axis pairs
    angles = acos(R)*180/pi;
    for r = 1:3
      for c = 1:3
        %hMatrixElement(r,c).String = sprintf('%5.1f',angles(r,c));
        hMatrixElement(r,c).String = sprintf('%5.4f',R(r,c));
      end
    end
    
    % Set tensor ellipsoid orientation
    settensor();
    if showTensor
        st = 'on';
    else
        st = 'off';
    end
    set(findobj('Tag','tensor'),'Visible',st);
   
    v = {xA,yA,zA,xB,yB,zB,y1,-y1,x1,x2};
    for kk = 1:10
      p_ = v{kk};
      set(hline(kk),'XData',[0 p_(1)],'YData',[0 p_(2)],'ZData',[0 p_(3)]);
    end
    for kk = 1:6
      set(t(kk),'Position',v{kk}*1.1);
    end
    set(pz,'XData',zzplane(1,:),'YData',zzplane(2,:),'ZData',zzplane(3,:));
    set(pA,'XData',xyAplane(1,:),'YData',xyAplane(2,:),'ZData',xyAplane(3,:));
    set(pB,'XData',xyBplane(1,:),'YData',xyBplane(2,:),'ZData',xyBplane(3,:));
  end

  function rotview_updateangle(iangle,src)
    if strcmp(src.Style,'slider')
      newval = src.Value;
    elseif strcmp(src.Style,'edit')
      newval = str2double(src.String);
    end
    if newval<minangle(iangle)
      newval = minangle(iangle);
    end
    if newval>maxangle(iangle)
      newval = maxangle(iangle);
    end
    switch iangle
      case 1
        alpha = newval;
      case 2
        beta = newval;
      case 3
        gamma = newval;
    end
    calculateframes;
    updateplot;
  end

  function rotview_labelpopupmenu(~,~)
    h = findobj('Tag','labelpopupmenu');
    val = get(h,'Value');
    labelID = val;
    calculateframes;
    updateplot;
  end

  function rotview_toggletensor(~,~)
    h = findobj('Tag','tensorcheck');
    showTensor = h.Value;
    calculateframes;
    updateplot;
  end
  
  function calculateframes
    [xA,yA,zA] = erot([0 0 0]*pi/180,'rows');
    [x1,y1, ~] = erot([alpha 0 0]*pi/180,'rows');
    [x2, ~, ~] = erot([alpha beta 0]*pi/180,'rows');
    [xB,yB,zB] = erot([alpha beta gamma]*pi/180,'rows');
    R = [xB yB zB].';
    
    % (xA,yA), (xB,yB) and (zA,zB) planes
    nPhi = 101;
    phi = linspace(0,2*pi,nPhi);
    for kk = 1:numel(phi)
      zzplane(1:3,kk)  = cos(phi(kk))*zA + sin(phi(kk))*x1;
      xyAplane(1:3,kk) = cos(phi(kk))*xA + sin(phi(kk))*yA;
      xyBplane(1:3,kk) = cos(phi(kk))*xB + sin(phi(kk))*yB;
    end
  end
  
  function setaxislabels()
    switch labelID
      case 1, FrameBName = 'B'; FrameAName = 'A';
      case 2, FrameBName = 'g'; FrameAName = 'M';
      case 3, FrameBName = 'A'; FrameAName = 'M';
      case 4, FrameBName = 'D'; FrameAName = 'M';
      case 5, FrameBName = 'Q'; FrameAName = 'M';
      case 6, FrameBName = 'ee'; FrameAName = 'M';
      case 7, FrameBName = 'Diff'; FrameAName = 'M';
      case 8, FrameBName = 'M'; FrameAName = 'C';
      case 9, FrameBName = 'C'; FrameAName = 'L';
    end
    t(1).String = ['{\itx}_{' FrameAName '}'];
    t(2).String = ['{\ity}_{' FrameAName '}'];
    t(3).String = ['{\itz}_{' FrameAName '}'];
    t(4).String = ['{\itx}_{' FrameBName '}'];
    t(5).String = ['{\ity}_{' FrameBName '}'];
    t(6).String = ['{\itz}_{' FrameBName '}'];
    hRow(1).String = ['x_' FrameBName];
    hRow(2).String = ['y_' FrameBName];
    hRow(3).String = ['z_' FrameBName];
    hCol(1).String = ['x_' FrameAName];
    hCol(2).String = ['y_' FrameAName];
    hCol(3).String = ['z_' FrameAName];
  end

  function settensor()
    sz = size(xTensor);
    vecT = R*[cTensor(1)*xTensor(:) cTensor(2)*yTensor(:) cTensor(3)*zTensor(:)].';
    x_ = reshape(vecT(1,:),sz);
    y_ = reshape(vecT(2,:),sz);
    z_ = reshape(vecT(3,:),sz);
    h = findobj('Tag','tensor');
    set(h,'XData',x_,'YData',y_,'ZData',z_);
  end



end
