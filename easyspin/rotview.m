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
%------------------------------------------------------------------
% Three Euler angles, in degrees
% for R = erot([alpha beta gamma]*pi/180);
switch nargin
  case 1
    if all(size(varargin{1})==[3 3])
      R = varargin{1};
      [alpha,beta,gamma] = eulang(R);
      alpha = alpha*180/pi;
      beta = beta*180/pi;
      gamma = gamma*180/pi;
    elseif numel(varargin{1})==3
      angles_in = varargin{1}*180/pi; % rad -> deg
      alpha = angles_in(1);
      beta = angles_in(2);
      gamma = angles_in(3);
    else
      error('For one input argument, it must be a 3x3 rotation matrix or a 3-element vector of Euler angles.');
    end
  case 3
    alpha = varargin{1}*180/pi;
    beta = varargin{2}*180/pi;
    gamma = varargin{3}*180/pi;
  otherwise
    alpha = 30;
    beta = 20;
    gamma = 40;
end

showTensor = true;

% Initialize quantities derived from the Euler angles
xA = 0; yA = 0; zA = 0; xB = 0; yB = 0; zB = 0;
x1 = 0; x2 = 0; y1 = 0; zzplane = 0; xyAplane = 0; xyBplane = 0;
R = [];

calculateframes();


% UI setup
%------------------------------------------------------------------
hFig = figure(443123);
clf(hFig)
set(hFig,'Name','Rotation viewer [EasySpin]',...
  'NumberTitle','off',...
  'WindowStyle','normal',...
  'MenuBar','none',...
  'Color','w');
p = get(gcf,'Position');
set(gcf,'Position',[p(1:2)-[120 100] 800 500]);

x0 = 10;
y0 = 10;
uicontrol('Style','slider','Tag','alphaslider',...
  'Position',[x0+50 y0 120 15],...
  'Min',0,'Max',360,'Value',alpha,'SliderStep',[1 10]/360,...
  'callback',@rotview_alphaslider);
uicontrol('Style','text','Position',[x0,y0,50,15],...
  'String','alpha (deg)','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','edit','Tag','alphatext',...
  'Position',[x0+180,y0,40,15],...
  'String','','HorizontalA','left');

uicontrol('Style','slider','Tag','betaslider',...
  'Position',[x0+50 y0+20 120 15],...
  'Min',0,'Max',180,'Value',beta,'SliderStep',[1 10]/180,...
  'callback',@rotview_betaslider);
uicontrol('Style','text','Position',[x0,y0+20,50,15],...
  'String','beta (deg)','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','edit','Tag','betatext',...
  'Position',[x0+180,y0+20,40,15],...
  'String','','HorizontalA','left');

uicontrol('Style','slider','Tag','gammaslider',...
  'Position',[x0+50 y0+40 120 15],...
  'Min',0,'Max',360,'Value',gamma,'SliderStep',[1 10]/360,...
  'callback',@rotview_gammaslider);
uicontrol('Style','text','Position',[x0,y0+40,50,15],...
  'String','gamma (deg)','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','edit','Tag','gammatext',...
  'Position',[x0+180,y0+40,40,15],...
  'String','','HorizontalA','left');

labelStrings = {'gFrame','AFrame','DFrame','QFrame','eeFrame','DiffFrame','MolFrame'};
labelID = 1;
uicontrol('Style','text','Position',[x0,y0+65,50,15],...
  'String','scenario','HorizontalA','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','popupmenu','Tag','labelpopupmenu',...
  'Position',[x0+50 y0+70 120 15],...
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


xR0 = x0+20;
yR0 = y0+110;
dx = 50;
dy = 20;
for a = 1:3
    for b = 1:3
        hAngle(a,b) = uicontrol(...
            'Style','text',...
            'HorizontalAl','right',...
            'Position',[xR0+(a-1)*dx,yR0+(3-b)*dy,dx,dy]);
    end
end
set(hAngle,'Background','w','ForegroundColor','k');

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
h.FaceColor = [1 1 1]*0.4;
h.FaceAlpha = 0.1;
light

% Planes
pz = patch('FaceColor','none');
pA = patch('FaceColor','r');
pB = patch('FaceColor','g');
set([pz pA pB],'FaceAlpha',0.05);

% Adjust view
axis equal
axis([-1 1 -1 1 -1 1]);
box on
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
box off
axis off

set(gca,'Position',[0 0 1 1]);
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

    % Names of the start and end frame
    setaxislabels();
    
    % Angles between all axis pairs
    angles = acos(R)*180/pi;
    for a = 1:3
        for b = 1:3
            hAngle(a,b).String = sprintf('%5.1f',angles(a,b));
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
    for k=1:10
      p_ = v{k};
      set(hline(k),'XData',[0 p_(1)],'YData',[0 p_(2)],'ZData',[0 p_(3)]);
    end
    for k=1:6
      set(t(k),'Position',v{k}*1.1);
    end
    set(pz,'XData',zzplane(1,:),'YData',zzplane(2,:),'ZData',zzplane(3,:));
    set(pA,'XData',xyAplane(1,:),'YData',xyAplane(2,:),'ZData',xyAplane(3,:));
    set(pB,'XData',xyBplane(1,:),'YData',xyBplane(2,:),'ZData',xyBplane(3,:));
  end

  function rotview_alphaslider(~,~)
    h = findobj('Tag','alphaslider');
    val = get(h,'Value');
    alpha = val;
    calculateframes;
    updateplot;
  end

  function rotview_betaslider(~,~)
    h = findobj('Tag','betaslider');
    val = get(h,'Value');
    beta = val;
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

  function rotview_gammaslider(~,~)
    h = findobj('Tag','gammaslider');
    val = get(h,'Value');
    gamma = val;
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
    R = [xB yB zB];
    
    % (xA,yA), (xB,yB) and (zA,zB) planes
    nPhi = 101;
    phi = linspace(0,2*pi,nPhi);
    for k=1:numel(phi)
      zzplane(1:3,k)  = cos(phi(k))*zA + sin(phi(k))*x1;
      xyAplane(1:3,k) = cos(phi(k))*xA + sin(phi(k))*yA;
      xyBplane(1:3,k) = cos(phi(k))*xB + sin(phi(k))*yB;
    end
  end
  
  function setaxislabels()
    switch labelID
      case 1, FrameBName = 'g'; FrameAName = 'M';
      case 2, FrameBName = 'A'; FrameAName = 'M';
      case 3, FrameBName = 'D'; FrameAName = 'M';
      case 4, FrameBName = 'Q'; FrameAName = 'M';
      case 5, FrameBName = 'ee'; FrameAName = 'M';
      case 6, FrameBName = 'Diff'; FrameAName = 'M';
      case 7, FrameBName = 'M'; FrameAName = 'C';
    end
    t(1).String = ['{\itx}_{' FrameAName '}'];
    t(2).String = ['{\ity}_{' FrameAName '}'];
    t(3).String = ['{\itz}_{' FrameAName '}'];
    t(4).String = ['{\itx}_{' FrameBName '}'];
    t(5).String = ['{\ity}_{' FrameBName '}'];
    t(6).String = ['{\itz}_{' FrameBName '}'];
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
