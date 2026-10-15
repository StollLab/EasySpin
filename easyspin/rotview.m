% rotview    Display rotation between two frames
%
%    rotview()
%    rotview([a b c])
%    rotview(a,b,c)
%    rotview(R)
%
% Plots two frames and allows adjustment of their relative orientation.
% Euler angles and the rotation matrix are shown. The matrix elements can
% also be shown as angles between axis pairs. If inputs are given, they are
% used as the starting Euler angles.
%
% A rotation matrix from the workspace can be imported by selecting it in the
% "import R" menu and clicking "Import". The menu lists all real 3x3 rotation
% matrices in the workspace. Click "Refresh" after defining new ones.
%
% Input:
%    a, b, c:    Euler angles alpha, beta, and gamma (in radians)
%    R:          rotation matrix, 3x3
%

function rotview(varargin)

% Parameters
%-------------------------------------------------------------------------------
switch nargin
  case 0
    angles = [30 20 40]*pi/180;
  case 1
    if isequal(size(varargin{1}),[3 3])
      angles = eulang(varargin{1});
    elseif numel(varargin{1})==3
      angles = varargin{1}(:).';
    else
      error('For one input argument, it must be a 3x3 rotation matrix or a 3-element vector of Euler angles.');
    end
  case 3
    angles = [varargin{1} varargin{2} varargin{3}];
  otherwise
    error('rotview takes 0, 1, or 3 input arguments.');
end

% Euler angles in degrees, limited to the slider range
minAngle = -360;
maxAngle = 360;
angles = min(max(angles*180/pi,minAngle),maxAngle);

showTensor = true;
showAngles = false;
labelID = 1;

% Quantities derived from the Euler angles, calculated in calculateframes()
xA = 0; yA = 0; zA = 0; xB = 0; yB = 0; zB = 0;
x1 = 0; x2 = 0; y1 = 0; zzPlane = 0; xyAPlane = 0; xyBPlane = 0;
R = [];
calculateframes();


% UI setup
%-------------------------------------------------------------------------------
figSize = [1000 600];
screenSize = get(0,'ScreenSize');
scaleFactor = min([0.9*screenSize(3:4)./figSize 1]);
figSize = figSize*scaleFactor;
xPos = ceil((screenSize(3)-figSize(1))/2); % center the figure on the screen
yPos = ceil((screenSize(4)-figSize(2))/2);
hFig = figure(443123);
clf(hFig)
set(hFig,'Units','pixels','Position',[xPos yPos figSize],'Resize','on',...
  'Name','Rotation viewer','NumberTitle','off','WindowStyle','normal',...
  'MenuBar','none','Color','w');

x0 = 10;
y0 = 10;

% Euler angle sliders and edit fields
angleNames = {'alpha','beta','gamma'};
hSlider = gobjects(1,3);
hEdit = gobjects(1,3);
for iAngle = 1:3
  yi = y0+(3-iAngle)*20;
  uicontrol('Style','text','Position',[x0 yi 50 15],...
    'String',[angleNames{iAngle} ' (deg)'],'HorizontalAlignment','left',...
    'Background','w','ForegroundColor','k');
  hSlider(iAngle) = uicontrol('Style','slider',...
    'Position',[x0+50 yi 120 15],...
    'Min',minAngle,'Max',maxAngle,'Value',angles(iAngle),...
    'SliderStep',[1 10]/(maxAngle-minAngle),...
    'Callback',@(src,~) rotview_updateangle(iAngle,src));
  hEdit(iAngle) = uicontrol('Style','edit',...
    'Position',[x0+180 yi 50 15],...
    'String','','HorizontalAlignment','left',...
    'Callback',@(src,~) rotview_updateangle(iAngle,src));
end

% Frame labeling
labelStrings = {'general (A->B)','Sys.gFrame (M->g)','Sys.AFrame (M->A)',...
  'Sys.DFrame (M->D)','Sys.QFrame (M->Q)','Sys.eeFrame (M->ee)',...
  'Sys.DiffFrame (M->Diff)','Exp.MolFrame (C->M)','Exp.SampleFrame (L->C)'};
uicontrol('Style','text','Position',[x0 y0+65 50 15],...
  'String','labeling','HorizontalAlignment','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','popupmenu',...
  'Position',[x0+50 y0+70 180 15],...
  'String',labelStrings,'Value',labelID,...
  'Callback',@rotview_labelpopupmenu);

% Ellipsoid display and matrix display mode
uicontrol('Style','text','Position',[x0 y0+90 50 15],...
  'String','ellipsoid','HorizontalAlignment','left',...
  'Background','w','ForegroundColor','k');
uicontrol('Style','checkbox',...
  'Position',[x0+50 y0+90 20 15],...
  'Value',showTensor,'Background','w',...
  'Callback',@rotview_toggletensor);
uicontrol('Style','togglebutton',...
  'Position',[x0+130 y0+88 100 19],...
  'String','show angles (deg)','Value',showAngles,...
  'Callback',@rotview_toggleangles);

% Rotation matrix display
%-------------------------------------------------------------------------------
xR0 = x0+50;
yR0 = y0+110;
dx = 50;
dy = 20;
% Row and column labels are text objects in an invisible axes, since text
% uicontrols do not support TeX subscripts
hLabelAxes = axes(hFig,'Units','pixels','Position',[xR0-dx yR0 4*dx 4*dy],...
  'XLim',[0 4*dx],'YLim',[0 4*dy],'Visible','off','HitTest','off');
hRowLabel = gobjects(1,3);
hColLabel = gobjects(1,3);
for k = 1:3
  hRowLabel(k) = text(hLabelAxes,dx-4,(3-k+0.5)*dy,'');
  hColLabel(k) = text(hLabelAxes,(k+1)*dx-4,3.5*dy,'');
end
set([hRowLabel hColLabel],'HorizontalAlignment','right','Color','k',...
  'FontSize',get(0,'DefaultUicontrolFontSize')+1,'HitTest','off');

hMatrixElement = gobjects(3,3);
for iRow = 1:3
  for iCol = 1:3
    hMatrixElement(iRow,iCol) = uicontrol('Style','text',...
      'Position',[xR0+(iCol-1)*dx yR0+(3-iRow)*dy dx dy],...
      'HorizontalAlignment','right','Background','w','ForegroundColor','k');
  end
end

% Import of rotation matrix from workspace
%-------------------------------------------------------------------------------
uicontrol('Style','text','Position',[x0 y0+215 50 15],...
  'String','import R','HorizontalAlignment','left',...
  'Background','w','ForegroundColor','k');
hImportMenu = uicontrol('Style','popupmenu',...
  'Position',[x0+50 y0+220 100 15],...
  'String',{''});
hImportButton = uicontrol('Style','pushbutton',...
  'Position',[x0+155 y0+217 40 20],...
  'String','Import',...
  'Callback',@rotview_importmatrix);
uicontrol('Style','pushbutton',...
  'Position',[x0+200 y0+217 45 20],...
  'String','Refresh',...
  'Callback',@(~,~) rotview_refreshlist());
rotview_refreshlist();

% Plotting
%-------------------------------------------------------------------------------
hAx = axes(hFig);

% Frame axes: A dashed, B solid, plus nodal line etc. in black
axisColors = {'r',[0 0.6 0],'b'};
hLine = gobjects(1,10);
for k = 1:3
  hLine(k) = line(hAx,'Color',axisColors{k},'LineStyle','--');
  hLine(k+3) = line(hAx,'Color',axisColors{k},'LineStyle','-');
end
set(hLine(1:6),'LineWidth',2);
for k = 7:10
  hLine(k) = line(hAx,'Color','k');
end

% Axis labels
hAxisLabel = gobjects(1,6);
for k = 1:6
  hAxisLabel(k) = text(hAx,0,0,'');
end

% Tensor ellipsoid
[xTensor,yTensor,zTensor] = sphere(100);
cTensor = [0.2 0.3 0.4];
hTensor = surface(hAx,cTensor(1)*xTensor,cTensor(2)*yTensor,cTensor(3)*zTensor,...
  'EdgeColor','none','FaceColor','interp');

% Enclosing sphere
surface(hAx,xTensor,yTensor,zTensor,...
  'EdgeColor','none','FaceColor',[1 1 1]*0.95,'FaceAlpha',0.05);
light(hAx);

% Planes
hzzPlane = patch(hAx,'FaceColor','none');
hxyAPlane = patch(hAx,'FaceColor','r');
hxyBPlane = patch(hAx,'FaceColor','g');
set([hzzPlane hxyAPlane hxyBPlane],'FaceAlpha',0.1);

% Adjust view
axis(hAx,'equal');
axis(hAx,[-1 1 -1 1 -1 1]);
axis(hAx,'off');
set(hAx,'CameraViewAngleMode','manual','CameraViewAngle',8);
% Controls stay fixed at the bottom left, the plot fills the rest of the window
hFig.SizeChangedFcn = @(~,~) rotview_resize();
rotview_resize();
view(hAx,[167 16]);
hRotate = rotate3d(hFig);
setAllowAxesRotate(hRotate,hLabelAxes,false);
hRotate.Enable = 'on';

updateplot();

  function updateplot()

    for iAng = 1:3
      hEdit(iAng).String = sprintf('%2.3f',angles(iAng));
      hSlider(iAng).Value = angles(iAng);
    end

    setaxislabels();

    % Matrix elements, or angles between all axis pairs (in degrees)
    if showAngles
      values = acosd(min(max(R,-1),1)); % clamp to avoid complex values from rounding
      fmt = '%5.1f';
    else
      values = R;
      fmt = '%5.4f';
    end
    for iElement = 1:9
      hMatrixElement(iElement).String = sprintf(fmt,values(iElement));
    end

    % Tensor ellipsoid, with principal axes along the axes of frame B
    if showTensor
      sz = size(xTensor);
      vecT = R.'*[cTensor(1)*xTensor(:) cTensor(2)*yTensor(:) cTensor(3)*zTensor(:)].';
      set(hTensor,'XData',reshape(vecT(1,:),sz),'YData',reshape(vecT(2,:),sz),...
        'ZData',reshape(vecT(3,:),sz),'Visible','on');
    else
      hTensor.Visible = 'off';
    end

    v = {xA,yA,zA,xB,yB,zB,y1,-y1,x1,x2};
    for iLine = 1:10
      set(hLine(iLine),'XData',[0 v{iLine}(1)],'YData',[0 v{iLine}(2)],...
        'ZData',[0 v{iLine}(3)]);
    end
    for iLabel = 1:6
      hAxisLabel(iLabel).Position = v{iLabel}*1.1;
    end
    set(hzzPlane,'XData',zzPlane(1,:),'YData',zzPlane(2,:),'ZData',zzPlane(3,:));
    set(hxyAPlane,'XData',xyAPlane(1,:),'YData',xyAPlane(2,:),'ZData',xyAPlane(3,:));
    set(hxyBPlane,'XData',xyBPlane(1,:),'YData',xyBPlane(2,:),'ZData',xyBPlane(3,:));
  end

  function rotview_updateangle(iAng,src)
    if strcmp(src.Style,'slider')
      newValue = src.Value;
    else
      newValue = str2double(src.String);
    end
    if ~isnan(newValue)
      angles(iAng) = min(max(newValue,minAngle),maxAngle);
      calculateframes();
    end
    updateplot(); % also restores the edit field after invalid input
  end

  function rotview_resize()
    plotLeft = 200; % pixels
    figPos = getpixelposition(hFig);
    hAx.Units = 'pixels';
    hAx.Position = [plotLeft 0 max(figPos(3)-plotLeft,1) max(figPos(4),1)];
  end

  function rotview_labelpopupmenu(src,~)
    labelID = src.Value;
    updateplot();
  end

  function rotview_toggletensor(src,~)
    showTensor = logical(src.Value);
    updateplot();
  end

  function rotview_toggleangles(src,~)
    showAngles = logical(src.Value);
    if showAngles
      src.String = 'show elements';
    else
      src.String = 'show angles (deg)';
    end
    updateplot();
  end

  function rotview_refreshlist()
    % Keep the currently selected matrix selected if it is still available
    oldNames = cellstr(hImportMenu.String);
    oldName = oldNames{hImportMenu.Value};
    names = findrotmatrices();
    if isempty(names)
      hImportMenu.String = {'(none found)'};
      hImportMenu.Value = 1;
      hImportButton.Enable = 'off';
    else
      idx = find(strcmp(names,oldName),1);
      if isempty(idx), idx = 1; end
      hImportMenu.String = names;
      hImportMenu.Value = idx;
      hImportButton.Enable = 'on';
    end
  end

  function rotview_importmatrix(~,~)
    names = cellstr(hImportMenu.String);
    name = names{hImportMenu.Value};
    % Re-check the variable, since it might have changed since the last refresh
    if ~any(strcmp(findrotmatrices(),name))
      rotview_refreshlist();
      warndlg(sprintf('''%s'' is no longer a rotation matrix in the workspace.',name),...
        'Rotation viewer','modal');
      return
    end
    angles = eulang(double(evalin('base',name)))*180/pi;
    calculateframes();
    updateplot();
  end

  function names = findrotmatrices()
    % Find all real 3x3 rotation matrices in the base workspace
    vars = evalin('base','whos');
    names = {};
    for iVar = 1:numel(vars)
      vinfo = vars(iVar);
      if ~isequal(vinfo.size,[3 3]) || ~any(strcmp(vinfo.class,{'double','single'})) || ...
          vinfo.complex || vinfo.sparse
        continue
      end
      M = double(evalin('base',vinfo.name));
      % same thresholds as in eulang
      if all(isfinite(M(:))) && norm(M.'*M-eye(3))<1e-2 && det(M)>0
        names{end+1} = vinfo.name; %#ok<AGROW>
      end
    end
  end

  function calculateframes()
    anglesRad = angles*pi/180;
    [xA,yA,zA] = erot([0 0 0],'rows');
    [x1,y1,~] = erot([anglesRad(1) 0 0],'rows');
    [x2,~,~] = erot([anglesRad(1:2) 0],'rows');
    [xB,yB,zB] = erot(anglesRad,'rows');
    R = [xB yB zB].';

    % (zA,zB), (xA,yA) and (xB,yB) planes
    phi = linspace(0,2*pi,101);
    zzPlane  = zA*cos(phi) + x1*sin(phi);
    xyAPlane = xA*cos(phi) + yA*sin(phi);
    xyBPlane = xB*cos(phi) + yB*sin(phi);
  end

  function setaxislabels()
    switch labelID
      case 1, frameA = 'A'; frameB = 'B';
      case 2, frameA = 'M'; frameB = 'g';
      case 3, frameA = 'M'; frameB = 'A';
      case 4, frameA = 'M'; frameB = 'D';
      case 5, frameA = 'M'; frameB = 'Q';
      case 6, frameA = 'M'; frameB = 'ee';
      case 7, frameA = 'M'; frameB = 'Diff';
      case 8, frameA = 'C'; frameB = 'M';
      case 9, frameA = 'L'; frameB = 'C';
    end
    xyz = 'xyz';
    for iAxis = 1:3
      labelA = ['{\it' xyz(iAxis) '}_{' frameA '}'];
      labelB = ['{\it' xyz(iAxis) '}_{' frameB '}'];
      hAxisLabel(iAxis).String = labelA;
      hAxisLabel(iAxis+3).String = labelB;
      hColLabel(iAxis).String = labelA;
      hRowLabel(iAxis).String = labelB;
    end
  end

end
