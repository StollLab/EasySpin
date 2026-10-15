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
% used as the starting Euler angles. The "all zero" and "random" buttons set all
% Euler angles to zero or to a random orientation, respectively.
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

% Quantities derived from the Euler angles, calculated in calculateframes()
xA = 0; yA = 0; zA = 0; xB = 0; yB = 0; zB = 0;
x1 = 0; x2 = 0; y1 = 0; zzPlane = 0; xyAPlane = 0; xyBPlane = 0;
R = [];
calculateframes();


% Figure and layout
%-------------------------------------------------------------------------------
% Close an existing rotview window, so there is only one
delete(findall(groot,'Type','figure','Tag','rotviewFigure'));

figSize = [1000 600];
screenSize = get(groot,'ScreenSize');
figSize = figSize*min([0.9*screenSize(3:4)./figSize 1]);
figPos = [ceil((screenSize(3:4)-figSize)/2) figSize]; % centered on the screen
hFig = uifigure('Tag','rotviewFigure','Name','Rotation viewer',...
  'Position',figPos);

% Fixed-width control panel on the left, plot panel filling the rest
controlWidth = 300;
mainGrid = uigridlayout(hFig,[1 2],'ColumnWidth',{controlWidth,'1x'},...
  'Padding',[5 5 5 5],'ColumnSpacing',5);
controlPanel = uipanel(mainGrid);
plotPanel = uipanel(mainGrid,'BackgroundColor','w');

% Controls, stacked from the top
%-------------------------------------------------------------------------------
rowHeight = 22;
nRows = 9;
controlGrid = uigridlayout(controlPanel,[nRows+1 3],...
  'RowHeight',[num2cell(repmat(rowHeight,1,nRows)) {'1x'}],...
  'ColumnWidth',{75,'1x',70});
controlGrid.RowHeight{8} = 5*rowHeight;  % rotation matrix display

% Euler angle sliders and edit fields
angleNames = {'alpha','beta','gamma'};
hSlider = gobjects(1,3);
hEdit = gobjects(1,3);
for iAngle = 1:3
  hl = uilabel(controlGrid,'Text',[angleNames{iAngle} ' (deg)']);
  hl.Layout.Row = iAngle;
  hl.Layout.Column = 1;
  hSlider(iAngle) = uislider(controlGrid,...
    'Limits',[minAngle maxAngle],'MajorTicks',[],'MinorTicks',[],...
    'ValueChangingFcn',@(~,evt) rotview_setangle(iAngle,evt.Value),...
    'ValueChangedFcn',@(src,~) rotview_setangle(iAngle,src.Value));
  hSlider(iAngle).Layout.Row = iAngle;
  hSlider(iAngle).Layout.Column = 2;
  hEdit(iAngle) = uieditfield(controlGrid,'numeric',...
    'Limits',[minAngle maxAngle],'ValueDisplayFormat','%.3f',...
    'ValueChangedFcn',@(src,~) rotview_setangle(iAngle,src.Value));
  hEdit(iAngle).Layout.Row = iAngle;
  hEdit(iAngle).Layout.Column = 3;
end

% Buttons for zero and random Euler angles
buttonGrid = uigridlayout(controlGrid,[1 2],'Padding',[0 0 0 0]);
buttonGrid.Layout.Row = 4;
buttonGrid.Layout.Column = [2 3];
uibutton(buttonGrid,'Text','all zero','ButtonPushedFcn',@(~,~) rotview_setangles([0 0 0]));
uibutton(buttonGrid,'Text','random','ButtonPushedFcn',@(~,~) rotview_randomangles());

% Frame labeling
labelStrings = {'general (A->B)','Sys.gFrame (M->g)','Sys.AFrame (M->A)',...
  'Sys.DFrame (M->D)','Sys.QFrame (M->Q)','Sys.eeFrame (M->ee)',...
  'Sys.DiffFrame (M->Diff)','Exp.MolFrame (C->M)','Exp.SampleFrame (L->C)'};
hl = uilabel(controlGrid,'Text','labeling');
hl.Layout.Row = 5;
hl.Layout.Column = 1;
hLabelMenu = uidropdown(controlGrid,'Items',labelStrings,...
  'ItemsData',1:numel(labelStrings),'Value',1,...
  'ValueChangedFcn',@(~,~) updateplot());
hLabelMenu.Layout.Row = 5;
hLabelMenu.Layout.Column = [2 3];

% Ellipsoid display
hc = uicheckbox(controlGrid,'Text','show ellipsoid','Value',showTensor,...
  'ValueChangedFcn',@rotview_toggletensor);
hc.Layout.Row = 6;
hc.Layout.Column = [1 3];

% Matrix display mode
hb = uibutton(controlGrid,'state','Text','show angles (deg)','Value',showAngles,...
  'ValueChangedFcn',@rotview_toggleangles);
hb.Layout.Row = 7;
hb.Layout.Column = [1 3];

% Rotation matrix display, with row and column labels
matrixGrid = uigridlayout(controlGrid,[4 4],'Padding',[0 0 0 0],...
  'RowHeight',{'1x','1x','1x','1x'},'ColumnWidth',{40,'1x','1x','1x'},...
  'RowSpacing',0);
matrixGrid.Layout.Row = 8;
matrixGrid.Layout.Column = [1 3];
hRowLabel = gobjects(1,3);
hColLabel = gobjects(1,3);
hMatrixElement = gobjects(3,3);
for k = 1:3
  hColLabel(k) = uilabel(matrixGrid,'Interpreter','tex');
  hColLabel(k).Layout.Row = 1;
  hColLabel(k).Layout.Column = k+1;
  hRowLabel(k) = uilabel(matrixGrid,'Interpreter','tex');
  hRowLabel(k).Layout.Row = k+1;
  hRowLabel(k).Layout.Column = 1;
end
set([hRowLabel hColLabel],'HorizontalAlignment','right');
for iRow = 1:3
  for iCol = 1:3
    hMatrixElement(iRow,iCol) = uilabel(matrixGrid,'HorizontalAlignment','right');
    hMatrixElement(iRow,iCol).Layout.Row = iRow+1;
    hMatrixElement(iRow,iCol).Layout.Column = iCol+1;
  end
end

% Import of rotation matrix from workspace
importGrid = uigridlayout(controlGrid,[1 4],'Padding',[0 0 0 0],...
  'ColumnWidth',{75,'1x',55,60});
importGrid.Layout.Row = 9;
importGrid.Layout.Column = [1 3];
uilabel(importGrid,'Text','import R');
hImportMenu = uidropdown(importGrid,'Items',{''});
hImportButton = uibutton(importGrid,'Text','Import',...
  'ButtonPushedFcn',@(~,~) rotview_importmatrix());
uibutton(importGrid,'Text','Refresh','ButtonPushedFcn',@(~,~) rotview_refreshlist());
rotview_refreshlist();

% Plotting
%-------------------------------------------------------------------------------
plotGrid = uigridlayout(plotPanel,[1 1],'Padding',[0 0 0 0],'BackgroundColor','w');
hAx = uiaxes(plotGrid);
hAx.Toolbar.Visible = 'off';
hAx.Interactions = rotateInteraction;

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
view(hAx,[167 16]);

% Light fixed relative to the viewer: reposition it whenever the view changes
hLight = camlight(hAx,'right');
addlistener(hAx,'View','PostSet',@(~,~) camlight(hLight,'right'));

% Low-level graphics functions (line, surface, patch) do not enable the
% mouse interactions, so enable them explicitly
enableDefaultInteractivity(hAx);

updateplot();

  function updateplot()

    for iAng = 1:3
      hEdit(iAng).Value = angles(iAng);
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
      hMatrixElement(iElement).Text = sprintf(fmt,values(iElement));
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

  function rotview_setangle(iAng,newValue)
    angles(iAng) = min(max(newValue,minAngle),maxAngle);
    calculateframes();
    updateplot();
  end

  function rotview_setangles(newAngles)
    angles = newAngles;
    calculateframes();
    updateplot();
  end

  function rotview_randomangles()
    % Uniformly distributed random orientation: cos(beta) is uniform in [-1,1]
    rotview_setangles([360*rand acosd(2*rand-1) 360*rand]);
  end

  function rotview_toggletensor(src,~)
    showTensor = src.Value;
    updateplot();
  end

  function rotview_toggleangles(src,~)
    showAngles = src.Value;
    if showAngles
      src.Text = 'show elements';
    else
      src.Text = 'show angles (deg)';
    end
    updateplot();
  end

  function rotview_refreshlist()
    % Keep the currently selected matrix selected if it is still available
    oldName = hImportMenu.Value;
    names = findrotmatrices();
    if isempty(names)
      hImportMenu.Items = {'(none found)'};
      hImportButton.Enable = 'off';
    else
      hImportMenu.Items = names;
      if any(strcmp(names,oldName))
        hImportMenu.Value = oldName;
      end
      hImportButton.Enable = 'on';
    end
  end

  function rotview_importmatrix()
    name = hImportMenu.Value;
    % Re-check the variable, since it might have changed since the last refresh
    if ~any(strcmp(findrotmatrices(),name))
      rotview_refreshlist();
      uialert(hFig,sprintf('''%s'' is no longer a rotation matrix in the workspace.',name),...
        'Rotation viewer','Icon','warning');
      return
    end
    rotview_setangles(eulang(double(evalin('base',name)))*180/pi);
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
    switch hLabelMenu.Value
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
      hColLabel(iAxis).Text = labelA;
      hRowLabel(iAxis).Text = labelB;
    end
  end

end
