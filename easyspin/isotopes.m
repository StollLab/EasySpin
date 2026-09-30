% isotopes   Graphical interface for nuclear isotope data
%
%   isotopes
%
%   Opens a window with a periodic table of the elements and a table of
%   nuclear isotope data: natural abundance, spin, nuclear g factor,
%   gyromagnetic ratio, electric quadrupole moment, and NMR frequency at
%   a given magnetic field.
%
%   Click an element to list its isotopes, or "all" to list all isotopes.
%   Checkboxes control whether unstable and nonmagnetic isotopes are shown.
%   The magnetic field can be entered directly or set to typical X-, Q-,
%   or W-band values using the buttons next to it.
%
%   Requires MATLAB R2021b or later.

function isotopes()

% Check MATLAB version (isMATLABReleaseOlderThan was introduced in R2020b)
if ~exist('isMATLABReleaseOlderThan','file') || isMATLABReleaseOlderThan("R2021b")
  error('isotopes requires MATLAB R2021b or later.');
end

% Settings
figTag = 'isotopesfig';
XbandField = 340;  % mT
QbandField = 1200;  % mT
WbandField = 3400;  % mT
Field = XbandField;
buttonFontSize = 14;

% Check for existing figure and raise
hFig = findall(0,'Type','figure','Tag',figTag);
if ~isempty(hFig)
  figure(hFig);
  return
end

% Read isotope data
data = readIsotopeDataFile;

figdata.fullData = data;

% GUI dimensions (pixels)
elementWidth = 36;
elementHeight = elementWidth;
border = 10;
spacing = 5;
xSpacing = elementWidth+spacing;
ySpacing = elementHeight+spacing;
classSpacing = 5;
tableSpacing = 8;  % extra space above table
tableHeight = 200;
bottomHeight = 30;
brightenFactor = 0.6;  % fraction of blending with white for selected element button

% Calculate window size
windowWidth = border + 17*xSpacing + elementWidth + 2*classSpacing + border;
windowHeight = border + 7*ySpacing + border + 2*ySpacing + border + ...
  tableSpacing + tableHeight + border + bottomHeight;

% Initialize figure window (hidden until fully built)
hFig = uifigure('Visible','off');
set(hFig,...
  'Position',[1 1 windowWidth windowHeight],...
  'Tag',figTag,...
  'Name','Nuclear isotopes',...
  'Toolbar','none',...
  'Menubar','none',...
  'NumberTitle','off',...
  'Resize','off');

% Add element buttons
ordNumber = 0;
yOff = windowHeight-elementHeight-border;
xOff = border;
for k = 1:numel(data.gn)
  if data.Z(k)<=ordNumber, continue; end
  ordNumber = data.Z(k);
  p = [0 0 elementWidth elementHeight];
  [period,group,cl] = elementclass(ordNumber);
  if cl==2
    p(2) = yOff - (period+1)*ySpacing - classSpacing;
    p(1) = xOff + (group-1)*xSpacing + classSpacing;
  else
    p(2) = yOff - (period-1)*ySpacing;
    p(1) = xOff + (group-1)*xSpacing;
    if group>2, p(1) = p(1) + classSpacing; end
    if group>12, p(1) = p(1) + classSpacing; end
  end
  hButton = uibutton(hFig);
  set(hButton,...
    'Position',p,...
    'ButtonPushedFcn',@elementButtonPushedCallback,...
    'Text',data.element{k},...
    'FontSize',buttonFontSize,...
    'Tooltip',[' ' data.name{k} ' ']);
  switch cl
    case 0
      if group<3
        bgcol = [99 154 255]/255;
      else
        bgcol = [255 207 0]/255;
      end
    case 1
      bgcol = [255 154 156]/255;
    case 2
      bgcol = [0 207 49]/255;
  end
  if data.N(k)<=0
    % No isotope data: keep default color and disable
    hButton.Enable = 'off';
  else
    hButton.BackgroundColor = bgcol;
  end
  hButton.UserData = hButton.BackgroundColor;  % unselected color
end

% Add selection button for all elements
hAll = uibutton(hFig);
p = [xOff+16*xSpacing+2*classSpacing yOff-8*ySpacing-classSpacing ...
     elementWidth+xSpacing elementHeight+ySpacing];
set(hAll,...
  'Position',p,...
  'Text','all',...
  'BackgroundColor',[1 1 1]*0.9,...
  'ButtonPushedFcn',@elementButtonPushedCallback,...
  'FontSize',buttonFontSize,...
  'Tooltip','all elements');
hAll.UserData = hAll.BackgroundColor;  % unselected color

% Mark "all" as selected initially
hAll.BackgroundColor = hAll.UserData + (1-hAll.UserData)*brightenFactor;
hAll.FontWeight = 'bold';
figdata.hSelectedButton = hAll;
figdata.brightenFactor = brightenFactor;

% Add checkbox for unstable isotopes
xpos = xOff;
hUnstableCheckbox = uicheckbox(hFig);
set(hUnstableCheckbox,...
  'Position',[xpos border 160 22],...
  'Text','Show unstable isotopes',...
  'Value',0,...
  'ValueChangedFcn',@(~,~)updateTable(hFig));
figdata.hUnstableCheckbox = hUnstableCheckbox;

% Add checkbox for nonmagnetic isotopes
xpos = xpos+170;
hNonmagneticCheckbox = uicheckbox(hFig);
set(hNonmagneticCheckbox,...
  'Position',[xpos border 180 22],...
  'Text','Show nonmagnetic isotopes',...
  'Value',1,...
  'ValueChangedFcn',@(~,~)updateTable(hFig));
figdata.hNonmagneticCheckbox = hNonmagneticCheckbox;

% Magnetic field edit box with label (right-aligned with band buttons)
fieldControlsWidth = 110 + 10 + 100 + 5 + 3*30 + 2*5;  % label, edit box, band buttons
xpos = windowWidth - border - fieldControlsWidth;
uilabel(hFig,...
  'Text','Magnetic field (mT)',...
  'VerticalAlignment','center',...
  'Position',[xpos border 110 22]);
hFieldEdit = uieditfield(hFig,'numeric');
xpos = xpos+120;
set(hFieldEdit,...
  'BackgroundColor','white',...
  'Position',[xpos border 100 22],...
  'Value',Field,...
  'Limits',[0 Inf],...
  'HorizontalAlignment','left',...
  'ValueChangedFcn',@(~,~)updateTable(hFig));
figdata.hFieldEdit = hFieldEdit;

% Add convenience buttons for X, Q and W band fields
xpos = xpos + 105;
hXbandButton = uibutton(hFig);
set(hXbandButton,...
  'Position',[xpos border 30 22],...
  'Text','X',...
  'ButtonPushedFcn',@(~,~)setField(hFieldEdit,XbandField));
xpos = xpos + 35;
hQbandButton = uibutton(hFig);
set(hQbandButton,...
  'Position',[xpos border 30 22],...
  'Text','Q',...
  'ButtonPushedFcn',@(~,~)setField(hFieldEdit,QbandField));
xpos = xpos + 35;
hWbandButton = uibutton(hFig);
set(hWbandButton,...
  'Position',[xpos border 30 22],...
  'Text','W',...
  'ButtonPushedFcn',@(~,~)setField(hFieldEdit,WbandField));

% Table of isotope data
hTable = uitable(hFig);
set(hTable,...
  'Position',[xOff border+bottomHeight windowWidth-2*border tableHeight]);
figdata.hTable = hTable;

tabledata = data(:,{'isotope','abundance','spin','gn','gamma','qm'});
tabledata.NMRfreq = zeros(height(tabledata),1);

hTable.Data = tabledata;
hTable.ColumnWidth = 'auto';
hTable.SelectionType = 'row';
hTable.ColumnName = {'Isotope','Abundance (%)','Spin',...
  'gn value','γ/2π (MHz/T)','Q (barn)','Frequency (MHz)'};
hTable.ColumnSortable = true;

figdata.Element = '';
figdata.tableData = tabledata;

guidata(hFig,figdata);
updateTable(hFig);

movegui(hFig,'center');
hFig.Visible = 'on';

end


%-------------------------------------------------------------------------------
function elementButtonPushedCallback(src,~)
Element = src.Text;
hFig = src.Parent;
data = guidata(hFig);

if Element=="all", Element = ''; end
data.Element = Element;

% Highlight selected button, restore previously selected one
data.hSelectedButton.BackgroundColor = data.hSelectedButton.UserData;
data.hSelectedButton.FontWeight = 'normal';
src.BackgroundColor = src.UserData + (1-src.UserData)*data.brightenFactor;
src.FontWeight = 'bold';
data.hSelectedButton = src;

guidata(hFig,data);

updateTable(hFig);

end


%-------------------------------------------------------------------------------
function updateTable(hFig)
data = guidata(hFig);
hTable = data.hTable;

element = data.Element;

% Update NMR frequencies
B0 = data.hFieldEdit.Value;
data.tableData.NMRfreq = B0*1e-3*nmagn*data.tableData.gn/planck/1e6;

% Filter table by element
if isempty(element)
  idx = true(height(data.tableData),1);
else
  idx = data.fullData.element==string(element);
end

% Hide unstable isotopes if desired
if ~data.hUnstableCheckbox.Value
  idx = idx & data.fullData.radioactive=="-";
end

% Hide nonmagnetic isotopes if desired
if ~data.hNonmagneticCheckbox.Value
  idx = idx & data.fullData.spin~=0;
end

% Hide elements without any isotopes
idx = idx & data.fullData.spin>=0;

if ~any(idx)
  hTable.Data = [];
else
  hTable.Data = data.tableData(idx,:);
end

end


%-------------------------------------------------------------------------------
function data = readIsotopeDataFile

% Determine full data file name
esPath = fileparts(which(mfilename));
DataFile = [esPath filesep 'private' filesep 'isotopedata.txt'];

% Load data
fh = fopen(DataFile);
if fh<0
  error('Could not open nuclear isotopes data file %s',DataFile);
end
C = textscan(fh,'%f %f %s %s %s %f %f %f %f','commentstyle','%');
fclose(fh);

% Calculate gyromagnetic ratios (MHz/T)
magmom = C{7};
% Nuclear g factors: Infer from the 1H entry whether gn or gn*I is listed in
% the data file. If gn*I is listed, then divide out I.
gn = magmom;
spin = C{6};
if magmom(1)<5
  nzidx = spin~=0;
  gn(nzidx) = magmom(nzidx)./spin(nzidx);
end
C{7} = gn;
C{10} = gn*nmagn/planck/1e6;

% Assemble isotope symbols
N = C{2};
element = C{4};
radioactive = C{3};
for k = numel(N):-1:1
  isostr = sprintf('%d%s',N(k),element{k});
  if radioactive{k}=='*'
    isoSymbols{k} = [isostr '*'];
  else
    isoSymbols{k} = isostr;
  end
end
C{11} = isoSymbols(:);

% Construct table
vnames = {'Z','N','radioactive','element','name','spin',...
  'gn','abundance','qm','gamma','isotope'};
data = table(C{:},'VariableNames',vnames);

end


%-------------------------------------------------------------------------------
function [Period,Group,Class] = elementclass(N)

Class = 0;

periodLimits = [0 2 10 18 36 54 86 1000];

% Determine period of element
for Period = 1:8
  if N<=periodLimits(Period), break; end
end
Period = Period - 1;

%Determine group and class of element
%Class 0 - main groups, 1 - transition metals, 2 - rare earths
Group = N - periodLimits(Period);
switch Period
case 1
  Class = 0;
  if Group~=1, Group=18; end
case {2,3}
  Class = 0;
  if Group>2, Group = Group + 10; end
case {4,5}
  Class = 1;
  if Group<3 || Group>12, Class = 0; end
case {6,7}
  if Group<3 || Group>26
    Class = 0;
  else
    if Group>16
      Class = 1;
    else
      Class = 2;
    end
  end
  if Class<2 && Group>16
    Group = Group - 14;
  end
end

end


%-------------------------------------------------------------------------------
function setField(hFieldEdit,Field)

hFieldEdit.Value = Field;  % mT

updateTable(ancestor(hFieldEdit,'figure'));

end
