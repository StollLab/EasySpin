% isotopes   Graphical interface for nuclear isotope data 
%
%   Displays a periodic table of the elements with
%   a table of nuclear isotopes.
%

function isotopes()

% Settings
figTag = 'isotopesfig';
Field = 340;  % mT
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
figdata.DefaultField = Field;

% GUI dimensions (pixels)
elementWidth = 36;
elementHeight = elementWidth;
border = 10;
spacing = 5;
xSpacing = elementWidth+spacing;
ySpacing = elementHeight+spacing;
classSpacing = 5;
labelHeight = 15;
tableHeight = 200;
bottomHeight = 30;

% Calculate and set window size
screensize = get(0, 'ScreenSize');
screenWidth = screensize(3);
screenHeight = screensize(4);
windowWidth = border + 18*xSpacing + 2*classSpacing + border;
windowHeight = border + 7*ySpacing + border + 2*ySpacing + border + ...
  + labelHeight/2 + tableHeight + border + bottomHeight;
figPos(1) = (screenWidth-windowWidth)/2;
figPos(2) = (screenHeight-windowHeight)/2;
figPos(3) = windowWidth;
figPos(4) = windowHeight;

% Initialize figure window
hFig = uifigure();
set(hFig,...
  'Position',figPos,...
  'Tag',figTag,...
  'Name','Nuclear isotopes',...
  'Toolbar','none',...
  'Menubar','none',...
  'NumberTitle','off',...
  'Resize','off');

% Add element buttons
ordNumber = 0;
yOff = figPos(4)-elementHeight-border;
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
    'FontSize',buttonFontSize);
  if ~verLessThan('matlab','9.5')  % R2018b = 9.5
    set(hButton,'Tooltip',[' ' data.name{k} ' ']);
  end
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
    bgcol = get(hButton,'BackgroundColor');
  end
  hButton.BackgroundColor = bgcol;
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
  'FontSize',buttonFontSize);

if ~verLessThan('Matlab','9.5')  % R2018b = 9.5
  set(hAll,'Tooltip','all elements');
end

% Add checkbox for unstable isotopes
xpos = xOff;
hUnstableCheckbox = uicheckbox(hFig);
set(hUnstableCheckbox,...
  'Position',[xpos border 160 22],...
  'Text','Show unstable isotopes',...
  'Value',0,...
  'ValueChangedFcn',@(~,~)updateTable);
figdata.hUnstableCheckbox = hUnstableCheckbox;

% Add checkbox for nonmagnetic isotopes
xpos = xpos+170;
hNonmagneticCheckbox = uicheckbox(hFig);
set(hNonmagneticCheckbox,...
  'Position',[xpos border 180 22],...
  'Text','Show nonmagnetic isotopes',...
  'Value',1,...
  'ValueChangedFcn',@(~,~)updateTable);
figdata.hNonmagneticCheckbox = hNonmagneticCheckbox;

% Magnetic field edit box with label
xpos = xpos+255;
uilabel(hFig,...
  'Text','Magnetic field (mT)',...
  'Position',[xpos border 110 19]);
hFieldEdit = uieditfield(hFig,'numeric');
xpos = xpos+120;
set(hFieldEdit,...
  'BackgroundColor','white',...
  'Position',[xpos border 100 22],...
  'Value',Field,...
  'HorizontalAlignment','left',...
  'ValueChangedFcn',@(~,~)updateTable);
figdata.hFieldEdit = hFieldEdit;

% Add convenience buttons for X, Q and W band fields
xpos = xpos + 105;
hXbandButton = uibutton(hFig);
set(hXbandButton,...
  'Position',[xpos border 30 22],...
  'Text','X',...
  'ButtonPushedFcn',@XbandButtonPushedFcn);
xpos = xpos + 35;
hQbandButton = uibutton(hFig);
set(hQbandButton,...
  'Position',[xpos border 30 22],...
  'Text','Q',...
  'ButtonPushedFcn',@QbandButtonPushedFcn);
xpos = xpos + 35;
hWbandButton = uibutton(hFig);
set(hWbandButton,...
  'Position',[xpos border 30 22],...
  'Text','W',...
  'ButtonPushedFcn',@WbandButtonPushedFcn);

% Table of isotope data
hTable = uitable(hFig);
set(hTable,...
  'Position',[xOff border+bottomHeight figPos(3)-2*border tableHeight]);
figdata.hTable = hTable;

tabledata = data(:,{'isotope','abundance','spin','gn','gamma','qm'});
tabledata.NMRfreq = zeros(height(tabledata),1);

hTable.Data = tabledata;
hTable.ColumnWidth = 'auto';
if ~verLessThan('matlab','9.11')  % 9.11 = R2021b
  hTable.SelectionType = 'row';
end
hTable.ColumnName = {'Isotope','Abundance (%)','Spin',...
  'gn value','γ/2π (MHz/T)','Q (barn)','Frequency (MHz)'};

if ~verLessThan('matlab','9.7')  % 9.7 = R2019b
  hTable.ColumnSortable = true;
end

figdata.Element = '';
figdata.tableData = tabledata;

guidata(hFig,figdata);
updateTable;

end


%-------------------------------------------------------------------------------
function elementButtonPushedCallback(src,~)
Element = src.Text;
hFig = src.Parent;
data = guidata(hFig);

if Element=="all", Element = ''; end
data.Element = Element;
guidata(hFig,data);

updateTable;

end


%-------------------------------------------------------------------------------
function updateTable()
hFig = findall(0,'Type','figure','Tag','isotopesfig');
data = guidata(hFig);
hTable = data.hTable;

element = data.Element;

% Update NMR frequencies
B0 = data.hFieldEdit.Value;
if isempty(B0)
  errordlg('Invalid magnetic field value!');
  B0 = data.DefaultField;
  data.hFieldEdit.Value = data.DefaultField;
end
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
if ~exist(DataFile,'file')
  error('Could not open nuclear isotopes data file %s',DataFile);
end

% Load data
fh = fopen(DataFile);
C = textscan(fh,'%f %f %s %s %s %f %f %f %f','commentstyle','%');

% Calculate gyromagnetic ratioes (MHz/T)
gn = C{7};
C{10} = gn*nmagn/planck/1e6;

% Assemble isotope symbols
N = C{2};
element = C{4};
radioactive = C{3};
for k = numel(N):-1:1
  isostr = sprintf('%d%s',N(k),element{k});
  if radioactive{k}=='*'
    isotopes{k} = [isostr '*'];
  else
    isotopes{k} = isostr;
  end
end
C{11} = isotopes(:);

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
function XbandButtonPushedFcn(~,~)

hFig = findall(0,'Type','figure','Tag','isotopesfig');
data = guidata(hFig);

data.hFieldEdit.Value = 340;  % mT

updateTable;

end


%-------------------------------------------------------------------------------
function QbandButtonPushedFcn(~,~)

hFig = findall(0,'Type','figure','Tag','isotopesfig');
data = guidata(hFig);

data.hFieldEdit.Value = 1200;  % mT

updateTable;

end


%-------------------------------------------------------------------------------
function WbandButtonPushedFcn(~,~)

hFig = findall(0,'Type','figure','Tag','isotopesfig');
data = guidata(hFig);

data.hFieldEdit.Value = 3400;  % mT

updateTable;

end
