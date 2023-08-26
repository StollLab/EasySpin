% isotopes   Graphical interface for nuclear isotope data 
%
%   Displays a PSE (periodic system of the elements) with
%   a list of nuclear isotope data. The columns in the table
%   give the following data
%
%    1. Mass number
%    2. Element symbol
%    3. Nuclear spin
%    4. gn value
%    5. ENDOR frequency in MHz, using the magnetic field
%       entered below the table
%    6. electric quadrupole moment, in barn
%    7. Natural abundance in %

function isotopes()

% Settings
figTag = 'isotopesfig';
Field = 350;  % mT
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
    bgcol = get(hButton,'BackgroundColor');
  end  
  hButton.BackgroundColor = bgcol;
end

% Add selection button
hAll = uibutton(hFig);
p = [xOff+16*xSpacing+2*classSpacing yOff-8*ySpacing-classSpacing ...
     elementWidth+xSpacing elementHeight+ySpacing];
set(hAll,...
  'Position',p,...
  'Text','all',...
  'BackgroundColor',[1 1 1]*0.9,...
  'Tooltip','all elements',...
  'ButtonPushedFcn',@elementButtonPushedCallback,...
  'FontSize',buttonFontSize);

% Magnetic field edit box with label
hFieldEdit = uieditfield(hFig,'numeric');
set(hFieldEdit,...
  'BackgroundColor','white',...
  'Position',[xOff+200 border 100 22],...
  'Value',Field,...
  'HorizontalAlignment','left',...
  'ValueChangedFcn',@updateList);
figdata.hFieldEdit = hFieldEdit;
uilabel(hFig,...
  'Text','Magnetic field (mT)',...
  'Position',[xOff+80 border 150 19]);

% Table of isotope data
hTable = uitable(hFig);
set(hTable,...
  'Position',[xOff border+bottomHeight figPos(3)-2*border tableHeight]);
figdata.hTable = hTable;

%vnames = {'Z','N','radioactive','element','name','spin',...
%  'gn','abundance','qm','gamma','isotope'};

tabledata = data(:,{'Z','N','isotope','abundance','spin','gn','gamma','qm'});

hTable.Data = tabledata;
hTable.ColumnSortable = true;
hTable.ColumnWidth = 'auto';
hTable.SelectionType = 'row';
hTable.ColumnName = {'Z','N','Isotope','Abundance (%)','Spin','gn value','γ/2π (MHz/T)','Q (barn)'};

figdata.Element = '';
figdata.tableData = tabledata;

guidata(hFig,figdata);

end


%-------------------------------------------------------------------------------
function elementButtonPushedCallback(src,~)
Element = src.Text;
hFig = src.Parent;
data = guidata(hFig);

if Element=="all", Element = ''; end
data.Element = Element;
guidata(hFig,data);

updateList;

end


%-------------------------------------------------------------------------------
function updateList()
hFig = findall(0,'Type','figure','Tag','isotopesfig');
data = guidata(hFig);
hTable = data.hTable;

element = data.Element;

B0 = data.hFieldEdit.Value;
if isempty(B0)
  errordlg('Invalid magnetic field value!');
  B0 = data.DefaultField;
  data.hFieldEdit.Value = data.DefaultField;
end

if isempty(element)
  hTable.Data = data.tableData;
else
  idx = data.fullData.element==string(element);
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
  str = sprintf('%d%s',N(k),element{k});
  if radioactive{k}=='*'
    str = [str '*'];
  end
  isotopes{k} = str;
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
