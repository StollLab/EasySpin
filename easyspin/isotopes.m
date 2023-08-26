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

% Read data
data = readDataFile;
figdata.Data = data;
figdata.DefaultField = Field;

% Initialize figure window
hFig = uifigure();
set(hFig,...
  'Tag',figTag,...
  'Name','Nuclear isotopes',...
  'Toolbar','none',...
  'Menubar','none',...
  'NumberTitle','off',...
  'Resize','off');

% GUI dimensions (pixels)
elementWidth = 36;
elementHeight = elementWidth;
border = 10;
spacing = 5;
xSpacing = elementWidth+spacing;
ySpacing = elementHeight+spacing;
classSpacing = 5;
labelHeight = 15;
listHeight = 150;
bottomHeight = 30;

% Calculate and set window size
windowWidth = border + 18*xSpacing + 2*classSpacing + border;
windowHeight = border + 7*ySpacing + border + 2*ySpacing + border + ...
  + labelHeight/2 + listHeight + border + bottomHeight;
figPos = hFig.Position;
figPos(3) = windowWidth;
figPos(4) = windowHeight;
hFig.Position = figPos;

% Add element buttons
ordNumber = 0;
yOff = figPos(4)-elementHeight-border;
xOff = border;
for k = 1:numel(data.gn)
  if data.Protons(k)<=ordNumber, continue; end
  ordNumber = data.Protons(k);
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
    'Text',data.Element{k},...
    'FontSize',buttonFontSize,...
    'Tooltip',[' ' data.Name{k} ' ']);
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
  if data.Nucleons(k)<=0
    bgcol = get(hButton,'BackgroundColor');
  end  
  hButton.BackgroundColor = bgcol;
end

hAll = uibutton(hFig);
p = [xOff+16*xSpacing+2*classSpacing yOff-8*ySpacing-classSpacing ...
     elementWidth+xSpacing elementHeight+ySpacing];
set(hAll,...
  'Position',p,...
  'Text','all',...
  'BackgroundColor',[1 1 1]*0.9,...
  'Tooltip','all elements',...
  'ButtonPushedFcn',@ElementCallback,...
  'FontSize',buttonFontSize);

% Add labels above data columns
%{
labels = ['Mass Sym.    I              Nuclear g             '...
    'ENDOR Freq.          Elec. Quadrupole      Abundance']; 
hListLabel = uicontrol('Style','text','String',labels);
set(hListLabel,'HorizontalAlignment','left');
set(hListLabel,'Position',[xOff Border+bottomHeight+listHeight figPos(3)-2*Border labelHeight]);
figdata.hListLabel = hListLabel;
%}

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
  'Position',[xOff border+bottomHeight figPos(3)-2*border listHeight],...
  'FontName',get(0,'FixedWidthFontName'),...
  'BackgroundColor',[1 1 1]);
figdata.hTable = hTable;

figdata.Element = 'all';

guidata(hFig,figdata);
%updateList;

end


%-------------------------------------------------------------------------------
function elementButtonPushedCallback(src,~)
Element = src.Text;
hFig = src.Parent;
data = guidata(hFig);
data.Element = Element;
guidata(hFig,data);

updateList;

end


%-------------------------------------------------------------------------------
function updateList()
hFig = findall(0,'Type','figure','Tag','isotopesfig');
data = guidata(hFig);

newField = data.hFieldEdit.Value;
if isempty(newField)
  errordlg('Invalid magnetic field value!');
  newField = data.DefaultField;
  data.hFieldEdit.Value = data.DefaultField;
end

element = data.Element;
if strcmp(element,'all')
  element = '';
end

lines = IsotopeTable(element,data.Data,newField);
set(data.hTable,'Value',1);
set(data.hTable,'String',lines);

end


%-------------------------------------------------------------------------------
function Lines = IsotopeTable(Element,Data,Field)
Lines = [];
if isempty(Element)
  Isotopes = 1:numel(Data.gn);
else
  Isotopes = find(strcmp(Element,Data.Element));
end
for k = 1:numel(Isotopes)
  iIso = Isotopes(k);
  if Data.Spin(iIso)<=0, continue; end
  if Data.Nucleons(iIso)>0
    Symbol = sprintf('%3d %s   ',Data.Nucleons(iIso),Data.Element{iIso});
  else
    Symbol = sprintf('    %s   ',Data.Element{iIso});
  end
  Symbol = Symbol(1:6);
  SpinStr = sprintf('%1.1f',Data.Spin(iIso));
  Freq = nmagn*Field*abs(Data.gn(iIso))/planck/1e9;
  Abundance = '       ';
  if Data.Radioactive{iIso}=='-'
    a = Data.Abundance(iIso);
    Abund = sprintf('%g%%',a);
    if (a<100), Abund = [' ' Abund]; end
    if (a<10), Abund = [' ' Abund]; end
    Abundance(1:numel(Abund)) = Abund;
  end
  QuadMoment = '         ';
  if Data.qm(iIso)~=0
    if isnan(Data.qm(iIso))
      d = ' n/a';
    else
      d = sprintf('%+gb',Data.qm(iIso));
    end
    QuadMoment(1:numel(d)) = d;
  end
  if Data.Spin(iIso)>0
    Format = '%s   %s   %-+10g  %-8gMHz   %s    %s';
    Lines{end+1} = sprintf(Format,Symbol,SpinStr,...
                  Data.gn(iIso),Freq,QuadMoment,Abundance);
  else
    Format = '%s                           %s';
    Lines{end+1} = sprintf(Format,Symbol,Abundance);
  end
end
end


%-------------------------------------------------------------------------------
function data = readDataFile

% Determine full data file name
esPath = fileparts(which(mfilename));
DataFile = [esPath filesep 'private' filesep 'isotopedata.txt'];
if ~exist(DataFile,'file')
  error('Could not open nuclear isotopes data file %s',DataFile);
end

% Load data
fh = fopen(DataFile);
C = textscan(fh,'%f %f %s %s %s %f %f %f %f','commentstyle','%');
[data.Protons,data.Nucleons,data.Radioactive,...
 data.Element,data.Name,data.Spin,data.gn,data.Abundance,data.qm] = C{:};

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
