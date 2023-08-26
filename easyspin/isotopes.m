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

Field = 350;  % mT
FigTag = 'isotopesfig';

hFig = uifigure();
clf

Data = ReadDataFile;
figdata.Data = Data;
figdata.DefaultField = Field;
guidata(hFig,figdata);

set(hFig,...
  'Tag',FigTag,...
  'Toolbar','none',...
  'Menubar','none',...
  'NumberTitle','off',...
  'Resize','off',...
  'Name','Nuclear isotopes');

% Dimensions
elementWidth = 40;
elementHeight = elementWidth;
Border = 10;
xSpacing = elementWidth+3;
ySpacing = elementHeight+3;
ClassSpacing = 5;
LabelHeight = 15;
ListHeight = 100;
BottomHeight = 30;

windowWidth = Border + 18*xSpacing + 2*ClassSpacing + Border;
windowHeight = Border + 7*ySpacing + Border + 2*ySpacing + Border + ...
  + LabelHeight/2 + ListHeight + Border + BottomHeight;
figPos = hFig.Position;
figPos(3) = windowWidth;
figPos(4) = windowHeight;
hFig.Position = figPos;

yOff = figPos(4)-elementHeight-Border;
xOff = Border;

% Add element buttons
OrdNumber = 0;
for k = 1:numel(Data.gn)
  if Data.Protons(k)<=OrdNumber, continue; end
  OrdNumber = Data.Protons(k);
  p = [0 0 elementWidth elementHeight];
  [period,group,cl] = elementclass(OrdNumber);
  if cl==2
    p(2) = yOff - (period+1)*ySpacing - ClassSpacing;
    p(1) = xOff + (group-1)*xSpacing + ClassSpacing;
  else
    p(2) = yOff - (period-1)*ySpacing;
    p(1) = xOff + (group-1)*xSpacing;
    if group>2, p(1) = p(1) + ClassSpacing; end
    if group>12, p(1) = p(1) + ClassSpacing; end
  end
  hButton = uibutton(hFig);
  set(hButton,'Position',p,...
    'ButtonPushedFcn',@ElementCallback,...
    'Text',Data.Element{k},...
    'FontSize',14,...
    'Tooltip',[' ' Data.Name{k} ' ']);
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
  if Data.Nucleons(k)<=0
    bgcol = get(hButton,'BackgroundColor');
  end  
  hButton.BackgroundColor = bgcol;
end

hAll = uicontrol('Style','pushbutton','Units','pixels',...
'Position',[xOff+16*xSpacing+2*ClassSpacing yOff-8*ySpacing-ClassSpacing elementWidth+xSpacing elementHeight+ySpacing]);
set(hAll,'String','all','BackgroundColor',[1 1 1]*0.9,...
'ToolTipstring','all elements',...
'Callback','ElementCallback','FontSize',10);

% Add labels above data columns
labels = ['Mass Sym.    I              Nuclear g             '...
    'ENDOR Freq.          Elec. Quadrupole      Abundance']; 
hListLabel = uicontrol('Style','text','String',labels);
set(hListLabel,'HorizontalAlignment','left');
set(hListLabel,'Position',[xOff Border+BottomHeight+ListHeight figPos(3)-2*Border LabelHeight]);
figdata.hListLabel = hListLabel;

% Magnetic field value
uicontrol('Style','text','String','Magnetic field (mT)',...
'Position',[xOff+80 Border 150 19],'HorizontalAlignment','left');

hFieldEdit = uieditfield(hFig,'numeric');
set(hFieldEdit,...
  'BackgroundColor','white',...
  'Position',[xOff+200 Border 100 22],...
  'Value',Field,...
  'HorizontalAlignment','left',...
  'ValueChangedFcn',@updateList);
figdata.hFieldEdit = hFieldEdit;

% Sort selection
hSort = uicontrol('Style','popupmenu','String',...
{'nucleons','gn value','spin','abundance'},'Position',[xOff Border 100 22],...
'HorizontalAlignment','left','Callback',@resort);
set(hSort,'Visible','off','Enable','inactive');
figdata.hSort = hSort;

% List of isotope data
hList = uicontrol('Style','listbox','HorizontalAlignment','left');
set(hList,'FontName',get(0,'FixedWidthFontName'));
set(hList,'Units','pixels','Position',[xOff Border+BottomHeight figPos(3)-2*Border ListHeight]);
set(hList,'BackgroundColor','white');
figdata.hList = hList;

figdata.Element = 'all';

guidata(hFig,figdata);
updateList;

set(hFig,'HandleVisibility','callback');

end

%-------------------------------------------------------------------------------
function resort
hFig = findobj('Tag','isotopesfig');
d = guidata(hFig);
v = get(d.hSort,'Value');
switch v
case 1, [~,idx] = sort(d.Data.Nucleons);
case 2, [~,idx] = sort(abs(d.Data.gn),'descend');
case 3, [~,idx] = sort(d.Data.Spin,'descend');
case 4, [~,idx] = sort(d.Data.Abundance,'descend');
end

d.Data.Protons = d.Data.Protons(idx);
d.Data.Nucleons = d.Data.Nucleons(idx);
d.Data.Radioactive = d.Data.Radioactive(idx);
d.Data.Element = d.Data.Element(idx);
d.Data.Name = d.Data.Name(idx);
d.Data.Spin = d.Data.Spin(idx);
d.Data.gn = d.Data.gn(idx);
d.Data.Abundance = d.Data.Abundance(idx);
guidata(hFig,d);
updateList;

end

%-------------------------------------------------------------------------------
function ElementCallback(~,~)

Element = get(gcbo,'Text');
hFig = findobj('Tag','isotopesfig');
d = guidata(hFig);
d.Element = Element;
guidata(hFig,d);

updateList;

end

%-------------------------------------------------------------------------------
function updateList
hFig = findobj('Tag','isotopesfig');
d = guidata(hFig);

newField = d.hFieldEdit.Value;
if isempty(newField)
  errordlg('Invalid magnetic field value!');
  newField = d.DefaultField;
  set(d.hField,'String',sprintf('%g',newField));
end

Element = d.Element;
if strcmp(Element,'all')
  Element = '';
end

Lines = IsotopeTable(Element,d.Data,newField);
set(d.hList,'Value',1);
set(d.hList,'String',Lines);

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
function Data = ReadDataFile

% Determine data file name
esPath = fileparts(which(mfilename));
DataFile = [esPath filesep 'private' filesep 'isotopedata.txt'];
%DataFile = [esPath filesep 'nucmoments.txt'];
if ~exist(DataFile,'file')
  error('Could not open nuclear data file %s',DataFile);
end

fh = fopen(DataFile);
C = textscan(fh,'%f %f %s %s %s %f %f %f %f','commentstyle','%');
[Data.Protons,Data.Nucleons,Data.Radioactive,...
 Data.Element,Data.Name,Data.Spin,Data.gn,Data.Abundance,Data.qm] = C{:};

% idx = Data.Spin<=0;
% Data.Protons(idx) = [];
% Data.Nucleons(idx) = [];
% Data.Radioactive(idx) = [];
% Data.Element(idx) = [];
% Data.Name(idx) = [];
% Data.Spin(idx) = [];
% Data.gn(idx) = [];
% Data.Abundance(idx) = [];

end


%-------------------------------------------------------------------------------
function [Period,Group,Class] = elementclass(N)

Class = 0;

PeriodLimits = [0 2 10 18 36 54 86 1000];

% Determine period of element
for Period = 1:8
  if N<=PeriodLimits(Period), break; end
end
Period = Period - 1;

%Determine group and class of element
%Class 0 - main groups, 1 - transition metals, 2 - rare earths
Group = N - PeriodLimits(Period);
switch Period
case 1
  Class = 0;
  if Group~=1, Group=18; end
case {2,3}
  Class = 0;
  if Group>2, Group = Group + 10; end
case {4,5}
  Class = 1;
  if (Group<3) || (Group>12), Class = 0; end
case {6,7}
  if (Group<3) || (Group>26)
    Class = 0;
  else
    if (Group>16)
      Class = 1;
    else
      Class = 2;
    end
  end
  if (Class<2) && (Group>16)
    Group = Group - 14;
  end
end

end
