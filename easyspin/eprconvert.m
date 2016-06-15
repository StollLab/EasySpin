% eprconvert Frequency/field/g value conversion utility 
%
%   eprconvert
%
%   A dialog for conversion of frequency and field units and
%   for computation of one of frequeny/field/g value from the
%   two other quantities using the equation
%
%      planck*frequency = bmagn*field*g

function varargout = eprconvert(varargin)

% Last Modified by GUIDE v2.5 20-Jun-2004 17:32:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eprconvert_OpeningFcn, ...
                   'gui_OutputFcn',  @eprconvert_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = esguimain(gui_State, varargin{:});
else
    esguimain(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before eprconvert is made visible.
function eprconvert_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eprconvert (see VARARGIN)

% Choose default command line output for eprconvert
handles.output = hObject;

Col = get(0,'defaultUicontrolBackgroundColor');
set(handles.FreqList,'BackgroundColor',Col);
set(handles.FieldList,'BackgroundColor',Col);
set(hObject,'Color',get(0,'defaultUicontrolBackgroundColor'));
set(hObject,'WindowStyle','normal','Resize','off');

handles.FieldSI = [];
handles.FreqSI = 9.8e9;
handles.gVal = gfree;
UpdateDisplay(handles);

% Update handles structure
guidata(handles.eprconvert, handles);

% UIWAIT makes eprconvert wait for user response (see UIRESUME)
% uiwait(handles.eprconvert);


% --- Outputs from this function are returned to the command line.
function varargout = eprconvert_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function FreqEdit_CreateFcn(hObject, eventdata, handles)

function FreqEdit_Callback(hObject, eventdata, handles)
handles.FreqSI = GetFreqInSIUnits(handles);
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

function FreqUnits_CreateFcn(hObject, eventdata, handles)

function FreqUnits_Callback(hObject, eventdata, handles)
handles.FreqSI = GetFreqInSIUnits(handles);
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

function gEdit_CreateFcn(hObject, eventdata, handles)

function gEdit_Callback(hObject, eventdata, handles)
handles.gVal = GetgValue(handles);
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

function gfreeButton_Callback(hObject, eventdata, handles)
handles.gVal = gfree;
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

function FieldEdit_CreateFcn(hObject, eventdata, handles)

function FieldEdit_Callback(hObject, eventdata, handles)
handles.FieldSI = GetFieldInSIUnits(handles);
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

function FieldUnits_CreateFcn(hObject, eventdata, handles)

function FieldUnits_Callback(hObject, eventdata, handles)
handles.FieldSI = GetFieldInSIUnits(handles);
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

%-----------------------------------------------------------------
function ComputeFieldButton_Callback(hObject, eventdata, handles)
if isempty(handles.FreqSI),
  errordlg('Invalid frequency value!');
  return;
end
if isempty(handles.gVal),
  errordlg('Invalid g value!');
  return;
end
handles.FieldSI = planck*handles.FreqSI/bmagn/handles.gVal;
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

function ComputegButton_Callback(hObject, eventdata, handles)
if isempty(handles.FreqSI),
  errordlg('Invalid frequency value!');
  return;
end
if isempty(handles.FieldSI),
  errordlg('Invalid field value!');
  return;
end
handles.gVal = planck*handles.FreqSI/bmagn/handles.FieldSI;
guidata(handles.eprconvert,handles);
UpdateDisplay(handles);

function ComputeFreqButton_Callback(hObject, eventdata, handles)
if isempty(handles.FieldSI),
  errordlg('Invalid field value!');
  return;
end
if isempty(handles.gVal),
  errordlg('Invalid g value!');
  return;
end
handles.FreqSI = bmagn*handles.gVal*handles.FieldSI/planck;
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

function FieldList_CreateFcn(hObject, eventdata, handles)
function FieldList_Callback(hObject, eventdata, handles)
function FreqList_CreateFcn(hObject, eventdata, handles)
function FreqList_Callback(hObject, eventdata, handles)

function [Unit,prefactor] = GetFreqUnit(handles)
UnitStr = get(handles.FreqUnits,'String');
iUnit = get(handles.FreqUnits,'Value');
Unit = UnitStr{iUnit};
switch Unit
  case 'MHz', prefactor = 1e6;
  case 'GHz', prefactor = 1e9;
  case 'cm^-1', prefactor = 100*clight;
  case 'eV', prefactor = planck/echarge;
end
return

function [Unit,prefactor] = GetFieldUnit(handles)
UnitStr = get(handles.FieldUnits,'String');
iUnit = get(handles.FieldUnits,'Value');
Unit = UnitStr{iUnit};
switch Unit
  case 'mT', prefactor = 1e-3;
  case 'T', prefactor = 1;
  case 'G', prefactor = 1e-4;
  case 'kG', prefactor = 1e-1;
end
return

function Freq = GetFreqInSIUnits(handles)
Freq = str2num(get(handles.FreqEdit,'String'));
[Unit,prefactor] = GetFreqUnit(handles);
Freq = prefactor*Freq;
return

function Field = GetFieldInSIUnits(handles)
Field = str2num(get(handles.FieldEdit,'String'));
[Unit,prefactor] = GetFieldUnit(handles);
Field = prefactor*Field;
return

function g = GetgValue(handles)
g = str2num(get(handles.gEdit,'String'));
return

function SetgVal(g,handles)
set(handles.gEdit,'String',sprintf('%g',g));
return

%---------------------------------------------------------------
function UpdateDisplay(handles)

FreqSI = handles.FreqSI;
[Unit,prefactor] = GetFreqUnit(handles);
set(handles.FreqEdit,'String',sprintf('%g',FreqSI/prefactor));
if ~isempty(FreqSI)
  FrqStr{1} = sprintf('%g MHz',FreqSI/1e6);
  FrqStr{2} = sprintf('%g GHz',FreqSI/1e9);
  FrqStr{3} = sprintf('%g cm^-1',FreqSI/100/clight);
  FrqStr{4} = sprintf('%g eV',FreqSI*planck/echarge);
else
  FrqStr = '';
end
set(handles.FreqList,'String',FrqStr,'FontSize',10);

FieldSI = handles.FieldSI;
[Unit,prefactor] = GetFieldUnit(handles);
FieldString = sprintf('%g',FieldSI/prefactor);
set(handles.FieldEdit,'String',FieldString);
if ~isempty(FieldSI)
  FieldStr{1} = sprintf('%g mT',FieldSI*1e3);
  FieldStr{2} = sprintf('%g T',FieldSI);
  FieldStr{3} = sprintf('%g G',FieldSI*1e4);
  FieldStr{4} = sprintf('%g kG',FieldSI*10);
else
  FieldStr = '';
end
set(handles.FieldList,'String',FieldStr,'FontSize',10);

g = handles.gVal;
set(handles.gEdit,'String',sprintf('%g',g));


% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)

handles.FreqSI = [];
handles.FieldSI = [];
handles.gVal = [];
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);


% --- Executes on button press in WbandButton.
function WbandButton_Callback(hObject, eventdata, handles)
% hObject    handle to WbandButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.FreqSI = 94e9;
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

% --- Executes on button press in QbandButton.
function QbandButton_Callback(hObject, eventdata, handles)
% hObject    handle to QbandButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.FreqSI = 34e9;
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

% --- Executes on button press in XbandButton.
function XbandButton_Callback(hObject, eventdata, handles)
% hObject    handle to XbandButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.FreqSI = 9.8e9;
UpdateDisplay(handles);
guidata(handles.eprconvert,handles);

