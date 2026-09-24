%-------------------------------------------------------------------------------
%   Read XML file format of Magnettech spectrometers (MS5000)
%-------------------------------------------------------------------------------

function [Data, Abscissa, Parameters] = eprload_MagnettechXML(FileName)

% Read XML file and convert to Matlab structure for easy access
Document = xml2struct(FileName);

% Assert it's an XML file from a Magnettech spectrometer
if ~isfield(Document,'ESRXmlFile')
  error('File %s is not a Magnettech xml file.',FileName);
elseif ~isfield(Document.ESRXmlFile,'Data')
  error('ESRXmlFile.Data node not found in xml file.');
elseif ~isfield(Document.ESRXmlFile.Data,'Measurement')
  error('ESRXmlFile.Data.Measurement node not found in xml file.');
elseif ~isfield(Document.ESRXmlFile.Data.Measurement,'DataCurves')
  error('ESRXmlFile.Data.Measurement.DataCurves node not found in xml file.');
elseif ~isfield(Document.ESRXmlFile.Data.Measurement.DataCurves,'Curve')
  error('No data. No <Curve> node found in xml file.');
end

% Add attributes from Measurement node to Paramater structure
Data = Document.ESRXmlFile.Data;
Measurement = Data.Measurement;
Parameters = Measurement.Attributes;
CurveList = Measurement.DataCurves.Curve;
% xml2struct returns a single element as struct, multiple ones as cell array
if isstruct(CurveList), CurveList = {CurveList}; end

% Add all children Param nodes from Parameters node to Parameter structure
% (if Recipe is present - it's absent for a dip sweep)
if isfield(Measurement,'Recipe')
  ParameterList = Measurement.Recipe.Parameters.Param;
  if isstruct(ParameterList), ParameterList = {ParameterList}; end
  for p = 1:numel(ParameterList)
    PName = ParameterList{p}.Attributes.Name;
    P_ = ParameterList{p}.Text;
    Parameters.(PName) = P_;
  end
end

% Determine what type of sweep this is (field sweep or other)
if isfield(Measurement,'Recipe')
  switch Measurement.Recipe.Attributes.Type
    case 'single', xAxis = 'field';
    case 'kinetic', xAxis = 'time';
    case 'temperature', xAxis = 'field';
    case 'goniometric', xAxis = 'field';
    case 'modulationSweep', xAxis = 'field';
    case 'powerSweep', xAxis = 'field';
    case 'xRay', xAxis = 'field';
    otherwise
      error('Unknown Recipe.Type in the xml file.');
  end
  if isfield(Parameters,'BHoldEnabled')
    if strcmp(Parameters.BHoldEnabled,'True')
      xAxis = 'time';
    end
  end
end

% Get data (stored in a series of <Curve ...> </Curve> nodes under <DataCurves>)
for iCurve = 1:numel(CurveList)
  thisCurve = CurveList{iCurve};
  Name = thisCurve.Attributes.YType;
  Mode = thisCurve.Attributes.Mode;
  
  % Avoid duplicate names (e.g. BField can be stored twice in the same file, once
  % with Mode='Raw' and once with Mode='Pre')
  if strcmp(Name,'BField') && strcmp(Mode,'Raw')
    Name = [Name '_' Mode];  %#ok
  end
  
  % Read curve data (if they are base64 encoded)
  data_txt = thisCurve.Text;
  if ~isempty(data_txt)
    % Check whether it is base64 compression
    if ~strcmp(thisCurve.Attributes.Compression,'Base64')
      error('Data is not Base64 encoded. Cannot read file.');
    end

    % Append 6 zero bits to each double (64 bits) to pad to 9 base64
    % chars per entry (9*6=72 bits=9 chars)
    data_txt(data_txt=='=') = 'A';

    % Decode (9 bytes per double)
    bytestream_ = base64decode(data_txt);
    
    bytestream_(9:9:end) = []; % remove every ninth byte (double has 8 bytes)
    data = typecast(bytestream_,'double'); % typecast without changing the underlying data
    Curves.(Name).data = data;

    % Check whether all chars got decoded and typecast.
    if numel(data)~=numel(data_txt)/12
      error('Number of base64 chars and number of data points are inconsistent.');
    end

    % Horizontal axis (equal to time for most data types [except for Frequency[Raw],
    % ADC_24bit[Raw])
    XOffset = sscanf(thisCurve.Attributes.XOffset,'%f');
    XSlope = sscanf(thisCurve.Attributes.XSlope,'%f');
    Curves.(Name).x = XOffset + (0:numel(data)-1)*XSlope;
  end
  
end

if isfield(Curves,'BField')
  % Field sweeps, transients, and parameter sweeps
  
  switch xAxis
    case 'field'
      Abscissa = interp1(Curves.BField.x,Curves.BField.data,Curves.MW_Absorption.x);
      Abscissa = Abscissa(:);
      Data = Curves.MW_Absorption.data(:);
      % Remove data outside desired field range
      xmin = str2double(Parameters.Bfrom);
      xmax = str2double(Parameters.Bto);
      idx = Abscissa>=xmin & Abscissa<=xmax;
      Abscissa = Abscissa(idx);
      Data = Data(idx);
    case 'time'
      Abscissa = Curves.MW_Absorption.x(:);
      Abscissa = Abscissa(:);
      Data = Curves.MW_Absorption.data(:);
    otherwise
      error('Cannot construct abscissa for this measurement type');
  end
  
elseif isfield(Curves,'Frequency')
  % Capture of IQ raw data (dip sweep)
  
  Abscissa = Curves.Frequency.data(:);
  Data = Curves.ADC_24bit.data(:);
  
end

Parameters = parsefieldparams(Parameters);

% Store all curves from the file in the parameters
% (incl. sin and cos MW absorption data)
Parameters.Curves = Curves;

end

%-------------------------------------------------------------------------------
function bytes = base64decode(str)
% Decode Base64 string into column vector of bytes (uint8). Characters
% outside the Base64 alphabet (e.g. whitespace) are ignored.

alphabet = ['A':'Z' 'a':'z' '0':'9' '+/'];
lut = zeros(1,256,'uint8');
lut(double(alphabet)+1) = 0:63;
valid = false(1,256);
valid(double([alphabet '='])+1) = true;

str = str(valid(min(double(str),255)+1));
if mod(numel(str),4)~=0
  error('Invalid Base64 data: number of characters is not a multiple of 4.');
end
nPad = numel(str) - numel(regexprep(str,'=+$',''));

% Convert groups of 4 characters (4x6 bits) to 3 bytes (3x8 bits)
v = reshape(lut(double(str)+1),4,[]);
b1 = bitor(bitshift(v(1,:),2),bitshift(v(2,:),-4));
b2 = bitor(bitshift(v(2,:),4),bitshift(v(3,:),-2));
b3 = bitor(bitshift(v(3,:),6),v(4,:));
bytes = reshape([b1;b2;b3],[],1);
bytes(end-nPad+1:end) = [];

end
