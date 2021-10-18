%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_MagnettechXML(FileName)
%-------------------------------------------------------------------------------
%   XML file format of newer Magnettech spectrometers (MS5000)
%-------------------------------------------------------------------------------
% Preparation for Base64 decoding: Use Java class depending on Matlab version
if exist('org.apache.commons.codec.binary.Base64','class')
  % seen on R2012b and R2017b
  base64 = org.apache.commons.codec.binary.Base64;
  oldJavaClass = false;
elseif exist('org.apache.axis.encoding.Base64','class')
  % seen on R2007b
  base64 = org.apache.axis.encoding.Base64; 
  oldJavaClass = true;
else
  error('No Java Base64 decoder available to read Magnettech XML data.');
end

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

% Add all children Param nodes from Parameters node to Parameter structure
% (if Recipe is present - it's absent for a dip sweep)
if isfield(Measurement,'Recipe')
  ParameterList = Measurement.Recipe.Parameters.Param;
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
    Name = [Name '_' Mode];
  end
  
  % Read curve data (if they are base64 encoded)
  data = thisCurve.Text;
  if ~isempty(data)
    % Check whether it is base64 compression
    if ~strcmp(thisCurve.Attributes.Compression,'Base64')
      error('Data is not Base64 encoded. Cannot read file.');
    end
    if ~oldJavaClass
      data = typecast(int8(data),'uint8'); % typecast without changing the underlying data
      bytestream_ = base64.decode(data); % decode
      bytestream_(9:9:end) = []; % remove termination zeros
    else
      bytestream_ = base64.decode(data); % decode
    end
    data = typecast(bytestream_,'double'); % typecast without changing the underlying data
  end
  Curves.(Name).data = data;
  
  % Horizontal axis (equal to time for most data types [except for Frequency[Raw],
  % ADC_24bit[Raw])
  XOffset = sscanf(thisCurve.Attributes.XOffset,'%f');
  XSlope = sscanf(thisCurve.Attributes.XSlope,'%f');
  Curves.(Name).x = XOffset + (0:numel(data)-1)*XSlope;
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

return
%-------------------------------------------------------------------------------
