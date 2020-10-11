%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_BrukerBES3T(FullBaseName,FileExtension,Scaling)
%-------------------------------------------------------------------------------
% BES3T file processing
% (Bruker EPR Standard for Spectrum Storage and Transfer)
%    .DSC: descriptor file
%    .DTA: data file
% used on Bruker ELEXSYS and EMX machines
% Code based on BES3T version 1.2 (Xepr 2.1)
%-------------------------------------------------------------------------------

if ismember(FileExtension,{'.DSC','.DTA'})
  ParExtension = '.DSC';
  SpcExtension = '.DTA';
else
  ParExtension = '.dsc';
  SpcExtension = '.dta';
end

% Read descriptor file (contains key-value pairs)
[Parameters,err] = readDSCfile([FullBaseName ParExtension]);
error(err);

% IKKF: Complex-data Flag
% CPLX indicates complex data, REAL indicates real data.
if isfield(Parameters,'IKKF')
  parts = regexp(Parameters.IKKF,',','split');
  nDataValues = numel(parts); % number of data values per parameter point
  for k = nDataValues:-1:1
    switch parts{k}
      case 'CPLX', isComplex(k) = 1;
      case 'REAL', isComplex(k) = 0;
      otherwise, error('Unknown value for keyword IKKF in .DSC file!');
    end
  end
else
  warning('Keyword IKKF not found in .DSC file! Assuming IKKF=REAL.');
  isComplex = 0;
  nDataValues = 1;
end

% XPTS: X Points   YPTS: Y Points   ZPTS: Z Points
% XPTS, YPTS, ZPTS specify the number of data points in
%  x, y and z dimension.
if isfield(Parameters,'XPTS'), nx = sscanf(Parameters.XPTS,'%f'); else, error('No XPTS in DSC file.'); end
if isfield(Parameters,'YPTS'), ny = sscanf(Parameters.YPTS,'%f'); else, ny = 1; end
if isfield(Parameters,'ZPTS'), nz = sscanf(Parameters.ZPTS,'%f'); else, nz = 1; end
Dimensions = [nx,ny,nz];

% BSEQ: Byte Sequence
% BSEQ describes the byte order of the data. BIG means big-endian,
% LIT means little-endian. Sun and Motorola-based systems are
% big-endian (MSB first), Intel-based system little-endian (LSB first).
if isfield(Parameters,'BSEQ')
  switch Parameters.BSEQ
    case 'BIG', ByteOrder = 'ieee-be';
    case 'LIT', ByteOrder = 'ieee-le';
    otherwise, error('Unknown value for keyword BSEQ in .DSC file!');
  end
else
  warning('Keyword BSEQ not found in .DSC file! Assuming BSEQ=BIG.');
  ByteOrder = 'ieee-be';
end

% IRFMT: Item Real Format
% IIFMT: Item Imaginary Format
% Data format tag of BES3T is IRFMT for the real part and IIFMT
% for the imaginary part.
if isfield(Parameters,'IRFMT')
  parts = regexp(Parameters.IRFMT,',','split');
  if numel(parts)~=nDataValues
    error('Problem in BES3T DSC file: inconsistent IKKF and IRFMT fields.');
  end
  for k = 1:nDataValues
    switch upper(parts{k})
      case 'C', NumberFormat = 'int8';
      case 'S', NumberFormat = 'int16';
      case 'I', NumberFormat = 'int32';
      case 'F', NumberFormat = 'float32';
      case 'D', NumberFormat = 'float64';
      case 'A', error('Cannot read BES3T data in ASCII format!');
      case {'0','N'}, error('No BES3T data!');
      otherwise
        error('Unknown value for keyword IRFMT in .DSC file!');
    end
  end
else
  error('Keyword IRFMT not found in .DSC file!');
end

% We enforce IRFMT and IIFMT to be identical.
if isfield(Parameters,'IIFMT')
  if any(upper(Parameters.IIFMT)~=upper(Parameters.IRFMT))
    error('IRFMT and IIFMT in DSC file must be identical.');
  end
end

% Construct abscissa vectors
AxisNames = {'X','Y','Z'};
for a = 3:-1:1
  if Dimensions(a)<=1, continue; end
  AxisType = Parameters.([AxisNames{a} 'TYP']);
  if strcmp(AxisType,'IGD')
    % Nonlinear axis -> Try to read companion file (.XGF, .YGF, .ZGF)
    companionFileName = [FullBaseName '.' AxisNames{a} 'GF'];
    % Determine data format form XFMT/YMFT/ZFMT
    DataFormat = Parameters.([AxisNames{a} 'FMT']);
    switch DataFormat
      case 'D', sourceFormat = 'float64';
      case 'F', sourceFormat = 'float32';
      case 'I', sourceFormat = 'int32';
      case 'S', sourceFormat = 'int16';
      otherwise
        error('Cannot read data format %s for companion file %s',DataFormat,companionFileName);
    end
    % Open and read companion file
    fg = fopen(companionFileName,'r',ByteOrder);
    if fg>0
      Abscissa{a} = fread(fg,Dimensions(a),sourceFormat,ByteOrder);
      fclose(fg);
    else
      warning('Could not read companion file %s for nonlinear axis. Assuming linear axis.',companionFileName);
      AxisType = 'IDX';
    end
  end
  if strcmp(AxisType,'IDX')
    Minimum(a) = sscanf(Parameters.([AxisNames{a} 'MIN']),'%f');
    Width(a) = sscanf(Parameters.([AxisNames{a} 'WID']),'%f');
    if (Width(a)==0)
      fprintf('Warning: %s range has zero width.\n',AxisNames{a});
      Minimum(a) = 1;
      Width(a) = Dimensions(a)-1;
    end
    Abscissa{a} = Minimum(a) + linspace(0,Width(a),Dimensions(a));
  end
  if strcmp(AxisType,'NTUP')
    error('Cannot read data with NTUP axes.');
  end
end
if (numel(Abscissa)==1)
  Abscissa = Abscissa{1}(:);
end

% Read data matrix.
Data = getmatrix([FullBaseName,SpcExtension],Dimensions,NumberFormat,ByteOrder,isComplex);

% Scale spectrum/spectra
if ~isempty(Scaling)
  
  % #SPL/EXPT: type of experiment
  cwExperiment = strcmp(Parameters.EXPT,'CW');
  
  % #DSL/signalChannel/SctNorm: indicates whether CW data are already scaled
  if ~isfield(Parameters,'SctNorm')
    %error('Missing SctNorm field in the DSC file. Cannot determine whether data is already scaled')
    DataPreScaled = 0;
  else
    DataPreScaled = strcmpi(Parameters.SctNorm,'true');
  end
  
  % Number of scans
  if any(Scaling=='n')
    % #SPL/AVGS: number of averages
    if ~isfield(Parameters,'AVGS')
      error('Missing AVGS field in the DSC file.')
    end
    nAverages = sscanf(Parameters.AVGS,'%d');
    if DataPreScaled
      error('Scaling by number of scans not possible,\nsince data in DSC/DTA are already averaged\nover %d scans.',nAverages);
    else
      Data = Data/nAverages;
    end
  end
  
  % Receiver gain
  if cwExperiment
    if any(Scaling=='G')
      % #SPL/RCAG: receiver gain in decibels
      if ~isfield(Parameters,'RCAG')
        error('Cannot scale by receiver gain, since RCAG in the DSC file is missing.');
      end
      ReceiverGaindB = sscanf(Parameters.RCAG,'%f');
      ReceiverGain = 10^(ReceiverGaindB/20);
      % Xenon (according to Feb 2011 manual) uses 20*10^(RCAG/20)
      Data = Data/ReceiverGain;
    end
  end
  
  % Conversion/sampling time
  if cwExperiment && any(Scaling=='c')
    %if ~DataPreScaled
    % #SPL/SPTP: sampling time in seconds
    if ~isfield(Parameters,'SPTP')
      error('Cannot scale by sampling time, since SPTP in the DSC file is missing.');
    end
    % Xenon (according to Feb 2011 manual) already scaled data by ConvTime if
    % normalization is specified (SctNorm=True). Question: which units are used?
    % Xepr (2.6b.2) scales by conversion time even if data normalization is
    % switched off!
    ConversionTime = sscanf(Parameters.SPTP,'%f'); % in seconds
    ConversionTime = ConversionTime*1000; % s -> ms
    Data = Data/ConversionTime;
    %else
    %error('Scaling by conversion time not possible,\nsince data in DSC/DTA are already scaled.');
    %end
  end
  
  % Microwave power
  if cwExperiment
    if any(Scaling=='P')
      % #SPL/MWPW: microwave power in watt
      if ~isfield(Parameters,'MWPW')
        error('Cannot scale by power, since MWPW is absent in parameter file.');
      end
      mwPower = sscanf(Parameters.MWPW,'%f')*1000; % in milliwatt
      Data = Data/sqrt(mwPower);
    end
  else
    if any(Scaling=='P')
      error('Cannot scale by microwave power, since these are not CW-EPR data.');
    end
  end
  
  % Temperature
  if any(Scaling=='T')
    % #SPL/STMP: temperature in kelvin
    if ~isfield(Parameters,'STMP')
      error('Cannot scale by temperature, since STMP in the DSC file is missing.');
    end
    Temperature = sscanf(Parameters.STMP,'%f');
    Data = Data*Temperature;
  end
  
end

Parameters = parsefieldparams(Parameters);
return
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
function [Parameters,err] = readDSCfile(DSCFileName)

Parameters = [];
err = [];

if exist(DSCFileName,'file')
  fh = fopen(DSCFileName);
  MatlabVersion = sscanf(version,'%f',1);
  if MatlabVersion<8.6  % prior to R2015b
    warning('off','MATLAB:textscan:BufSizeDeprecation');
    allLines = textscan(fh,'%s','whitespace','','delimiter','\n','bufsize',500000); %#ok<BUFSIZE>
  else
    allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
  end
  allLines = allLines{1};
  fclose(fh);
else
  err = sprintf('Cannot find the file %s.',DSCFileName);
  return
end

for k=1:numel(allLines)
  
  line = allLines{k};
  
  % Go to next if line is empty
  if isempty(line); continue; end
  
  % If line is terminated by \, append next line
  if (line(end)=='\')
    k2 = k+1;
    while (allLines{k2}(end)=='\')
      line = [line(1:end-1) allLines{k2}];
      allLines{k2} = '';
      k2 = k2 + 1;
    end
    line(end) = '';
    % Replace all \n with newline character
    line = sprintf(line);
  end
  
  [Key,Value] = strtok(line);
  if isempty(Key); continue; end
  
  % If key is not valid, go to next line.
  if ~isletter(Key(1))
    % Stop reading when Manipulation History Layer is reached.
    if strcmpi(Key,'#MHL'); break; end
    continue;
  end
  
  Value = deblank(Value(end:-1:1)); Value = deblank(Value(end:-1:1));
  
  if ~isempty(Value)
    if Value([1 end])==''''
      Value([1 end]) = [];
    end
  end
  if ~isvarname(Key)
    Key = genvarname(Key);
  end
  % Set field in output structure.
  Parameters.(Key) = Value;
  
end

return
