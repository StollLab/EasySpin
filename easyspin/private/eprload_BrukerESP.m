%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_BrukerESP(FullBaseName,FileExtension,Scaling)
%-------------------------------------------------------------------------------
% ESP data file processing
%   Bruker ECS machines
%   Bruker ESP machines
%   Bruker WinEPR, Simfonia
%-------------------------------------------------------------------------------

% Read parameter file (contains key-value pairs)
ParExtension = '.par';
SpcExtension = '.spc';
if ismember(FileExtension,{'.PAR','.SPC'})
  ParExtension = upper(ParExtension);
  SpcExtension = upper(SpcExtension);
end
[Parameters,err] = eprload_readPARfile([FullBaseName,ParExtension]);
error(err);

% FileType: flag for specific file format
% w   Windows machines, WinEPR
% c   ESP machines, cw EPR data
% p   ESP machines, pulse EPR data
FileType = 'c';

TwoD = false; % Flag for two-dimensional data
isComplex = false; % Flag for complex data
nx = 1024;
ny = 1;

% For DOS ByteOrder is ieee-le, in all other cases ieee-be
if isfield(Parameters,'DOS')
  Endian = 'ieee-le';
  FileType = 'w';
else
  Endian = 'ieee-be';
end

% Analyse data type flags stored in JSS.
if isfield(Parameters,'JSS')
  Flags = sscanf(Parameters.JSS,'%f');
  isComplex = bitget(Flags,5);
  TwoD = bitget(Flags,13);
end

% If present, SSX contains the number of x points.
if isfield(Parameters,'SSX')
  if TwoD
    if FileType=='c', FileType='p'; end
    nx = sscanf(Parameters.SSX,'%f');
    if isComplex, nx = nx/2; end
  end
end

% If present, SSY contains the number of y points.
if isfield(Parameters,'SSY')
  if TwoD
    if FileType=='c', FileType='p'; end
    ny = sscanf(Parameters.SSY,'%f');
  end
end

% If present, ANZ contains the total number of points.
if isfield(Parameters,'ANZ')
  nAnz = sscanf(Parameters.ANZ,'%f');
  if ~TwoD
    if FileType=='c', FileType='p'; end
    nx = nAnz;
    if isComplex, nx = nx/2; end
  else
    if nx*ny~=nAnz
      error('Two-dimensional data: SSX, SSY and ANZ in .par file are inconsistent.');
    end
  end
end

% If present, RES contains the number of x points.
if isfield(Parameters,'RES')
  nx = sscanf(Parameters.RES,'%f');
end
% If present, REY contains the number of y points.
if isfield(Parameters,'REY')
  ny = sscanf(Parameters.REY,'%f');
end

% If present, XPLS contains the number of x points.
if isfield(Parameters,'XPLS')
  nx = sscanf(Parameters.XPLS,'%f');
end

% Set number format
switch FileType
  case 'c', NumberFormat = 'int32';
  case 'w', NumberFormat = 'float'; % WinEPR/Simfonia: single float
  case 'p', NumberFormat = 'int32'; % old: 'float'
end

% Construct abscissa vector
if nx>1
  
  % Get experiment type
  if ~isfield(Parameters,'JEX'), Parameters.JEX = 'field-sweep'; end
  if ~isfield(Parameters,'JEY'), Parameters.JEY = ''; end
  JEX_Endor = strcmp(Parameters.JEX,'ENDOR');
  JEX_TimeSweep = strcmp(Parameters.JEX,'Time-Sweep');
  JEY_PowerSweep = strcmp(Parameters.JEY,'mw-power-sweep');
  
  % Convert values of all possible range keywords
  %-------------------------------------------------------
  GST = []; GSI = []; HCF = []; HSW = [];
  XXLB = []; XXWI = []; XYLB = []; XYWI = [];
  if isfield(Parameters,'HCF')
    HCF = sscanf(Parameters.HCF,'%f');
  end
  if isfield(Parameters,'HSW')
    HSW = sscanf(Parameters.HSW,'%f');
  end
  if isfield(Parameters,'GST')
    GST = sscanf(Parameters.GST,'%f');
  end
  if isfield(Parameters,'GSI')
    GSI = sscanf(Parameters.GSI,'%f');
  end
  
  % XXLB, XXWI, XYLB, XYWI
  % In files from the pulse S-band spectrometer at ETH,
  % both HSW/HCF and GST/GSI are absent.
  if isfield(Parameters,'XXLB')
    XXLB = sscanf(Parameters.XXLB,'%f');
  end
  if isfield(Parameters,'XXWI')
    XXWI = sscanf(Parameters.XXWI,'%f');
  end
  if isfield(Parameters,'XYLB')
    XYLB = sscanf(Parameters.XYLB,'%f');
  end
  if isfield(Parameters,'XYWI')
    XYWI = sscanf(Parameters.XYWI,'%f');
  end
  
  % Determine which abscissa range parameters to take
  %-----------------------------------------------------------
  TakeGH = 0; % 1: take GST/GSI,  2: take HCF/HSW, 3: try XXLB
  if JEX_Endor
    % Endor experiment: take GST/GSI
    TakeGH = 1;
  elseif ~isempty(XXLB) && ~isempty(XXWI) && ~isempty(XYLB) && ~isempty(XYWI)
    % EMX 2D data -> use XXLB/XXWI/XYLB/XYWI
    TakeGH = 3;
  elseif ~isempty(HCF) && ~isempty(HSW) && ~isempty(GST) && ~isempty(GSI)
    % All fields present: take GST/GSI (even if inconsistent
    % with HCF/HSW) (not sure this is correct in all cases)
    TakeGH = 1;
  elseif ~isempty(HCF) && ~isempty(HSW)
    % Only HCF and HSW given: take them
    TakeGH = 2;
  elseif ~isempty(GST) && ~isempty(GSI)
    TakeGH = 1;
  elseif isempty(GSI) && isempty(HSW)
    HSW = 50;
    TakeGH = 2;
  elseif isempty(HCF)
    TakeGH = 3;
  end
  
  % Construct abscissa vector
  %----------------------------------------------------
  if JEX_TimeSweep
    if isfield(Parameters,'RCT')
      ConversionTime = sscanf(Parameters.RCT,'%f');
    else
      ConversionTime = 1;
    end
    Abscissa = (0:nx-1)*ConversionTime/1e3;
  else
    Abscissa = [];
    if TakeGH==1
      Abscissa = GST + GSI*linspace(0,1,nx);
    elseif TakeGH==2
      Abscissa = HCF + HSW/2*linspace(-1,1,nx);
    elseif TakeGH==3
      if ~isempty(XXLB) && ~isempty(XXWI)
        if ~isempty(XYLB) && ~isempty(XYWI)
          Abscissa{1} = XXLB + linspace(0,XXWI,nx);
          Abscissa{2} = XYLB + linspace(0,XYWI,ny);
        else
          Abscissa = XXLB + linspace(0,XXWI,nx);
        end
      end
    else
      error('Could not determine abscissa range from parameter file!');
    end
  end
  
end

% Slice of 2D data, as saved by WinEPR: RES/REY refer to
% original 2D size, but JSS 2D flag is not set -> 1D data
if ~TwoD && ny>1, ny = 1; end

% Read data file.
nz = 1;
Dimensions = [nx ny nz];
Data = getmatrix([FullBaseName,SpcExtension],Dimensions,NumberFormat,Endian,isComplex);

% Convert to a row vector in the case of 1D dataset
if ny==1, Data = Data(:).'; end
  
% Scale spectrum/spectra
if ~isempty(Scaling)
  
  % Number of scans
  if any(Scaling=='n')
    if ~isfield(Parameters,'JSD')
      error('Cannot scale by number of scans, since JSD is absent in parameter file.');
    end
    nScansDone = sscanf(Parameters.JSD,'%f');
    Data = Data/nScansDone;
  end
  
  % Receiver gain
  if any(Scaling=='G')
    if ~isfield(Parameters,'RRG')
      %Parameters.RRG = '2e4'; % default value on UC Davis ECS106
      error('Cannot scale by gain, since RRG is absent in parameter file.');
    end
    ReceiverGain = sscanf(Parameters.RRG,'%f');
    Data = Data/ReceiverGain;
  end
  
  % Microwave power
  if any(Scaling=='P')
    if ~isfield(Parameters,'MP')
      error('Cannot scale by power, since MP is absent in parameter file.');
    end
    if ~JEY_PowerSweep
      mwPower = sscanf(Parameters.MP,'%f'); % in milliwatt
      Data = Data/sqrt(mwPower);
    else
      % 2D power sweep, power along second dimension
      nPowers = size(Data,2);
      dB = XYLB+linspace(0,XYWI,nPowers);
      mwPower = sscanf(Parameters.MP,'%f'); % in milliwatt
      mwPower = mwPower.*10.^(-dB/10);
      for iPower = 1:nPowers
        Data(:,iPower) = Data(:,iPower)/sqrt(mwPower(iPower));
      end
    end
  end
  
  % Temperature
  if any(Scaling=='T')
    if ~isfield(Parameters,'TE')
      error('Cannot scale by temperature, since TE is absent in parameter file.');
    end
    Temperature = sscanf(Parameters.TE,'%f'); % in kelvin
    if Temperature==0
      error('Cannot scale by temperature, since TE is zero in parameter file.');
    end
    Data = Data*Temperature;
  end
  
  % Conversion/sampling time
  if any(Scaling=='c')
    if ~isfield(Parameters,'RCT')
      error('Cannot scale by sampling time, since RCT in the .par file is missing.');
    end
    ConversionTime = sscanf(Parameters.RCT,'%f'); % in milliseconds
    Data = Data/ConversionTime;
  end
  
end

Parameters = parsefieldparams(Parameters);
return
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
function [Parameters,err] = eprload_readPARfile(PARFileName)

Parameters = [];
err = [];

if exist(PARFileName,'file')
  fh = fopen(PARFileName);
  MatlabVersion = sscanf(version,'%f',1);
  if MatlabVersion<8.6 % prior to R2015b
    warning('off','MATLAB:textscan:BufSizeDeprecation');
    allLines = textscan(fh,'%s','whitespace','','delimiter','\n','bufsize',500000); %#ok<BUFSIZE>
  else
    allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
  end
  allLines = allLines{1};
  fclose(fh);
else
  err = sprintf('Cannot find the file %s.',PARFileName);
  return
end

for k = 1:numel(allLines)
  
  line = allLines{k};
  if isempty(line), continue; end
  
  [Key,n,err_,idx] = sscanf(line,'%s',1);
  if isempty(Key); continue; end
  
  if ~isletter(Key(1)), continue; end
  
  Value = deblank(line(end:-1:idx));
  Value = deblank(Value(end:-1:1));
  if ~isempty(Value)
    if Value([1 end])=='''' % remove leading and trailing quotes
      Value([1 end]) = [];
    end
  end
  if ~isvarname(Key)
    Key = genvarname(Key);
  end
  % set field in output structure
  Parameters.(Key) = Value;
  
end

return
