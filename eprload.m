% eprload  Load experimental EPR data 
%
%   y = eprload(FileName)
%   [x,y] = eprload(FileName)
%   [x,y,Pars] = eprload(FileName)
%   [x,y,Pars,FileN] = eprload(FileName)
%   ... = eprload(FileName,Scaling)
%   ... = eprload
%
%   Read spectral data from a file specified in the string
%   'FileName' into the arrays x (abscissa) and y (ordinate).
%   The structure Pars contains entries from the parameter
%   file, if present.
%
%   All strings in the parameter structure containing numbers
%   are converted to numbers for easier use.
%
%   If FileName is a directory, a file browser is
%   displayed. If FileName is omitted, the current
%   directory is used as default. eprload returns the
%   name of the loaded file (including its path) as
%   fourth parameter FileN.
%
%   For DSC/DTA data, x contains the vector or
%   the vectors specifying the abscissa or abscissae of the
%   spectral data array, i.e. magnetic field range
%   for cw EPR, RF range for ENDOR and time delays
%   for pulse EPR. Units are those specified in
%   the parameter file. See the fields XPTS, XMIN, XWID
%   etc. in the Pars structure.
%
%   Supported formats are identified via the extension
%   in 'FileName'. Extensions:
%
%     MAGRES:          .PLT
%     BES3T:           .DTA, .DSC
%     ESP, WinEPR:     .spc, .par
%     qese, tryscore:  .eco
%     Varian:          .spk, .ref
%     ESE:             .d00, .exp
%     SpecMan:         .d01, .exp
%
%     For reading general ASCII formats, use textread(...)
%
%   'Scaling' tells eprload to scale the data
%
%      'n':   divide by number of scans
%      'P':   divide by square root of microwave power in mW
%      'G':   divide by receiver gain
%      'T':   multiply by temperature in kelvin
%      'c':   divide by conversion/sampling time in milliseconds

function varargout = eprload(FileName,Scaling)

if (nargout<0) || (nargout>4)
  error('Please provide 1, 2, 3 or 4 output arguments!');
end

if (nargin<1), FileName = pwd; end

if (nargin<2)
  Scaling = '';
end

LocationType = exist(FileName,'file');

if (LocationType==7), % a directory
  CurrDir = pwd;
  cd(FileName);
  [uiFile,uiPath] = uigetfile({'*.DTA;*.dta;*.spc','Bruker (*.dta,*.spc)';'*.spk;*.ref','Varian (*.spk,*.ref)';'*.eco','qese/tryscore (*.eco)';...
  '*.d00','ETH/WIS (*.d00)';'*.d01','SpecMan (*.d01)';'*.plt','Magres (*.plt)'},'Load EPR data file...');
  cd(CurrDir);
  if (uiFile==0),
    varargout = cell(1,nargout);
    return;
  end
  FileName = [uiPath uiFile];
end

% Initialize output arguments
Data = [];
Parameters = [];
Abscissa = [];

% General remarks
%-----------------------------------------------------------
% No complete format specification for any format was available.
% Code for all formats except BES3T is a Matlab translation
% of cvt, a c program in use at ETH. Code for BES3T is
% built according to the "specification" of this format
% in the help file sharedHelp/BES3T.voc of Xepr Version 2.2b.
%-----------------------------------------------------------

% Decompose file name, supply default extension .DTA
[p,Name,Extension] = fileparts(FileName);
FullBaseName = fullfile(p,Name);

if isempty(Extension)
  if exist([FullBaseName '.dta'],'file'), Extension = '.dta'; end
  if exist([FullBaseName '.DTA'],'file'), Extension = '.DTA'; end
  if exist([FullBaseName '.spc'],'file'), Extension = '.spc'; end
end

FileName = [FullBaseName Extension];
LocationType = exist(FileName,'file');
if any(LocationType==[0 1 5 8]), % not a file/directory
  error('The file or directory %s does not exist!',FileName);
end

% Scaling works only for par/spc files (ECS106, ESP, etc)
if ~isempty(Scaling)
  S_ = Scaling;
  S_(S_=='n' | S_=='P' | S_=='G' | S_=='T' | S_=='c') = [];
  if ~isempty(S_)
    error('Scaling can only contain ''n'', ''P'', ''G'', ''T'', and ''c''.');
  end
end

% Determination of input file type
ParseParameters = 0;

switch Extension

case {'.DTA','.DSC','.dsc','.dta'}
  %--------------------------------------------------------
  % BES3T file processing
  % (Bruker EPR Standard for Spectrum Storage and Transfer)
  %    .DSC: descriptor file
  %    .DTA: data file
  % used on Bruker ELEXSYS and EMX machines
  % Code based on BES3T version 1.2 (Xepr 2.1)
  %--------------------------------------------------------
  ParseParameters = 1;

  if ismember(Extension,{'.DSC','.DTA'})
    ParExtension = '.DSC';
    SpcExtension = '.DTA';
  else
    ParExtension = '.dsc';
    SpcExtension = '.dta';
  end
  
  % Read descriptor file (contains key-value pairs)
  [Parameters,err] = readDSCfile([FullBaseName ParExtension]);
  error(err);
  
  % IKKF: Item Complex Flag
  % CPLX indicates complex data, REAL indicates real data.
  if isfield(Parameters,'IKKF')
    switch Parameters.IKKF
    case 'CPLX', isComplex = 1;
    case 'REAL', isComplex = 0;
    otherwise, error('Unknown value for keyword IKKF in .DSC file!');
    end
  else
    warning('Keyword IKKF not found in .DSC file! Assuming IKKF=REAL.');
    isComplex = 0;
  end
  
  % XPTS: X Points   YPTS: Y Points   ZPTS: Z Points
  % XPTS, YPTS, ZPTS specify the number of data points in
  %  x, y and z dimension.
  if isfield(Parameters,'XPTS'), nx = sscanf(Parameters.XPTS,'%f'); else error('No XPTS in DSC file.'); end
  if isfield(Parameters,'YPTS'), ny = sscanf(Parameters.YPTS,'%f'); else ny = 1; end
  if isfield(Parameters,'ZPTS'), nz = sscanf(Parameters.ZPTS,'%f'); else nz = 1; end
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
    switch upper(Parameters.IRFMT)
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
  for a=3:-1:1
    if (Dimensions(a)<=1), continue; end
    AxisType = Parameters.([AxisNames{a} 'TYP']);
    if strcmp(AxisType,'IGD')
      % Nonlinear axis -> Try to read companion file (.XGF, .YGF, .ZGF)
      fg = fopen([FullBaseName '.' AxisNames{a} 'GF'],'r',ByteOrder);
      if fg>0
        % Here we should check for the number format in
        % XFMT/YFMT/ZFMT instead of assuming 'float64'.
        Abscissa{a} = fread(fg,Dimensions(a),'float64',ByteOrder);
        fclose(fg);
      else
        warning('Could not read companion file for nonlinear axis.');
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
  Data = getmatrix([FullBaseName,SpcExtension],Dimensions,1:3,NumberFormat,ByteOrder,isComplex);

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
        error(sprintf('Scaling by number of scans not possible,\nsince data in DSC/DTA are already averaged\nover %d scans.',nAverages));
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
        ReceiverGain = 10^(ReceiverGaindB/10);
        % Xenon (according to Feb 2011 manual) uses 20*10^(RCAG/20)
        % ReceiverGain = 20*10^(ReceiverGaindB/20);
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


case {'.PAR','.SPC','.par','.spc'}
  %--------------------------------------------------
  % ESP data file processing
  %   Bruker ECS machines
  %   Bruker ESP machines
  %   Bruker WinEPR, Simfonia
  %--------------------------------------------------
  ParseParameters = 1;
  
  % Read parameter file (contains key-value pairs)
  ParExtension = '.par';
  SpcExtension = '.spc';
  if ismember(Extension,{'.PAR','.SPC'})
    ParExtension = upper(ParExtension);
    SpcExtension = upper(SpcExtension);
  end
  [Parameters,err] = readPARfile([FullBaseName,ParExtension]);
  error(err);
  
  % FileType: flag for specific file format
  % w   Windows machines, WinEPR
  % c   ESP machines, cw EPR data
  % p   ESP machines, pulse EPR data
  FileType = 'c';
  
  TwoD = 0; % Flag for two-dimensional data
  isComplex = 0; % Flag for complex data
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
    if TwoD,
      if FileType=='c', FileType='p'; end
      nx = sscanf(Parameters.SSX,'%f');
      if isComplex, nx = nx/2; end
    end
  end
  
  % If present, SSY contains the number of y points.
  if isfield(Parameters,'SSY')
    if TwoD,
      if FileType=='c', FileType='p'; end
      ny = sscanf(Parameters.SSY,'%f');
    end
  end
  
  % If present, ANZ contains the total number of points.
  if isfield(Parameters,'ANZ')
    nAnz = sscanf(Parameters.ANZ,'%f');
    if ~TwoD,
      if FileType=='c', FileType='p'; end
      nx = nAnz;
      if isComplex, nx = nx/2; end
    else
      if (nx*ny~=nAnz)
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
  if (nx>1)

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
      if (TakeGH==1)
        Abscissa = GST + GSI*linspace(0,1,nx);
      elseif (TakeGH==2)
        Abscissa = HCF + HSW/2*linspace(-1,1,nx);
      elseif (TakeGH==3)
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
  if ~TwoD && (ny>1), ny = 1; end
  
  % Read data file.
  nz = 1;
  Dimensions = [nx ny nz];
  Data = getmatrix([FullBaseName,SpcExtension],Dimensions,1:3,NumberFormat,Endian,isComplex);

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
        Parameters.RRG = '2e4'; % default value on ECS106
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
      Temperature = sscanf(Parameters.TE,'%f'); % in Kelvin
      if (Temperature==0)
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
  
case {'.ECO','.eco'}
  %--------------------------------------------------
  % ECO file processing
  %   qese     old ETH acquisition software
  %   tryscore Weizmann HYSCORE simulation program
  %--------------------------------------------------

  % open file
  fid = fopen(FileName,'r');
  if fid<0, error(['Could not open ' FileName]); end
  
  % read first line: nx ny Complex
  Data = sscanf(fgetl(fid),'%i%i%i',3)';
  
  % set dimensions and complex flag
  switch length(Data)
  case 3, Dims = Data([1 2]); isComplex = Data(3);
  case 2, Dims = Data; isComplex = 0;
  case 1, Dims = [Data 1]; isComplex = 0;
  end
  
  % read data
  Data = fscanf(fid,'%f',prod(Dims)*(isComplex+1));
  
  % combine to complex and reshape
  if isComplex
    Data = complex(Data(1:2:end),Data(2:2:end));
  end
  Data = reshape(Data,Dims);
  
  % close file
  St = fclose(fid);
  if St<0, error('Unable to close ECO file.'); end
  
  if ~isempty(Scaling)
    error('Scaling does not work for this file type.');
  end

case {'.PLT','.plt'}
  %--------------------------------------------------
  % PLT file processing
  %   MAGRES  Nijmegen EPR/ENDOR simulation program
  %--------------------------------------------------
  
  [Line,found] = findtagsMAGRES(FileName,{'DATA'});
  if found(1), nx = str2double(Line{1}); else nx=0; end
  if ~nx,
    error('Unable to determine number of x points in PLT file.');
  end
  
  fid = fopen(FileName,'r');
  if (fid<0), error(['Could not open ' FileName]); end
  
  for k=1:3, fgetl(fid); end
  
  % read data
  ny = 1;
  [Data,N] = fscanf(fid,'%f',[nx,ny]);
  if (N<nx*ny),
    warning('Could not read entire data set from PLT file.');
  end
  
  % close file
  St = fclose(fid);
  if St<0, error('Unable to close PLT file.'); end
  
  if ~isempty(Scaling)
    error('Scaling does not work for this file type.');
  end

case {'.SPK','.spk','.ref','.REF'}
  %--------------------------------------------------
  % SPK, REF file processing
  %   Varian E9 file format (ETH specific, home-built
  %   computer acquisition system written in 1991)
  %--------------------------------------------------
  fid = fopen(FileName,'r','ieee-le');
  if fid<0, error('Could not open %s.',FileName); end
  [RawData,N] = fread(fid,inf,'single');
  if fclose(fid)<0, error('Unable to close %s.',FileName); end
  
  K = [500 1e3 2e3 5e3 1e4];
  idx = find(N>K);
  if isempty(idx), error('File too small.'); end
  
  Data = RawData(N-K(idx(end))+1:end).';
  % No idea what the first part of such a file contains...
  % There is no documentation available...
  
  if ~isempty(Scaling)
    error('Scaling does not work for this file type.');
  end
  
case {'.D00','.d00'}
  %----------------------------------------------
  % d00 file processing
  %   ESE  Weizmann and ETH acquisition software
  %----------------------------------------------
  
  % Read parameter file
  % -> not implemented
  
  % open the .d00 file and error if unsuccessful
  h = fopen(FileName);
  if h<0, error(['Could not open ' FileName]); end
  
  % read in first three 16bit integers
  Dims = fread(h,3,'int16').';
  %nDims = sum(Dims>1);
  
  % read in data, complex
  Data = fread(h,[2,inf],'double');
  Data = complex(Data(1,:) ,Data(2,:));
  
  % and shape into correct array size
  Data = reshape(Data,Dims);
  
  % close data file
  St = fclose(h);
  if St<0, error('Unable to close D00 file.'); end
  %----------------------------------------------

  if ~isempty(Scaling)
    error('Scaling does not work for this file type.');
  end

case {'.d01','.D01'}
  %----------------------------------------------
  % d01 file processing
  %   SpecMan
  %----------------------------------------------  
  
  % Read parameter file
  % -> not implemented
  
  % Open the .d01 file and error if unsuccessful
  [h,ignore] = fopen(FileName,'r','ieee-le');
  if (h<0), error(['Could not open ' FileName]); end
  
  nDataSets = fread(h,1,'uint32');  % number of headers, re/im etc.
  ndim = 1;
  
  % Number format: 0-double(64bit),1-float(32bit)
  FormatID = fread(h,1,'uint32'); 
  switch FormatID
    case 1, DataFormat = 'float32';
    case 0, DataFormat = 'double';
    otherwise, error('Could not determine format in %s',FileName);
  end
 
  for iDataSet = 1:nDataSets
    ndim2(iDataSet) = fread(h,1,'int32');  % re/im ?
    dims(:,iDataSet) = fread(h,4,'int32');
    nTotal(iDataSet) = fread(h,1,'int32' );
  end
  dims(dims==0) = 1;

  Data = fread(h,sum(nTotal),DataFormat);

  %try
  switch (nDataSets)
  case 2,
    Data = complex(Data(1:nTotal),Data(nTotal+1:end));
    Data = reshape(Data,dims(:,1).');
  case 1,
    Data = reshape(Data,dims(:).');
  end
  %end
  
  % Close data file
  St = fclose(h);
  if (St<0), error(['Unable to close ' FileName]); end

  if ~isempty(Scaling)
    error('Scaling does not work for this file type.');
  end

otherwise
  
  error('Files with extension %s not supported.',Extension);
  
end


if ParseParameters
  Parameters = parseparams(Parameters);
else
  Parameters = [];
end

switch (nargout)
  case 1
    varargout = {Data};
  case 2
    varargout = {Abscissa, Data};
  case 3
    varargout = {Abscissa, Data, Parameters};
  case 4
    varargout = {Abscissa, Data, Parameters, FileName};
  case 0
    if min(size(Data))==1
      if isreal(Data)
        plot(Abscissa,Data);
      else
        plot(Abscissa,real(Data),'b',Abscissa,imag(Data),'r');
      end
      title(FileName,'Interpreter','none');
      axis tight
      if ~isreal(Data)
        legend('real','imag');
      end
    else
      pcolor(real(Data)); shading flat;
    end
end

return

%--------------------------------------------------
function out = getmatrix(FileName,Dims,DimOrder,NumberFormat,ByteOrder,isComplex)

% Open data file, error if fail.
FileID = fopen(FileName,'r',ByteOrder);
if (FileID<1), error('Unable to open data file %s',FileName); end

% Calculate expected number of elements and read in.
% Real and imaginary data are interspersed.
if isComplex, N = 2*prod(Dims); else N = prod(Dims); end
[x,effN] = fread(FileID,N,NumberFormat);
if (effN<N)
  error('Unable to read all expected data.');
end

% Close data file
CloseStatus = fclose(FileID);
if (CloseStatus<0), error('Unable to close data file %s',FileName); end

% Combine real and imaginary data to complex.
if isComplex
  x = complex(x(1:2:end),x(2:2:end));
end

% Reshape to matrix and permute dimensions if wanted.
out = ipermute(reshape(x(:),Dims(DimOrder)),DimOrder);

return
%--------------------------------------------------

%--------------------------------------------------
function [out,found] = findtagsMAGRES(FileName,TagList)

% open file
fid = fopen(FileName,'r');
if fid<0, error(['Could not open ' FileName]); end

found = zeros(1,length(TagList));
out = cell(1,length(TagList));
while ~feof(fid)
  Line = fgetl(fid);
  whitespace = find(isspace(Line)); % space or tab
  if ~isempty(whitespace),
    endTag = whitespace(1)-1;
    if endTag>0
      I = strcmp(Line(1:endTag),TagList);
      if ~isempty(I),
        out{I} = fliplr(deblank(Line(end:-1:endTag+1)));
        found(I) = 1;
      end
    end
  end
end

% close file
St = fclose(fid);
if St<0, error('Unable to close data file.'); end

return
%--------------------------------------------------


function [Parameters,err] = readPARfile(PARFileName)

Parameters = [];
err = [];

if exist(PARFileName,'file')
  allLines = textread(PARFileName,'%s','whitespace','','delimiter','\n');
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
    if Value([1 end])=='''', % remove leading and trailing quotes
      Value([1 end]) = [];
    end
  end
  
  % set field in output structure
  Parameters.(Key) = Value;
  
end

return

%---------------------------------------------------------------
function [Parameters,err] = readDSCfile(DSCFileName)

Parameters = [];
err = [];

if exist(DSCFileName,'file')
  allLines = textread(DSCFileName,'%s','whitespace','','delimiter','\n','bufsize',20000);
else
  err = sprintf('Cannot find the file %s.',DSCFileName);
  return;
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
  if ~isletter(Key(1)),
    % Stop reading when Manipulation History Layer is reached.
    if strcmpi(Key,'#MHL'); break; end
    continue;
  end
  
  Value = deblank(Value(end:-1:1)); Value = deblank(Value(end:-1:1));
  
  if ~isempty(Value)
    if Value([1 end])=='''',
       Value([1 end]) = [];
    end
  end
  
  % Set field in output structure.
  Parameters.(Key) = Value;
  
end

return

%-----------------------------------------------------------------
function P = parseparams(ParamsIn)

P = ParamsIn;

Fields = fieldnames(P);
for iField = 1:numel(Fields)
  v = P.(Fields{iField});
  if isempty(v), continue; end
  if strcmpi(v,'true')
    v_num = 1;
  elseif strcmpi(v,'false')
    v_num = 0;
  elseif isletter(v(1))
    v_num = '';
    continue;
  else
    [v_num,cnt,errormsg,nxt] = sscanf(v,'%e');
    % Converts '3345 G' to [3345] plus an error message...
    % Unclear whether conversion makes sense for the user. If not,
    % exclude such cases with
    if ~isempty(errormsg)
      v_num = '';
    end
  end
  if ~isempty(v_num)
    P.(Fields{iField}) = v_num(:).';
  end
end

return
