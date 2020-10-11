% eprload  Load experimental EPR data
%
%   y = eprload(FileName)
%   [x,y] = eprload(FileName)
%   [x,y,Pars] = eprload(FileName)
%   [x,y,Pars,FileN] = eprload(FileName)
%   ... = eprload(FileName,Scaling)
%   ... = eprload
%
%   Read spectral data from a file specified in the string 'FileName' into the
%   arrays x (abscissa) and y (ordinate). The structure Pars contains entries
%   from the parameter file, if present.
%
%   All strings in the parameter structure containing numbers are converted to
%   numbers for easier use.
%
%   If FileName is a directory, a file browser is displayed. If FileName is
%   omitted, the current directory is used as default. eprload returns the
%   name of the loaded file (including its path) as fourth parameter FileN.
%
%   For DSC/DTA data, x contains the vector or the vectors specifying the
%   abscissa or abscissae of the spectral data array, i.e. magnetic field range
%   for cw EPR, RF range for ENDOR and time delays for pulse EPR. Units are
%   those specified in the parameter file. See the fields XPTS, XMIN, XWID
%   etc. in the Pars structure.
%
%   Supported formats are identified via the extension in 'FileName':
%
%     Bruker BES3T:        .DTA, .DSC
%     Bruker ESP, WinEPR:  .spc, .par
%     SpecMan:             .d01, .exp
%     Magnettech:          .spe (binary), .xml (xml)
%     Active Spectrum:     .ESR
%     Adani:               .dat, .json
%     JEOL:                (no extension)
%
%     MAGRES:              .PLT
%     qese, tryscore:      .eco
%     Varian:              .spk, .ref
%     ESE:                 .d00, .exp
%
%     For reading general ASCII formats, use textscan(...)
%
%   'Scaling' tells eprload to scale the data (works only for Bruker files):
%
%      'n':   divide by number of scans
%      'P':   divide by square root of microwave power in mW
%      'G':   divide by receiver gain
%      'T':   multiply by temperature in kelvin
%      'c':   divide by conversion/sampling time in milliseconds

function varargout = eprload(FileName,Scaling)

if nargout<0 || nargout>4
  error('Please provide 1, 2, 3 or 4 output arguments!');
end

if nargin<1, FileName = pwd; end
if nargin<2, Scaling = ''; end

% Convert string to char array
if isstring(FileName)
  FileName = char(FileName);
end
LocationType = exist(FileName,'file');

if LocationType==7 % a directory
  CurrDir = pwd;
  cd(FileName);
  [uiFile,uiPath] = uigetfile({...
    '*.DTA;*.dta;*.spc','Bruker (*.dta,*.spc)';...
    '*.d01','SpecMan (*.d01)';...
    '*','All files, incl. JEOL (*.*)';...
    '*.spe;*.xml','Magnettech (*.spe,*.xml)';...
    '*.esr','Active Spectrum (*.esr)';...
    '*.spk;*.ref','Varian (*.spk,*.ref)';...
    '*.eco','qese/tryscore (*.eco)';...
    '*.d00','ETH/WIS (*.d00)';...
    '*.plt','Magres (*.plt)'},...
    'Load EPR data file...');
  cd(CurrDir);
  if uiFile==0
    varargout = cell(1,nargout);
    return;
  end
  FileName = [uiPath uiFile];
end

% Decompose file name, supply default extension .DTA
[p,Name,FileExtension] = fileparts(FileName);
FullBaseName = fullfile(p,Name);

if isempty(FileExtension) || length(FileExtension) < 2
  if exist([FullBaseName '.dta'],'file'), FileExtension = '.dta'; end
  if exist([FullBaseName '.DTA'],'file'), FileExtension = '.DTA'; end
  if exist([FullBaseName '.spc'],'file'), FileExtension = '.spc'; end
end

FileName = [FullBaseName FileExtension];
LocationType = exist(FileName,'file');
if any(LocationType==[0 1 5 8]) % not a file/directory
  error('The file or directory %s does not exist!',FileName);
end

% Determine file format from file extension
switch upper(strtrim(FileExtension))
  case {'.DTA','.DSC'}, FileFormat = 'BrukerBES3T';
  case {'.PAR','.SPC'}, FileFormat = 'BrukerESP';
  case '.D01', FileFormat = 'SpecMan';
  case '.SPE', FileFormat = 'MagnettechBinary';
  case '.XML', FileFormat = 'MagnettechXML';
  case '.ESR', FileFormat = 'ActiveSpectrum';
  case '.DAT', FileFormat = 'AdaniDAT';
  case '.JSON', FileFormat = 'AdaniJSON';
  case '.ECO', FileFormat = 'qese/tryscore';
  case '.PLT', FileFormat = 'MAGRES';
  case {'.SPK','.REF'}, FileFormat = 'VarianETH';
  case '.D00', FileFormat = 'WeizmannETH';
  otherwise
    % Test for JEOL file
    h = fopen(FileName,'r');
    if (h<0), error(['Could not open ' FileName]); end
    processType = fread(h,16,'*char');  % first 16 characters
    fclose(h);
    ok = regexp(processType.','^spin|^cAcqu|^endor|^pAcqu|^cidep|^sod|^iso|^ani','once');
    if ~isempty(ok)
      FileFormat = 'JEOL';
    else
      error('Unsupported file extension %s',FileExtension);
    end
end

% Scaling works only for par/spc files (ECS106, ESP, etc)
if ~isempty(Scaling) && (strcmp(FileFormat,'BrukerBES3T') || strcmp(FileFormat,'BrukerESP'))
  S_ = Scaling;
  S_(S_=='n' | S_=='P' | S_=='G' | S_=='T' | S_=='c') = [];
  if ~isempty(S_)
    error('Scaling can only contain ''n'', ''P'', ''G'', ''T'', and ''c''.');
  end
elseif ~isempty(Scaling)
  error('Scaling does not work for this file type.');
end

% Read data
switch FileFormat
  case 'BrukerBES3T'
    [Data,Abscissa,Parameters] = eprload_BrukerBES3T(FullBaseName,FileExtension,Scaling);
  case 'BrukerESP'
    [Data,Abscissa,Parameters] = eprload_BrukerESP(FullBaseName,FileExtension,Scaling);
  case 'SpecMan'
    [Data,Abscissa,Parameters] = eprload_specman(FileName);
  case 'MagnettechBinary'
    [Data,Abscissa,Parameters] = eprload_MagnettechBinary(FileName);
  case 'MagnettechXML'
    [Data,Abscissa,Parameters] = eprload_MagnettechXML(FileName);
  case 'ActiveSpectrum'
    [Data,Abscissa,Parameters] = eprload_ActiveSpectrum(FileName);
  case 'AdaniDAT'
    [Data,Abscissa,Parameters] = eprload_AdaniDAT(FileName);
  case 'AdaniJSON'    
    [Data,Abscissa,Parameters] = eprload_AdaniJSON(FileName);
  case 'JEOL'
    [Data,Abscissa,Parameters] = eprload_jeol(FileName);
  case 'qese/tryscore'
    [Data,Abscissa,Parameters] = eprload_qeseETH(FileName);
  case 'MAGRES'
    [Data,Abscissa,Parameters] = eprload_MAGRES(FileName);
  case 'VarianETH'
    [Data,Abscissa,Parameters] = eprload_VarianE9ETH(FileName);
  case 'WeizmannETH'
    [Data,Abscissa,Parameters] = eprload_d00WISETH(FileName);
  otherwise
    error('File format ''%s'' not implemented.',FileFormat);
end

switch nargout
  case 1
    varargout = {Data};
  case 2
    varargout = {Abscissa, Data};
  case 3
    varargout = {Abscissa, Data, Parameters};
  case 4
    varargout = {Abscissa, Data, Parameters, FileName};
end

doPlot = nargout==0;
if ~doPlot, return; end

if doPlot
  if isempty(Data), return; end
  if ~iscell(Data), Data = {Data}; end
  nDataSets = numel(Data);
  for k = 1:nDataSets
    subplot(nDataSets,1,k);
    if min(size(Data{k}))==1
      if isreal(Data{k})
        plot(Abscissa,Data{k});
      else
        plot(Abscissa,real(Data{k}),'b',Abscissa,imag(Data{k}),'r');
      end
      if (nDataSets>1)
        title([FileName sprintf(', dataset %d',k)],'Interpreter','none');
      else
        title(FileName,'Interpreter','none');
      end
      axis tight
      if ~isreal(Data{k})
        legend('real','imag');
        legend boxoff
      end
    else
      pcolor(real(Data{k})); shading flat;
    end
  end
end

return
