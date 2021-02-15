% eprsave  Save data in Bruker BES3T format
%
%   eprsave(FileName,x,data)
%   eprsave(FileName,x,data,TitleString)
%   eprsave(FileName,x,data,TitleString,mwFreq)
%
%   Saves the dataset in x and data in the Bruker BES3T
%   format in a .DTA and a .DSC file with the file
%   name given in FileName. x is the x axis data, and data
%   is the intensity data (real or complex). TitleString is
%   the name of the dataset that will be displayed in the
%   Bruker software. mwFreq is the microwave frequency, in GHz.
%
%   Two-dimensional data can be saved by giving a matrix in data,
%   and by supplying both axes in x as a cell array.
%
%   Examples:
%     eprsave(myFilename,B,spc);       % save 1D data spc, x axis = B 
%     eprsave(myFilename,{t1,t2},V);   % save 2D data matrix V,
%                                      %   x axis = t1, y axis = t2

function eprsave(filename,x,data,TitleString,mwFreq)

if nargin==0, help(mfilename); return; end

if nargin<4, TitleString = ''; end
if nargin<5, mwFreq = NaN; end

if ismatrix(data)
  TwoDimData = min(size(data))>1;
else
  error('Cannot save data with more than 2 dimensions.');
end

if TwoDimData
  if ~iscell(x) || numel(x)~=2
    error('For two-dimensional data, x must be a cell array containing the two axes.');
  end
  y = x{2};
  x = x{1};
else
  if iscell(x)
    x = x{1};
  end
end

complexData = ~isreal(data);

BES3TVersion = 1.2;

if BES3TVersion<1.2
  error('Cannot save in BES3T format versions older than 1.2.');
end

if length(filename)>3
  if strcmpi(filename(end-3:end),'.DSC') || ...
      strcmpi(filename(end-3:end),'.DTA')
    filename = filename(1:end-4);
  end
end

NumberFormat = 'float64';
ByteOrder = 'ieee-be'; % big-endian is default for XEPR (Linux)

% Save data in DTA file
%-----------------------------------------------------------
DTAfilename = [filename '.DTA'];
fDTA = fopen(DTAfilename,'w',ByteOrder);
if fDTA<1, error('Unable to open data file %s',DTAfilename); end
datalist = data(:);
if complexData
  datalist = [real(datalist) imag(datalist)].';
end
fwrite(fDTA,datalist(:),NumberFormat,0,ByteOrder);
CloseStatus = fclose(fDTA);
if CloseStatus<0, error('Unable to close data file %s',DTAfilename); end

% Save parameters in DSC file
%-----------------------------------------------------------
DSCfilename = [filename '.DSC'];
fDSC = fopen(DSCfilename,'w',ByteOrder);
if fDSC<1, error('Unable to open description file %s',DSCfilename); end
writedsckeyval = @(key,val)fprintf(fDSC,[key '\t' val '\n']);
writedsc = @(val)fprintf(fDSC,[val '\n']);

fprintf(fDSC,'* Exported from Matlab using EasySpin, %s\n',datestr(now));
if BES3TVersion==1.2
  VersionString = '1.2';
elseif BES3TVersion==1.3
  VersionString = '1.3';
elseif BES3TVersion==2.0
  VersionString = '2.0';
end
writedsckeyval('#DESC',[VersionString ' * DESCRIPTOR INFORMATION ***********************']);

writedsckeyval('DSRC','MAN');

if strcmp(ByteOrder,'ieee-be')
  writedsckeyval('BSEQ','BIG');
else
  writedsckeyval('BSEQ','LIT');
end

if complexData
  writedsckeyval('IKKF','CPLX');
else
  writedsckeyval('IKKF','REAL');
end

if strcmp(NumberFormat,'float64')
  NumberFormatCode = 'D';
else
  error('Unsupported number format ''%s''',NumberFormat);
end

writedsckeyval('IRFMT',NumberFormatCode);
if complexData
  writedsckeyval('IIFMT',NumberFormatCode);
end

% Axis information
%------------------------------------------------------------
if BES3TVersion<2.0
  GaugeFileExt = {'.XGF','.YGF','.ZGF'};
else
  GaugeFileExt = {'.GF1','.GF2','.GF3'};
end

% Determine if X axis is linear
deviationFromLinear = @(d) max(diff(d(:)))/min(diff(d(:)));
linThreshold = 1e-5;
isLinearAxis = @(d) abs(deviationFromLinear(d)-1)<linThreshold;

isLinearX = isLinearAxis(x);
if ~isLinearX, XType = 'IGD'; else, XType = 'IDX'; end

% Determine if Y axis is linear
if TwoDimData
  isLinearY = isLinearAxis(y);
  if ~isLinearY, YType = 'IGD'; else, YType = 'IDX'; end
else
  YType = 'NODATA';
end

ZType = 'NODATA';

% Write companion index-gauge file for X axis if necessary
if ~isLinearX
  fGF = fopen([filename GaugeFileExt{1}],'w',ByteOrder);
  fwrite(fGF,x(:),'float64',0,ByteOrder);
  fclose(fGF);
  if BES3TVersion<2.0
    writedsckeyval('XFMT','D');
  else
    writedsckeyval('AX1FMT','D');
  end
end

% Write companion index-gauge file for Y axis if necessary
if TwoDimData && ~isLinearY
  fGF = fopen([filename GaugeFileExt{2}],'w',ByteOrder);
  fwrite(fGF,y(:),'float64',0,ByteOrder);
  fclose(fGF);
  if BES3TVersion<2.0
    writedsckeyval('YFMT','D');
  else
    writedsckeyval('AX2FMT','D');
  end
end

% Write X/Y/Z axis types
if BES3TVersion<2.0
  writedsckeyval('XTYP',XType);
  writedsckeyval('YTYP',YType);
  writedsckeyval('ZTYP',ZType);
else
  writedsckeyval('AX1TYP',XType);
  %writedsckeyval('AX2TYP',YType);
  %writedsckeyval('AX3TYP',ZType);
end

% Write X/Y/Z axis characteristics
if BES3TVersion<2.0
  writedsckeyval('XPTS',sprintf('%d',numel(x)));
  writedsckeyval('XMIN',sprintf('%g',x(1)));
  writedsckeyval('XWID',sprintf('%g',x(end)-x(1)));
  if TwoDimData
    writedsckeyval('YPTS',sprintf('%d',numel(y)));
    writedsckeyval('YMIN',sprintf('%g',y(1)));
    writedsckeyval('YWID',sprintf('%g',y(end)-y(1)));
  end
else
  writedsckeyval('AX1PTS',sprintf('%d',numel(x)));
  writedsckeyval('AX1MIN',sprintf('%g',x(1)));
  writedsckeyval('AX1WID',sprintf('%g',x(end)-x(1)));
  if TwoDimData
    writedsckeyval('AX1PTS',sprintf('%d',numel(y)));
    writedsckeyval('AX1MIN',sprintf('%g',y(1)));
    writedsckeyval('AX1WID',sprintf('%g',y(end)-y(1)));
  end
end

% Write title
if ~isempty(TitleString)
  writedsckeyval('TITL',['''' TitleString '''']);
end

% Standard parameter layer (SPL)
%------------------------------------------------------------------
writedsc('*');
writedsc('************************************************************');
writedsc('*');
writedsckeyval('#SPL',[VersionString ' * STANDARD PARAMETER LAYER']);
writedsckeyval('OPER','');
writedsckeyval('DATE',datestr(now,'dd/mm/yy'));
writedsckeyval('TIME',datestr(now,'HH:MM:SS'));
writedsckeyval('CMNT','');
writedsckeyval('SAMP','');
writedsckeyval('SFOR','');
if ~isnan(mwFreq)
  mwFreqString = sprintf('%0.9g',mwFreq*1e9); % BES3T requires Hz
  writedsckeyval('MWFQ',mwFreqString);
end

CloseStatus = fclose(fDSC);
if CloseStatus<0
  error('Unable to close data file %s',DSCfilename);
end

