% eprsave  Save data in Bruker BES3T format
%
%   eprsave(FileName,x,y)
%   eprsave(FileName,x,y,TitleString)
%
%   Saves the dataset in x and y in the Bruker BES3T
%   format in a .DTA and a .DSC file with the file
%   name given in FileName. x is the x axis data, and y
%   is the intensity data. TitleString is the name of the
%   dataset that will be displayed in the Bruker software.

function eprsave(filename,x,y,titlestring)

if (nargin==0); help(mfilename); return; end

if nargin<4
  titlestring = '';
end

%{
if (nargin==0)
  x = 1:100;
  y = gaussian(x,mean(x),(max(x)-min(x))/5);
  x = x.^2;
  %x = rand(10);
  y = y + 1i*linspace(0,max(y),numel(x));
  filename = 'test6nonlin';
  titlestring = 'EasySpin test example';
  plot(x,real(y),'b',x,imag(y),'r');
  axis tight
end
%}

if ndims(y)>2
  error('Cannot save data with more than 2 dimensions.')
end
if ndims(y)==2
  if min(size(y))>1
    error('Cannot save two-dimensional data.')
  end
end

complexData = ~isreal(y);

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
fDTA = fopen([filename '.DTA'],'w',ByteOrder);
if complexData
  yy = [real(y(:)) imag(y(:))].';
else
  yy = y;
end
fwrite(fDTA,yy(:),NumberFormat,0,ByteOrder);
fclose(fDTA);

% Save parameters in DSC file
%-----------------------------------------------------------
fDSC = fopen([filename '.DSC'],'w',ByteOrder);
writedsc = @(key,val)fprintf(fDSC,[key '\t' val '\n']);

fprintf(fDSC,'* Exported from Matlab using EasySpin, %s\n',datestr(now));

if BES3TVersion==1.2
  VersionString = '1.2';
elseif BES3TVersion==1.3
  VersionString = '1.3';
elseif BES3TVersion==2.0
  VersionString = '2.0';
end
writedsc('#DESC',VersionString);

writedsc('DSRC','MAN');

if strcmp(ByteOrder,'ieee-be')
  writedsc('BSEQ','BIG');
else
  writedsc('BSEQ','LIT');
end

if complexData
  writedsc('IKKF','CPLX');
else
  writedsc('IKKF','REAL');
end

if strcmp(NumberFormat,'float64')
  NumberFormatCode = 'D';
else
  error('Unsupported number format ''%s''',NumberFormat);
end

writedsc('IRFMT',NumberFormatCode);
if complexData
  writedsc('IIFMT',NumberFormatCode);
end

% Axis information
%------------------------------------------------------------
if BES3TVersion<2.0
  GaugeFileExt = {'.XGF','.YGF','.ZGF'};
else
  GaugeFileExt = {'.GF1','.GF2','.GF3'};
end

% X axis: Determine if axis is linear
DeviationFromLinear = abs(x(:) - linspace(x(1),x(end),numel(x)).');
isLinearX = max(DeviationFromLinear)/abs(x(end)-x(1))<0.001;
if ~isLinearX, XType = 'IGD'; else XType = 'IDX'; end
YType = 'NODATA';
ZType = 'NODATA';

% Write companion file
if ~isLinearX
  fGF = fopen([filename GaugeFileExt{1}],'w',ByteOrder);
  fwrite(fGF,x(:),'float64',0,ByteOrder);
  fclose(fGF);
  if BES3TVersion<2.0
    writedsc('XFMT','D');
  else
    writedsc('AX1FMT','D');
  end
end

% Write axis types
if (BES3TVersion<2.0)
  writedsc('XTYP',XType);
  writedsc('YTYP',YType);
  writedsc('ZTYP',ZType);
else
  writedsc('AX1TYP',XType);
  %writedsc('AX2TYP',YType);
  %writedsc('AX3TYP',ZType);
end

% Write linear axis characteristics
if BES3TVersion<2.0
  writedsc('XPTS',sprintf('%d',numel(x)));
  writedsc('XMIN',sprintf('%g',x(1)));
  writedsc('XWID',sprintf('%g',x(end)-x(1)));
else
  writedsc('AX1PTS',sprintf('%d',numel(x)));
  writedsc('AX1MIN',sprintf('%g',x(1)));
  writedsc('AX1WID',sprintf('%g',x(end)-x(1)));
end

if ~isempty(titlestring)
  writedsc('TITL',['''' titlestring '''']);
end

fclose(fDSC);
