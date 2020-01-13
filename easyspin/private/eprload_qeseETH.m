%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_qeseETH(FileName)
%-------------------------------------------------------------------------------
% ECO file processing
%   qese     old ETH acquisition software
%   tryscore Weizmann HYSCORE simulation program
%-------------------------------------------------------------------------------

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

Abscissa = [];
Parameters = [];

return
%-------------------------------------------------------------------------------
