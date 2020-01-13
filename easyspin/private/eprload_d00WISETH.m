%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_d00WISETH(FileName)
%-------------------------------------------------------------------------------
% d00 file processing
%   ESE Weizmann and ETH acquisition software
%-------------------------------------------------------------------------------

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

Abscissa = [];
Parameters = [];
return
%-------------------------------------------------------------------------------

