function [Abscissa,Data,Parameters] = eprload_specman(FileName)

%----------------------------------------------
% Files from SpecMan, *.d01 extension
%----------------------------------------------
% Parameter file import is not implemented.

Abscissa = [];
Data = [];
Parameters = [];

% Open the .d01 file and error if unsuccessful
h = fopen(FileName,'r','ieee-le');
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

% Close data file
status = fclose(h);
if (status<0), error(['Unable to close ' FileName]); end

switch (nDataSets)
  case 2,
    Data = complex(Data(1:nTotal),Data(nTotal+1:end));
    Data = reshape(Data,dims(:,1).');
  case 1,
    Data = reshape(Data,dims(:).');
end

