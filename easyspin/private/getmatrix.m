%-------------------------------------------------------------------------------
function out = getmatrix(FileName,Dims,NumberFormat,ByteOrder,isComplex)

% Open data file, error if fail.
FileID = fopen(FileName,'r',ByteOrder);
if (FileID<1), error('Unable to open data file %s',FileName); end

% Calculate expected number of elements and read in.
% Real and imaginary data are interspersed.
nDataValuesPerPoint = numel(isComplex);
nRealsPerPoint = sum(isComplex+1);
N = nRealsPerPoint*prod(Dims);

[datalist,effN] = fread(FileID,N,NumberFormat);
if (effN~=N)
  error('Unable to read all expected data.');
end

% Close data file
CloseStatus = fclose(FileID);
if (CloseStatus<0), error('Unable to close data file %s',FileName); end

% Reshape data and combine real and imaginary data to complex.
datalist = reshape(datalist,nRealsPerPoint,prod(Dims));
for k = 1:nDataValuesPerPoint
  if isComplex(k)
    data{k} = complex(datalist(k,:),datalist(k+1,:)).';
    datalist(k+1,:) = [];
  else
    data{k} = datalist(k,:);
  end
end

% Reshape to matrix and permute dimensions if wanted.
for k = 1:nDataValuesPerPoint
  out{k} = reshape(data{k},Dims);
end

if numel(out)==1, out = out{1}; end

return
%-------------------------------------------------------------------------------
