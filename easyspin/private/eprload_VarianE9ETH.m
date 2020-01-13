%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_VarianE9ETH(FileName)
%-------------------------------------------------------------------------------
% SPK, REF file processing
%   Varian E9 file format (ETH specific, home-built
%   computer acquisition system written in 1991)
%-------------------------------------------------------------------------------
fid = fopen(FileName,'r','ieee-le');
if fid<0, error('Could not open %s.',FileName); end
[RawData,N] = fread(fid,inf,'float32');
if fclose(fid)<0, error('Unable to close %s.',FileName); end

K = [500 1e3 2e3 5e3 1e4];
idx = find(N>K);
if isempty(idx), error('File too small.'); end

Data = RawData(N-K(idx(end))+1:end).';
% No idea what the first part of such a file contains...
% There is no documentation available...

Abscissa = [];
Parameters = [];
return
%-------------------------------------------------------------------------------
