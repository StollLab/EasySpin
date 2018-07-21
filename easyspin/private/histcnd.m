% Copyright (c) 2011, Kota Yamaguchi
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [ n, bin ] = histcnd( x, edges )
%HISTCND Histogram count for n dimensional data.
% N = HISTCND(X,EDGES), for row vectors X, counts the number of values in
% X that fall between the grid defined by the cell array of EDGES, each
% of whose element is a vector that contain monotonically non-decreasing
% values. N is an N-D array each of whose dimension corresponds to
% LENGTH(EDGES{j}) and each element contains a count of data that falls
% into the edge.
%
% EDGES must have the same length to the number of columns of X.
% Alternatively, EDGES can be a numeric vector which gives a uniform
% grid for all dimensions of X.
% 
% N(k1,k2,...) will count the vector X(i,:) if for each dimension
% j = 1,2,..., EDGES{j}(kj) <= X(i,j) < EDGES{j}(kj+1). The last bin will
% count any values of X that match EDGES(end).  Values outside the values
% in EDGES are not counted. Use -inf and inf in EDGES to include all
% non-NaN values.
% 
% [N,BIN] = HISTCND(X,EDGES) also returns subscript indices BIN.
% BIN is zero for out of range values.
% 
% Example:
%     >> X = randn(100,2);          % 100-by-2 row vectors
%     >> edges = {-2:.4:2,-2:.5:2}; % ranges for each dimension
%     >> histcnd(X,edges)
% 
%     ans =
% 
%          0     0     0     1     1     1     0     0     0
%          0     1     1     2     1     1     0     0     0
%          1     0     3     4     0     3     0     0     0
%          0     1     2     1     3     0     1     0     0
%          0     3     1     4     2     1     3     1     0
%          1     1     2     3     3     4     1     0     0
%          0     1     1     2     1     1     4     0     0
%          0     1     2     2     2     1     0     0     0
%          1     2     0     3     2     0     0     1     0
%          0     1     1     0     0     0     0     0     0
%          0     0     0     0     0     0     0     0     0
% 
% Class support for inputs X:
%    float: double, single
% 
% See also histc.

% Validation
if ~isnumeric(x)
    error('Input array must be numeric.');
end
if ~iscell(edges)
    if isnumeric(edges) && isvector(edges)
        if isempty(x)
            n = [];
            bin = [];
            return;
        end
        tmp = cell(1,size(x,2));
        tmp(:) = {edges};
        edges = tmp;
    else
        error('Edges must be cell array.');
    end
end
if isempty(x)
    x = reshape(x,[0 length(edges)]);
end
if length(edges)~=size(x,2)
    error('Invalid cell array.');
end

% Compute dimension of the output histogram
dims = cellfun(@(c) length(c), edges(:)');
if length(dims)==1, dims = [1 dims]; end

% Quantize input along each dimension
bin = cell(1,size(x,2));
outlier = false(size(x,1),1);
for i = 1:size(x,2)
    [tmp,bin{i}] = histc(x(:,i),edges{i});
    outlier = outlier | (bin{i}==0);
end

% Remove outlier before computing a linear index
ind = bin;
for i = 1:size(x,2)
    ind{i} = ind{i}(~outlier);
end

% Compute count of the index
n = histc(sub2ind(dims,ind{:}),1:prod(dims));

% Format the output
n = reshape(n,dims);
bin = cat(2,bin{:});

end
