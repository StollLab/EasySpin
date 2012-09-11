% ewrls   Exponentially weighted recursive least squares adaptive filter
%
%   y = ewrls(data,p,lambda)
%   y = ewrls(data,p,lambda,nPreAvg)
%   y = ewrls(data,p,lambda,nPreAvg,delta)
%   y = ewrls(data,p,lambda,nPreAvg,delta,dir)
%   ewrls(...)
%
%   Given a series of scans in the 2D data array data, computes a filtered
%   average using the exponentially weighted recursive least squares (RLS)
%   adaptive filter method.
%   If you omit the output (y), ewrls plots the result.
%
%     data     2D data array, each column representing one scan
%     p        filter length, typically 5-50
%     lambda   memory factor, typically 0.96-0.99
%     nPreAvg  number of scans to obtain desired signal
%              if 0, all scans are averaged. default 0
%     delta    regularization parameter (default 100)
%     dir      filtering direction (default 'fb')
%              'f' forward, 'b' backward, 'fb' average of both
%     y        filtered averaged scans, a single column
%
%   Example:
%     y = ewrls(data,32,0.96)

function y = ewrls(data,p,lambda,nPreAvg,delta,direc)

nScans = size(data,2);
nPoints = size(data,1);

if nargin<6, direc = 'fb'; end
if nargin<5, delta = 100; end
if nargin<4, nPreAvg = 0; end
if nargin<3
  error('At least three inputs (data, p, lambda) are needed.');
end

if (nScans<=1), error('Multiple scans (columns) are needed. First input must be a 2D array with more than one column.'); end
if (nPoints<=1), error('Multiple points (rows) are needed. First input must be a 2D array with individual scans along columns.'); end
if (p>=nPoints), error('Filter length (%d) is longer than data (%d), must be shorter.',p,nPoints); end

if (lambda<0)||(lambda>1)
  error('lambda (third input) must be between 0 and 1. Typical values are around 0.96-0.98.');
end
if (p<=0)
  error('Filter length p (second input) must be positive. Typical values are 10-50.');
end

if strcmp(direc,'f')
  forwardsDirection = 1;
  backwardsDirection = 0;
elseif strcmp(direc,'b')
  forwardsDirection = 0;
  backwardsDirection = 1;  
elseif strcmp(direc,'fb') || strcmp(direc,'bf')
  forwardsDirection = 1;
  backwardsDirection = 1;  
else
  error('Unknown string ''%s''in direction parameter',direc);
end

% desired signal: average of scans
if (nPreAvg==0)
  d = mean(data,2);
else
  nPreAvg = min(nPreAvg,nScans);
  d = mean(data(:,1:nPreAvg),2);
end

% EWRLS filtering
%-------------------------------------------------

if forwardsDirection
  y = zeros(nPoints,nScans);
  P = delta*eye(p); % initialize P
  w = zeros(p,1); % initialize filter coefficients
  for iScan = nPreAvg+1:nScans
    for k = p:nPoints
      x = data(k-p+1:k,iScan);
      phi = P*x; % auxiliary vector
      g = phi/(lambda+x.'*phi);
      y(k,iScan) = x.'*w;
      e = d(k) - y(k,iScan);
      w = w + g*e;
      P = (P-g*phi.')/lambda;
      P = (P+P.')/2; % symmetrize, to prevent divergence due to numerical instability
    end
  end
  y(:,1:nPreAvg) = [];
  y_forward = mean(y,2);
end

if backwardsDirection
  data = data(end:-1:1,:);
  d = d(end:-1:1);
  y = zeros(nPoints,nScans);
  P = delta*eye(p); % initialize P
  w = zeros(p,1); % initialize filter coefficients
  for iScan = nPreAvg+1:nScans
    for k = p:nPoints
      x = data(k-p+1:k,iScan);
      phi = P*x; % auxiliary vector
      g = phi/(lambda+x.'*phi);
      y(k,iScan) = x.'*w;
      e = d(k) - y(k,iScan);
      w = w + g*e;
      P = (P-g*phi.')/lambda;
      P = (P+P.')/2; % symmetrize, to prevent divergence due to numerical instability
    end
  end
  y(:,1:nPreAvg) = [];
  y_backward = mean(y,2);
  y_backward = y_backward(end:-1:1);
  data = data(end:-1:1,:);
  d = d(end:-1,1);
end

if forwardsDirection
  y = y_forward;
  if backwardsDirection
    y = (y+y_backward)/2;
  end
else
  y = y_backward;
end

% Plot results
%-------------------------------------------------------------
if (nargout==0)
  clf
  yavg_f = y;
  yavg_c = mean(data,2);
  t = 1:nPoints;
  plot(t,yavg_c,'r',t,yavg_f,'b');
  title('System output') ;
  legend('conventional avg','filtered avg'); legend boxoff
  axis tight
end
