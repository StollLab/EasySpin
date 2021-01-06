% esspline1d  1D spline interpolation [EasySpin/private]
%
%    pp = esspline1d(y);
%    pp = esspline1d(y,EndCon);
%    yy = esspline1d(y,EndCon,nn);
%    yy = esspline1d(y,EndCon,xx);
%
%    The first two forms return a piecewise polynomial of the
%    spline interpolant of the columns of y over x = 1:size(y,2).
%    The last two forms return interpolated values. If
%    nn (a scalar) is given, xx is set to linspace(1,...
%    size(y,2),nn). If given explicitly, xx must be a vector.
%
%      EndCon = 0   not-a-knot (default)
%      EndCon = 1   endslopes zero
%                   First derivatives at first and last
%                   point are set to zero.
%      EndCon = 2   periodic
%                   Matches first and second derivative
%                   at first and last point.
%
%    Call spparms('autommd',0) before calling ESSLPINE1D (speed-up)!

% Matlab M-file dependencies
%   linspace, mkpp, unmkpp, spdiags

function out = esspline1d(y,EndCon,xx)
%--------------------------------------------------------
% This function is basically a contraction and
% specialization of Matlab's spline and ppval functions
% with a tuned interface. pp creation with esspline is
% ca 60% faster than spline in cases relevant to EasySpin
% functions like pepper and salt. yy creation is 30% faster.

% Supplement default for boundary condition setting.
if nargin<2
  EndCon = 0;
end

% Data is assumed to be along rows.
y = y.';
[n,nDims] = size(y);

if (n<2), error('Interpolation not possible: n < 2.'); end

% Compute the differences of adjacent data.
dy = diff(y);

if (n==2) % Two points
  
  if (EndCon~=1) % natural, periodic: the interpolant is a straight line
    %out = mkpp(1:n,[dy.' y(1,:).'],nDims);
    x = [1 2];
    c = [zeros(nDims,2) dy.' y(1,:).' ];
  else         % endslopes 0: the interpolant is the cubic Hermite polynomial
    dy2 = diff([zeros(1,nDims);dy; zeros(1,nDims)]);
    %out = mkpp((1:n).',...
    %  [(diff(dy2)).' ([2 -1]*dy2).' zeros(1,nDims).' y(1,:).'],nDims);
    x = 1:n;
    c = [(diff(dy2)).' ([2 -1]*dy2).' zeros(1,nDims).' y(1,:).'];
  end
  
elseif (n==3) && (EndCon==0) % Three points
  % The interpolant is a parabola.
 
  y(2:3,:) = dy;
  y(3,:) = diff(dy)/2;
  y(2,:) = y(2,:) - y(3,:);
  x = [1,3];
  c = [zeros(nDims,1) y([3 2 1],:).'];
  %out = mkpp([1,3],y([3 2 1],:).',nDims);
  
else
  
  % Set up linear equation system A*slp = b.
  % Most of A and b is independent of the end conditions.
  %A = spdiags([[on;0;0] [1;4*on;1] [0;0;on]],-1:1,n,n);
  A = sparse([2:n-1 2:n-1 1:n],[1:n-2 3:n 1:n],[ones(1,2*(n-2)+1) 4*ones(1,n-2) 1],n,n);
  b = zeros(n,nDims);
  b(2:n-1,:) = 3 * (dy(1:n-2,:) + dy(2:n-1,:));
  
  % Boundary-condition specific elements of A and b
  switch EndCon
  case 0 % natural
    A(1,2) = 2;
    A(n,n-1) = 2;
    b(1,:) = .5 * (5*dy(1,:)+dy(2,:));
    b(n,:) = .5 * (dy(n-2,:)+5*dy(n-1,:));
  case 1 % endslopes 0p
    A(1,2) = 0;
    A(n,n-1) = 0;
  case 2 % periodic
    A(1,n) = -1;
    A(n,[1 2 n-1 n]) = [2 1 1 2];
    b(n,:) = 3 * (dy(1,:) + dy(n-1,:));
  end
  
  % Sparse linear equation solution for the slopes.
  slp = A\b;
  
  % Compute 4th and 3rd order polynomial coefficients.
  c4 = slp(1:n-1,:) + slp(2:n,:) - 2*dy(1:n-1,:);
  c3 = dy(1:n-1,:) - slp(1:n-1,:) - c4;
  x = 1:n;
  c = reshape([c4.' c3.' slp(1:n-1,:).' ...
      y(1:n-1,:).'], (n-1)*nDims,4);
  
end


if (nargin<3) % return interpolated values
    % Collect piecewise polynomial stucture
  out = mkpp(x,c,nDims);
else

  % parse third parameter
  if (numel(xx)<=1)
    nn = xx;
    xs = linspace(1,n,nn);
  else
    nn = length(xx);
    xs = xx(:).';
  end
  
  npp = length(x)-1;
  ord = 4;
  
  %[x,c,npp,ord,nDims] = unmkpp(out);
  
  % do the interpolation
  
  % for each data point, compute its pp interval
  index = fix(xs);
  index(index>npp) = npp; % !!
  % go to local coordinates
  xs = xs - index;
  
  if (nDims>1) % replicate xs and index in case pp is vector-valued
    xs = reshape(xs(ones(nDims,1),:),nDims*nn,1);
    index = nDims*index;
    temp = (-nDims:-1).';
    index = reshape(1 + index(ones(nDims,1),:) + temp(:,ones(1,nn)), nDims*nn, 1 );
  else
    xs = xs(:);
  end
  
  % Evaluate polynomial
  c = c(index,:);
  yy = c(:,1);
  for i = 2:ord
    yy = xs.*yy + c(:,i);
  end 
  out = reshape(yy.',nDims,nn);
end

return

%-------------------------------------------------------------
% Test, comparing against spline and csape.

spparms('autommd',0);
n = 20; nd = 20;
x = 1:n;
y = rand(nd,n);
xx = 1:.02:n;
subplot(3,1,1);
plot(xx,spline(x,y,xx),'b.',xx,esspline1d(y,0,xx),'r');
subplot(3,1,2);
plot(xx,spline(x,[zeros(nd,1) y zeros(nd,1)],xx),'b.',xx,esspline1d(y,1,xx),'r');
subplot(3,1,3);
plot(xx,fnval(csape(x,y(:,[1:end-1,1]),'p'),xx),'b.',xx,esspline1d(y(:,[1:end-1,1]),2,xx),'r');

% Speed test against spline
%---------------------------------------
spparms('autommd',0);
n = 20; nd = 5;
x = (1:n).';
kk = 100;
xx = 1:.1:n;
y = rand(nd,n);
yz = [zeros(nd,1) y zeros(nd,1)];
t1 = cputime; for k=1:kk, p1=esspline1d(y,1); end; t1 = cputime-t1;
t2 = cputime; for k=1:kk, p2=spline(x,yz); end; t2 = cputime-t2; t1/t2
t1 = cputime; for k=1:kk, c1=esspline1d(y,1,xx); end; t1 = cputime-t1;
t2 = cputime; for k=1:kk, c2=spline(x,yz,xx); end; t2 = cputime-t2; t1/t2
