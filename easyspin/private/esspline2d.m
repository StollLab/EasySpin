% esspline2d  2D spline interpolation [EasySpin/private]
%
%    yy = esspline2d(y,rr,cc);
%    yy = esspline2d(y,rr,cc,EndCon);
%
%    y, a real 2D array, is interpolated using cubic
%    tensor product splines. y is assumed to defined over
%    r = 1:size(y,1) and c = 1:size(y,2). rr and cc are
%    the row and column vectors for the interpolation points.
%    EndCon is a 2-vector defining the boundary conditions as
%    in esspline. EndCon(1) is for interpolation along columns,
%    EndCon(2) along rows.
%
%      EndCon(.) = 0   not-a-knot (default)
%      EndCon(.) = 1   endslopes zero
%                      First derivatives at first and last
%                      point are set to zero.
%      EndCon(.) = 2   periodic
%                      Matches first and second derivative
%                      at first and last point.
%
%    Call spparms('autommd',0) before calling ESSLPINE2D (speed-up)!

% M-file dependencies
%   esspline1d
%   unmkpp, num2cell, repmat, sub2ind

%============================================================
function F = esspline2d(data,ri,ci,EndCon)
%============================================================
% Specialized from interp2, spline2 and splncore

if nargin<4, EndCon = [0 0]; end

points = [ri(:).'; ci(:).'];
n = length(ri);

% Construct spline representation
%------------------------------------------------------------

% (1) Interpolate along rows of data
% size(coefs) is [d*L,k]. d is dimension size(data(1)), k is
% polynom order (4), L is number of polynom pieces size(data(2))-1
[ignored,coefs,rL,k,d] = unmkpp(esspline1d(data,EndCon(2)));
% Prepend zero coefficients if order is not cubic
if k<4, coefs = [zeros(d,4-k) coefs]; end
% Reorder coefficients to [L*k,d]
values = reshape(coefs,d,rL*4).';
sizeg(2) = rL*4;

% Interpolate along rows of reordered coefficient array, ie
% interpolate coefficients along dimensions.
[ignored,coefs,cL,k,d] = unmkpp(esspline1d(values,EndCon(1)));
if k<4, coefs = [zeros(d,4-k) coefs]; end
values = reshape(coefs,d,cL*4).';
sizeg(1) = cL*4;

% Locate interpolation points
%------------------------------------------------------------
% for each data point, compute its pp interval
ipoint = fix(points);
ipoint(1,ipoint(1,:)>cL) = cL;
ipoint(2,ipoint(2,:)>rL) = rL;
rpoint = points - ipoint;

% Arrange coefficients
%--------------------------------------------------------------
ells = sizeg/4;
idx = ells(1)*[0:3].';
mm = repmat(ells(2),4,1);
idx = [idx mm*0; idx mm; idx mm*2; idx mm*3];
idx = num2cell(1+idx,1);
offset = reshape(sub2ind(sizeg,idx{:}),16,1);
base = ipoint(1,:) - 1 + sizeg(1)*(ipoint(2,:)-1);

coefs = reshape(values(base(ones(16,1),:)+offset(:,ones(1,n))),[4,4,n]);

% Compute interpolated values
%-------------------------------------------------------------------
for i = 2:-1:1
  s = reshape(rpoint(i,:),[1,1,n]);
  coefs = reshape(coefs,[4^(i-1),4,n]);
  ss = s(ones(4^(i-1),1),1,:);
  coefs = ((coefs(:,1,:).*ss + coefs(:,2,:)).*ss + coefs(:,3,:)).*ss + coefs(:,4,:);
end

F = reshape(coefs,size(ri));

return

%=============================================================================
% Test against csape from the Spline Toolbox

clear, spparms('autommd',0);
n = 10; m = 10;
nn = 5000;
r = 1:n; c = 1:m; z = rand(n,m);
rr = 1 + rand(1,nn)*(n-1);
cc = 1 + rand(1,nn)*(m-1);
% not-a-knot along both directions
t = cputime; zz1 = fnval(csape({r,c},z,'n'),[rr;cc]); t1 = cputime - t;
t = cputime; zz2 = esspline2d(z,rr,cc); t2 = cputime - t;
fprintf('(n,n): error %f, relative time %f\n',max(abs(zz1(:)-zz2(:))), t2/t1);
% zero-endslope along rows, natural along cols
zcs = zeros(n,m+1); zcs(:,2:end) = z; zcs(end,end+1) = 0;
t = cputime; zz1 = fnval(csape({r,c},zcs,{'n','c'}),[rr;cc]); t1 = cputime - t;
t = cputime; zz2 = esspline2d(z,rr,cc,[0 1]); t2 = cputime - t;
fprintf('(n,z): error %f, relative time %f\n',max(abs(zz1(:)-zz2(:))), t2/t1);
% periodic along rows, zero-endslopes along cols
z1 = z; z1(:,end) = z(:,1);
zcs = zeros(n+1,m); zcs(2:end,:) = z1; zcs(end+1,end) = 0;
t = cputime; zz1 = fnval(csape({r,c},zcs,{'c','p'}),[rr;cc]); t1 = cputime - t;
t = cputime; zz2 = esspline2d(z1,rr,cc,[1 2]); t2 = cputime - t;
fprintf('(z,p): error %f, relative time %f\n',max(abs(zz1(:)-zz2(:))), t2/t1);
