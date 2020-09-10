function [Triangulation,Areas] = triangles(nOctants,nKnots,Vecs)

% Vecs           3xn double array   orientation vectors
% nKnots         double             number of knots
% nOctants       double             number of octants
% Triangulation  3xn uint32 array   triangle index array
% Areas          1xn double array   triangle areas

if size(Vecs,1)~=3
  error('Vector array must have 3 rows!');
end

% Determine triangulation
Triangulation = sphtri(nOctants,nKnots,'f').';

if size(Triangulation,1)~=3
  error('Triangulation array must have 3 rows!');
end

nVertices = size(Vecs,2);
if nVertices~=max(Triangulation(:))
  error('Number of vertices and maximum triangulation index do not match!');
end

% Compute areas of spherical triangles
%--------------------------------------------
% Vertex vectors
x1 = Vecs(:,Triangulation(1,:));
x2 = Vecs(:,Triangulation(2,:));
x3 = Vecs(:,Triangulation(3,:));

% Edge arc lengths
a1 = acos(sum(x2.*x3));
a2 = acos(sum(x3.*x1));
a3 = acos(sum(x1.*x2));

% Formula of d'Huilier
s = (a1+a2+a3)/2; % triangle perimeter half
Areas = 4*atan(sqrt(tan(s/2).*tan((s-a1)/2).*tan((s-a2)/2).*tan((s-a3)/2)));

% Normalise to sum 4*pi
Areas = Areas/sum(Areas) * 4*pi;

if ~isreal(Areas)
  error('Complex triangle areas encountered! (nOctants %d, nKnots %d)!',nOctants,n);
end

return
