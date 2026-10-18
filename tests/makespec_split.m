function ok = test()

% Splitting between neighboring points keeps total amplitude and centroid

rng(1);
N = 1000;
Pos = rand(1,N);
Amp = rand(1,N);
[x,s] = makespec([0 1],51,Pos,Amp,'split');

tol = 1e-12;
ok(1) = areequal(sum(s),sum(Amp),tol,'rel');
ok(2) = areequal(sum(x.*s)/sum(s),sum(Pos.*Amp)/sum(Amp),tol,'rel');

% Single peak: shares go to the two neighboring points
[~,s] = makespec([0 1],11,0.425,1,'split');
ok(3) = areequal(s([5 6]),[0.75 0.25],tol,'abs') && nnz(s)==2;

% Peak exactly at the upper limit goes entirely to the last point
s = makespec([0 1],11,1,1,'split');
ok(4) = s(end)==1 && nnz(s)==1;
