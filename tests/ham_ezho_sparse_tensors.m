function ok = test()

% Assert that all elements of the tensors returned by ham_ezho are sparse
% if requested, and full otherwise.

rng(5);

Sys.S = 1;
Sys.Ham312 = rand(1,5);

Gs = ham_ezho(Sys,[],[],'sparse',3);
G3s = Gs{1};
Gf = ham_ezho(Sys,[],[],'',3);
G3f = Gf{1};

ok(1) = issparse(G3s{1,1,1}) && issparse(G3s{1,1,2}) && issparse(G3s{1,2,3});
ok(2) = ~issparse(G3f{1,1,1}) && ~issparse(G3f{1,1,2}) && ~issparse(G3f{1,2,3});
