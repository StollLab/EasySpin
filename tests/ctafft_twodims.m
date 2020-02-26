function ok = test()

% Check output size for 2D case

TD = rand(512,20);
FD = ctafft(TD,1:5);

ok = all(size(FD)==size(TD));
