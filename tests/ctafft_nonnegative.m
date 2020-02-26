function ok = test()

TD = rand(1,1024);
FD = ctafft(TD,1:5,2048);

ok = isreal(FD) && all(FD(:)>0);
