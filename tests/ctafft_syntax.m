function ok = test()

TD = rand(1,1024);
TD2 = rand(128,128);

fd = ctafft(TD,4);
fd = ctafft(TD,1:5:30);
fd = ctafft(TD,1:5:30,2048);
fd = ctafft(TD2,1:4:20);
fd = ctafft(TD2,1:4:20,512);

ok = true;

