function ok = test()

% Check whether ctafft operates along columns

TD = rand(512,2);
TD(:,3) = 0;
FD = ctafft(TD,1:5);
ok = all(FD(:,3)==0);
