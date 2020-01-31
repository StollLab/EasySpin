function err = test()

A = rand(2);
nI = 3;
B1 = runprivate('eyekron',nI,A,false);
B2 = runprivate('eyekron',nI,A,true);

err = ~areequal(B1,B2,1e-15,'abs');
