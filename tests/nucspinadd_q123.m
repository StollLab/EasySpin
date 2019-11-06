function err = test(opt)

% Test various combinations of Q tensor specifcations (eeQqh, eeQqh and eta,
% principal values)
% one = Q has one element only, eeQqh
% two = Q has two elements, [eeQqh eta]
% three = Q has three elements, [Q1 Q2 Q3]

Sys0.Nucs = '2H';
Sys0.A = 10;

% one + one, one + two, one + three
Sys0.Q = 2;
Sys = nucspinadd(Sys0,'2H',3,[],1);
err(1) = any(size(Sys.Q)~=[1 2]) || any(Sys.Q(2)~=1);
Sys = nucspinadd(Sys0,'2H',3,[],[1 0.1]);
err(2) = any(size(Sys.Q)~=[2 2]) || any(Sys.Q(2,:)~=[1 0.1]);
Sys = nucspinadd(Sys0,'2H',3,[],[0 -1 1]);
err(3) = any(size(Sys.Q)~=[2 3]) || any(Sys.Q(1,:)~=[-0.5 -0.5 1]);

% two + one, two + two, two + three
Sys0.Q = [2 0.2];
Sys = nucspinadd(Sys0,'2H',3,[],2);
err(4) = any(size(Sys.Q)~=[2 2]) || any(Sys.Q(2,:)~=[2 0]);
Sys = nucspinadd(Sys0,'2H',3,[],[1 0.1]);
err(5) = any(size(Sys.Q)~=[2 2]) || any(Sys.Q(2,:)~=[1 0.1]);
Sys = nucspinadd(Sys0,'2H',3,[],[0 -1 1]);
err(6) = any(size(Sys.Q)~=[2 3]) || any(Sys.Q(1,:)~=[-0.4 -0.6 1]);

% three + one, three + two
Sys0.Q = [-0.3 -0.7 1];
Sys = nucspinadd(Sys0,'2H',3,[],2);
err(7) = any(size(Sys.Q)~=[2 3]) || any(Sys.Q(2,:)~=[-0.5 -0.5 1]);
Sys = nucspinadd(Sys0,'2H',3,[],[2 0.2]);
err(8) = any(size(Sys.Q)~=[2 3]) || any(Sys.Q(2,:)~=[-0.4 -0.6 1]);

err = any(err);
