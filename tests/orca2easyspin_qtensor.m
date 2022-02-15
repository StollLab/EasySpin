function ok = test()

% Test whether Q tensor is read in correctly from ORCA property file

Sys = orca2easyspin('orca/hydroxyl_098.oof');
Q0 = [-0.0644363 -0.0896791 0.1541154]; % MHz
ok = all(abs(sort(Sys.Q)-sort(Q0))<1e-5);
