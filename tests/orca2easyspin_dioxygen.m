function ok = test()

% Test main ouput file against property output file.

Sys1 = orca2easyspin('orca/dioxygen_singlepoint.oof');
Sys2 = orca2easyspin('orca/dioxygen_singlepoint_property.txt');

ok = all(cellfun(@(a,b)isequal(a,b),fieldnames(Sys1),fieldnames(Sys2)));
