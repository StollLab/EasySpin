function ok = test()
% Check for match of Sys structures extracted from ORCA files for different
% versions and for different types of output files
% -------------------------------------------------------------------------------------------------

% Load file that contains a planar quinoline triplet with coordinates rotated
% with respect to the symmetry axes in the ORCA frame.
% ORCA v5, main output file
Sys{1} = orca2easyspin('orca\quinoline_triplet_tilted.out');
% ORCA v4, main output file
Sys{2} = orca2easyspin('orca\quinoline_triplet_tilted_v4.out');
% ORCA v4, binary property file
Sys{3} = orca2easyspin('orca\quinoline_triplet_tilted_v4.prop');

% Get transformation matrices from molecular frame to tensor frames
for i = 1:numel(Sys)
  R_M2g{i} = erot(Sys{i}.gFrame);
  R_M2D{i} = erot(Sys{i}.DFrame);
  R_M2A{i} = erot(Sys{i}.AFrame);
  R_M2Q{i} = erot(Sys{i}.QFrame);
end

% Check for agreement between Sys structures extracted from different files
id = [1 2; 2 3];
fieldnames = {'g','D','A','Q'};
numfields = numel(fieldnames);
for i = 1:size(id,1)

  % Principal values
  for j = 1:numfields
    oki(j) = areequal(Sys{id(i,1)}.(fieldnames{j}),Sys{id(i,2)}.(fieldnames{j}),1e-2,'rel');
  end

  % Transformation matrices
  oki(numfields+1) = areequal(abs(R_M2D{id(i,1)}),abs(R_M2D{id(i,2)}),1e-2,'rel');
  oki(numfields+2) = areequal(abs(R_M2A{id(i,1)}),abs(R_M2A{id(i,2)}),1e-2,'rel');
  oki(numfields+3) = areequal(abs(R_M2Q{id(i,1)}),abs(R_M2Q{id(i,2)}),1e-2,'rel');
  % oki(numfields+4) = areequal(abs(R_M2g{id(i,1)}),abs(R_M2g{id(i,2)}),1e-2,'rel'); %
  % differences due mainly to different ORCA versions

  ok(i) = all(oki==1);

end

end