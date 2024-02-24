function ok = test()
% Check that g-, D-, A- and Q-frames are loaded correctly for ORCA v4 from
% both the main output file and the binary property file
% -------------------------------------------------------------------------------------------------

% Main output file
% --------------------------------------------------------------------------------------
% Load file that contains a planar quinoline triplet with coordinates rotated
% with respect to the symmetry axes in the ORCA frame.
Sys = orca2easyspin('orca\quinoline_triplet_tilted_v4.out');

% The "molecular frame" in this script refers to the ORCA frame, i.e. the
% XYZ frame used to define the atom locations.

% Extract needed atomic positions (in molecular frame)
xyz = Sys.xyz.';
N = xyz(:,1);
C4 = xyz(:,4);
C9 = xyz(:,9);

% Calculate oop axis (in molecular frame; within sign)
oopaxis_M = cross(C4-N,C9-N);
oopaxis_M = oopaxis_M/norm(oopaxis_M);

% Get transformation matrices from molecular frame to tensor frames

% gx in oop direction
% -----------------------------------------------------------------------
x_T = [1;0;0];  % x axis of tensor, represented in tensor frame
R_M2g = erot(Sys.gFrame);
R_g2M = R_M2g.';
gx_M = R_g2M*x_T;  % z axis of g tensor, represented in molecular frame

% Assert that g tensor x axis is along out-of-plane axis
ok(1) = areequal(abs(gx_M.'*oopaxis_M),1,1e-2,'rel');

% Dz in oop direction
% -----------------------------------------------------------------------
z_T = [0;0;1];  % z axis of tensor, represented in tensor frame
R_M2D = erot(Sys.DFrame);
R_D2M = R_M2D.';
Dz_M = R_D2M*z_T;  % z axis of D tensor, represented in molecular frame

% Assert that D tensor z axis is along out-of-plane axis
ok(2) = areequal(abs(Dz_M.'*oopaxis_M),1,1e-3,'rel');

% Az in oop direction
% -----------------------------------------------------------------------
z_T = [0;0;1];  % z axis of tensor, represented in tensor frame
R_M2A = erot(Sys.AFrame);
R_A2M = R_M2A.';
Az_M = R_A2M*z_T;  % z axis of A tensor, represented in molecular frame

% Assert that A tensor z axis is along out-of-plane axis
ok(3) = areequal(abs(Az_M.'*oopaxis_M),1,1e-3,'rel');

% Qy in oop direction
% -----------------------------------------------------------------------
y_T = [0;1;0];  % y axis of tensor, represented in tensor frame
R_M2Q = erot(Sys.QFrame);
R_Q2M = R_M2Q.';
Qy_M = R_Q2M*y_T;  % z axis of Q tensor, represented in molecular frame

% Assert that Q tensor z axis is along out-of-plane axis
ok(4) = areequal(abs(Qy_M.'*oopaxis_M),1,1e-3,'rel');


% Binary property file
% --------------------------------------------------------------------------------------
% Load file that contains a planar quinoline triplet with coordinates rotated
% with respect to the symmetry axes in the ORCA frame.
Sys = orca2easyspin('orca\quinoline_triplet_tilted_v4.prop');

% The "molecular frame" in this script refers to the ORCA frame, i.e. the
% XYZ frame used to define the atom locations.

% Extract needed atomic positions (in molecular frame)
xyz = Sys.xyz.';
N = xyz(:,1);
C4 = xyz(:,4);
C9 = xyz(:,9);

% Calculate oop axis (in molecular frame; within sign)
oopaxis_M = cross(C4-N,C9-N);
oopaxis_M = oopaxis_M/norm(oopaxis_M);

% Get transformation matrices from molecular frame to tensor frames

% gx in oop direction
% -----------------------------------------------------------------------
x_T = [1;0;0];  % x axis of tensor, represented in tensor frame
R_M2g = erot(Sys.gFrame);
R_g2M = R_M2g.';
gx_M = R_g2M*x_T;  % z axis of g tensor, represented in molecular frame

% Assert that g tensor x axis is along out-of-plane axis
ok(5) = areequal(abs(gx_M.'*oopaxis_M),1,1e-2,'rel');

% Dz in oop direction
% -----------------------------------------------------------------------
z_T = [0;0;1];  % z axis of tensor, represented in tensor frame
R_M2D = erot(Sys.DFrame);
R_D2M = R_M2D.';
Dz_M = R_D2M*z_T;  % z axis of D tensor, represented in molecular frame

% Assert that D tensor z axis is along out-of-plane axis
ok(6) = areequal(abs(Dz_M.'*oopaxis_M),1,1e-3,'rel');

% Az in oop direction
% -----------------------------------------------------------------------
z_T = [0;0;1];  % z axis of tensor, represented in tensor frame
R_M2A = erot(Sys.AFrame);
R_A2M = R_M2A.';
Az_M = R_A2M*z_T;  % z axis of A tensor, represented in molecular frame

% Assert that A tensor z axis is along out-of-plane axis
ok(7) = areequal(abs(Az_M.'*oopaxis_M),1,1e-3,'rel');

% Qy in oop direction
% -----------------------------------------------------------------------
y_T = [0;1;0];  % y axis of tensor, represented in tensor frame
R_M2Q = erot(Sys.QFrame);
R_Q2M = R_M2Q.';
Qy_M = R_Q2M*y_T;  % z axis of Q tensor, represented in molecular frame

% Assert that Q tensor z axis is along out-of-plane axis
ok(8) = areequal(abs(Qy_M.'*oopaxis_M),1,1e-3,'rel');

end