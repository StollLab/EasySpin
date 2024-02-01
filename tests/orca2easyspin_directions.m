function [err,data] = test(opt,olddata)

% Load file that contains a planar methyl radical with C3 axis tilted away
% from the z axis of the ORCA frame.
Sys = orca2easyspin('orca\methyl_tilted_all_v4.prop');

% The "molecular frame" in this script refers to the ORCA frame, i.e. the
% XYZ frame used to define the atom locations.

% Extract needed atomic positions (in molecular frame)
xyz = Sys.xyz.';
C = xyz(:,1);
H1 = xyz(:,2);
H2 = xyz(:,3);

% Calculate C3 symmetry axis (in molecular frame; within sign)
C3axis_M = cross(H1-C,H2-C);
C3axis_M = C3axis_M/norm(C3axis_M);

% Get transformation matrices from molecular frame to A tensor frames
zA_A = [0;0;1];  % z axis of A tensor, represented in A frame
for k = 1:4
  R_M2A = erot(Sys.AFrame(k,:));
  R_A2M = R_M2A.';
  zA_M{k} = R_A2M*zA_A;  % z axis of A tensor, represented in molecular frame
end

% Assert that carbon A tensor z axis is along C3 symmetry axis
ok(1) = areequal(abs(zA_M{1}.'*C3axis_M),1,1e-3);

% Assert that hydrogen A tensor z axes are perpendicular to C3 symmetry axis
for k = 2:4
  ok(k) = areequal(zA_M{k}.'*C3axis_M,0,1e-3);
end

err = any(~ok);
data = [];