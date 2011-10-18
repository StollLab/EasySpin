% zfsframes   D tensor in various frame conventions
%
%    zfsframes(D1,D2,D3)
%    zfsframes([D1 D2 D3])
%
%    Determines and displays the zero-field splitting
%    parameters D and E for all possible choices of the
%    eigenframe of the D tensor, with the given principal
%    values D1, D2, and D3. Conventions are listed.
%
%    Example:
%      zfsframes(-1,-1,2)   % axial tensor
%      zfsframes(-1,0,1)    % rhombic tensor
%      zfsframes(-1,0,2)    % general tensor

function zfsframes(varargin)

switch nargin
  case 0, help(mfilename); return;
  case 1, D = varargin{1};
  case 3, D = [varargin{1:3}];
  otherwise
    error('Please provide three principal values of D tensor.');
end

fprintf('\n==================================================================\n');
fprintf('D tensor in various axis system conventions\n');
fprintf('==================================================================\n');
fprintf('Principal values (PVs) of D tensor in input order:\n  D1 = %+g     D2 = %+g     D3 = %+g\n',D(1),D(2),D(3));

if sum(D)~=0
  meanD = mean(D);
  D = D - meanD;
  fprintf('PVs don''t sum up to 0, subtracting mean (%f).\n',meanD);
  fprintf('PVs of traceless D tensor in input order:\n  D1 = %+g     D2 = %+g     D3 = %+g\n',D(1),D(2),D(3));
else
  fprintf('This D tensor is traceless (sum of PVs is zero).\n');
end

if all(D==0)
  fprintf('All D values are zero.\n');
  return
end

if (D(1)==D(2)) || (D(1)==D(3)) || D(2)==D(3)
  disp('This D tensor is exactly axial (two PVs are equal).');
elseif any(D==0)
  disp('This D tensor is exactly orthorhombic (one PV is zero).');
else
  disp('This D tensor is neither exactly axial nor exactly orthorhombic.');
end

disp('Principal axis system conventions:');
%disp('  Blumberg1 1967:   |E/D| has smallest possible value, and E>=0');
disp('  Poole 1974:    |Dz|>=|Dx|>=|Dy|    (used for organic triplets)');
disp('  Blumberg 1967: |Dz|>=|Dy|>=|Dx|    (used for transition metals)');
disp('  (for references see documentation)');
disp('Computation of D and E:');
disp('  D = 3/2*Dz      E = (Dx-Dy)/2');
disp('  Dx = -D/3+E     Dy = -D/3-E        Dz = 2/3*D');
%disp('  Gaffney 1993:     |Dx|>=|Dy|>=|Dz|    (very exotic)');

ax = perms([3 2 1]);
for k=1:6
  Dx = D(ax(k,1));
  Dy = D(ax(k,2));
  Dz = D(ax(k,3));
  D_(k) = 3/2*Dz;
  E_(k) = (Dx-Dy)/2;
  name{k} = '';
  %if abs(Dz)>=max(abs(Dx),abs(Dy)) & E_(k)>=0-eps, name{k} = [name{k} 'Blumberg ']; end
  if abs(Dz)>=abs(Dx) && abs(Dx)>=abs(Dy), name{k} = 'Poole '; end
  if abs(Dz)>=abs(Dy) && abs(Dy)>=abs(Dx), name{k} = [name{k} 'Blumberg ']; end
  %if abs(Dx)>=abs(Dy) & abs(Dy)>=abs(Dz), name{k} = [name{k} 'Gaffney ']; end
end

%ED = E_./D_;
%for k=1:6
%  if abs(ED(k))==min(abs(ED)) & E_(k)>=0, name{k} = [name{k} 'Blumberg1']; end
%end

warning('off','MATLAB:divideByZero');
ED_ = E_./D_;

disp('D and E in all possible principal axis systems:');
for k=1:6
  ax_ = ax(k,:);
  fprintf('  %d%d%d->xyz:  D = %-+5.3f    E = %-+5.3f    E/D = %-+5.2f    %s\n',...
    ax_(1),ax_(2),ax_(3),D_(k),E_(k),ED_(k),name{k});
end
disp('The Blumberg convention is the recommended one.')
disp(' ');

return
