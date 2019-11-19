% Converts a string shorthand for an direction to the corresponding unit vector.

function n = letter2vec(s)

if ~ischar(s)
  error('First argument must be a character array.');
end

switch s
  case 'x', n = [1;0;0];
  case 'y', n = [0;1;0];
  case 'z', n = [0;0;1];
  case 'xy', n = [1;1;0];
  case 'xz', n = [1;0;1];
  case 'yz', n = [0;1;1];
  case 'xyz', n = [1;1;1];
  case '-x', n = -[1;0;0];
  case '-y', n = -[0;1;0];
  case '-z', n = -[0;0;1];
  case '-xy', n = -[1;1;0];
  case '-xz', n = -[1;0;1];
  case '-yz', n = -[0;1;1];
  case '-xyz', n = -[1;1;1];
  otherwise
    error('Unknown value ''%s'' for rotation axis (2nd input argument).',s);
end

n = n/norm(n);
