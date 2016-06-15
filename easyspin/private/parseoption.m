function [out,err] = parseoption(Stru,Field,Values)

% Parses Fields specified as field Field in structure Stru, possibly
% setting default values. NOT case sensitive!

err = [];
out = [];

if isfield(Stru,Field)
  thisValue = Stru.(Field);
  found = strcmpi(thisValue,Values);
  nFound = sum(found);
  if nFound==0
    err = ['Unknown value ''',thisValue,''' for field ',Field,'.'];
  elseif nFound>1
    err = ['Ambiguous value ''',thisValue,''' for field ',Field,'.'];
  else
    out = find(found);
  end
else
  err = ['Missing field ''', Field, '.'];
end

if nargout<2
  error(err);
end

return
