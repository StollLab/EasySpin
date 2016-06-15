% adddefaults  Add default fields to user structure
%
%   SupplementedStruct = adddefaults(UserStruct,DefaultStruct)
%
%   Supplements missing fields in the user-given structure
%   UserStruct using fields and values from the default
%   stucture DefaultStruct.
%
%   Example:
%    UserStruct: a=1, b=1, d=1
%    DefaultStruct: a=2, c=2, e=4
%    SupplementedStruct: a=1, b=1, c=2, d=1, e=4

function SupplementedStruct = adddefaults(UserStruct,DefaultStruct)

SupplementedStruct = DefaultStruct;

if isempty(UserStruct), return; end

UserStruct = spellcheckfields(UserStruct,fieldnames(DefaultStruct));

% Transfer user fields to default structure, possibly overwriting defaults.
SupplementedStruct = DefaultStruct;
if ~isempty(UserStruct)
  Fields = fieldnames(UserStruct);
  for iField = 1:numel(Fields)
    Name = Fields{iField};
    SupplementedStruct.(Name) = UserStruct.(Name);
  end
end

return

% spellcheckfields   Checks the upper/lowercase spelling of fields
%
%   Out = spellcheckfields(Structure,field1,field2,field3,...)
%
%   In Structure, field names that match one of the
%   ones given as parameters (case insensitive), are
%   changed to the one matching (case sensitive).
%   This allows to adjust upper/lower case. The user
%   can write the fields in any upper/lower combination.

function Out = spellcheckfields(Structure,ExactFieldNames)

Out = [];

if isempty(Structure); return; end

Fields = fieldnames(Structure);

for iF = 1:numel(Fields)
  Name = Fields{iF};
  idxCorrect = find(strcmp(Name,ExactFieldNames)); % field present, correct capitalization
  idxAny = find(strcmpi(Name,ExactFieldNames)); % field present, any capitalization

  % If there is a field with incorrect capitalization, but none with the
  % correct one, issue an error.
  if isempty(idxCorrect) && ~isempty(idxAny)
    error('\n   Supplied field ''%s'' is similar to ''%s''.\n   Please correct spelling or remove field.',Name,ExactFieldNames{idxAny(1)});
  end

  % If there is a field with correct capitalization, take that
  % and discard all fields with other capializations
  if ~isempty(idxCorrect)
    Name = ExactFieldNames{idxCorrect};
  end
  Value = Structure.(Name);
  Out.(Name) = Value;
end

return
