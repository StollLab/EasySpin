%-------------------------------------------------------------------------------
% Parse fields of a structure and convert to numbers if possible
%-------------------------------------------------------------------------------

function Pout = parsefieldparams(ParamsIn)

Pout = ParamsIn;

Fields = fieldnames(Pout);
for iField = 1:numel(Fields)
  v = Pout.(Fields{iField});
  if isempty(v), continue; end
  if strcmpi(v,'true')
    v_num = true;
  elseif strcmpi(v,'false')
    v_num = false;
  elseif isletter(v(1))
    v_num = '';
    continue
  else
    [v_num,cnt,errormsg,nxt] = sscanf(v,'%e');
    % Converts '3345 G' to [3345] plus an error message...
    % Unclear whether conversion makes sense for the user. If not,
    % exclude such cases with
    if ~isempty(errormsg)
      v_num = '';
    end
  end
  if ~isempty(v_num)
    Pout.(Fields{iField}) = v_num(:)'; % don't use .' due to bug up to R2014a
  end
end

return
