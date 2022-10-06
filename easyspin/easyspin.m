% easyspin   Toolbox for Electron Paramagnetic Resonance (EPR)
%
%  easyspin
%  easyspin ?
%  easyspin doc
%  easyspin info
%  easyspin compile
%
%   If a parameter is given, various tasks are performed:
%     ?        show options
%     doc      display EasySpin documentation
%     info     display information about current EasySpin installation
%     compile  compile MEX files for EasySpin
%
%  If no parameter is given, 'info' is used by default.

function varargout = easyspin(str)

if nargin==0
  str = 'info';
end

switch str

  case '?'
    disp(' easyspin info       Display information about EasySpin');
    disp(' easyspin doc        Display EasySpin documentation');
    disp(' easyspin compile    Compile MEX files for EasySpin');
    varargout = {};

  case 'info'
    if nargout==0
      easyspin_info;
    else
      out = easyspin_info;
      varargout = {out};
    end

  case 'doc'
    if nargout==0
      easyspin_doc;
    else
      out = easyspin_doc;
      varargout = {out};
    end

  case 'compile'
    easyspin_compile;

  otherwise
    error('Unknown option ''%s''.',str);

end

end
