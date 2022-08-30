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

function easyspin(str)

if nargin==0
  str = 'info';
end

switch str
  case '?'
    dips(' easyspin info       Display information about EasySpin');
    disp(' easyspin doc        Display EasySpin documentation');
    disp(' easyspin compile    Compile MEX files for EasySpin');
  case 'info'
    easyspininfo;
  case 'doc'
    easyspindoc;
  case 'compile'
    easyspincompile;
  otherwise
    error('Unknown option ''%s''.',str);
end

end
