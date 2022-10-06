% easyspin_doc  Provide access to the EasySpin documentation

function varargout = easyspin_doc()

% Determine entry point for documentation
esPath = fileparts(which(mfilename));
esRoot = esPath(1:end-length('\easyspin\private'));
docEntry = [esRoot filesep 'documentation' filesep 'index.html'];
if ~exist(docEntry,'file')
  docEntry = [esRoot filesep 'docsrc' filesep 'index.html'];
end

if nargout==0
  fprintf('<a href="%s">EasySpin documentation</a>\n',docEntry);
  web(docEntry,'-new');
  varargout = {};
else
  varargout = {docEntry};
end

end
