function ok = test()

% File name handling: omitted extension (also for very short names),
% string input, and folder instead of file

folder = tempname;
mkdir(folder);
writexml('<a>1</a>',fullfile(folder,'x.xml'));
oldFolder = cd(folder);
cleanup = onCleanup(@()cleanupFolder(oldFolder,folder));

s1 = runprivate('xml2struct','x');
s2 = runprivate('xml2struct',"x.xml");

try
  runprivate('xml2struct',folder);
  folderError = false;
catch
  folderError = true;
end

ok = strcmp(s1.a.Text,'1') && isequal(s1,s2) && folderError;

end

%-------------------------------------------------------------------------------
function writexml(xml,fileName)
fid = fopen(fileName,'w');
fwrite(fid,unicode2native(xml,'UTF-8'));
fclose(fid);
end

%-------------------------------------------------------------------------------
function cleanupFolder(oldFolder,folder)
cd(oldFolder);
rmdir(folder,'s');
end
