function ok = test()

% Element with both text and child elements keeps both

fileName = [tempname '.xml'];
writexml('<a>hello<b>x</b></a>',fileName);
cleanup = onCleanup(@()delete(fileName));

s = runprivate('xml2struct',fileName);

ok = isfield(s.a,'Text') && strcmp(s.a.Text,'hello') && ...
     isfield(s.a,'b') && strcmp(s.a.b.Text,'x');

end

%-------------------------------------------------------------------------------
function writexml(xml,fileName)
fid = fopen(fileName,'w');
fwrite(fid,unicode2native(xml,'UTF-8'));
fclose(fid);
end
