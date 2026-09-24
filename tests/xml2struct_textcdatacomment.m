function ok = test()

% Text, CDATA, and comment content are stored in separate fields

fileName = [tempname '.xml'];
writexml('<a>t1<![CDATA[c]]><!--m-->t2</a>',fileName);
cleanup = onCleanup(@()delete(fileName));

s = runprivate('xml2struct',fileName);

ok = strcmp(s.a.Text,'t1t2') && strcmp(s.a.CDATA,'c') && strcmp(s.a.Comment,'m');

end

%-------------------------------------------------------------------------------
function writexml(xml,fileName)
fid = fopen(fileName,'w');
fwrite(fid,unicode2native(xml,'UTF-8'));
fclose(fid);
end
