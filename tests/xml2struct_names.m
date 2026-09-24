function ok = test()

% Element and attribute names that are not valid MATLAB field names

longName = repmat('a',1,namelengthmax+10);
xml = ['<r><Gr' char(246) 'sse ' longName '="2">1</Gr' char(246) 'sse><x-y.z>3</x-y.z></r>'];
fileName = [tempname '.xml'];
writexml(xml,fileName);
cleanup = onCleanup(@()delete(fileName));

s = runprivate('xml2struct',fileName);

ok = strcmp(s.r.Gr_sse.Text,'1') && ...
     strcmp(s.r.Gr_sse.Attributes.(longName(1:namelengthmax)),'2') && ...
     strcmp(s.r.x_dash_y_dot_z.Text,'3');

end

%-------------------------------------------------------------------------------
function writexml(xml,fileName)
fid = fopen(fileName,'w');
fwrite(fid,unicode2native(xml,'UTF-8'));
fclose(fid);
end
