% xml2struct       Convert XML file into a MATLAB structure
%
%     s = xml2struct(fileName)
%
%  Input:
%     fileName    file name of XML file
%
%  Output:
%     s           A hierarchical MATLAB structure representing the XML tree
%
%  A file containing:
%    <XMLname attrib1="Some value">
%      <Element>Some text</Element>
%      <Data attrib2="2">Some more text</Data>
%      <Data attrib3="8.3" attrib4="1">Even more text</Data>
%    </XMLname>
%
%  will produce:
%     s.XMLname.Attributes.attrib1 = "Some value";
%     s.XMLname.Element.Text = "Some text";
%     s.XMLname.Data{1}.Attributes.attrib2 = "2";
%     s.XMLname.Data{1}.Text = "Some more text";
%     s.XMLname.Data{2}.Attributes.attrib3 = "8.3";
%     s.XMLname.Data{2}.Attributes.attrib4 = "1";
%     s.XMLname.Data{2}.Text = "Even more text";
%
% The following characters are substituted:
%   '-'     '_dash_'
%   ':'     '_colon_'
%   '.'     '_dot_'
%
% XML elements with names #text, #comment, and #cdata-section are skipped.

% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012
% Modified by Stefan Stoll, University of Washington, Feb 2018

function s = xml2struct(fileName)

if nargin < 1
  help(mfilename);
  return
end

if isa(fileName, 'org.apache.xerces.dom.DeferredDocumentImpl') || ...
   isa(fileName, 'org.apache.xerces.dom.DeferredElementImpl')
  % Input is a java XML object
  xDoc = fileName;
else
  % Check for existence
  if exist(fileName,'file') == 0
    % Perhaps the xml extension was omitted from the file name.
    % Add the extension and try again.
    if ~strcmp(fileName(end-3:end),'.xml')
      fileName = [fileName '.xml'];
    end
    if exist(fileName,'file') == 0
      error(['The file ' fileName ' could not be found.']);
    end
  end
  % Read the xml file
  xDoc = xmlread(fileName);
end

%parse xDoc into a MATLAB structure
s = parseChildNodes(xDoc);

end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct;
textflag = 'Text';
if hasChildNodes(theNode)
  childNodes = getChildNodes(theNode);
  numChildNodes = getLength(childNodes);
  
  for count = 1:numChildNodes
    theChild = item(childNodes,count-1);
    [text,name,attr,childs,textflag] = getNodeData(theChild);
    
    if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
      %XML allows the same elements to be defined multiple times,
      %put each in a different cell
      if (isfield(children,name))
        if (~iscell(children.(name)))
          %put existsing element into cell format
          children.(name) = {children.(name)};
        end
        index = length(children.(name))+1;
        %add new element
        children.(name){index} = childs;
        if(~isempty(fieldnames(text)))
          children.(name){index} = text;
        end
        if(~isempty(attr))
          children.(name){index}.('Attributes') = attr;
        end
      else
        %add previously unknown (new) element to the structure
        children.(name) = childs;
        if(~isempty(text) && ~isempty(fieldnames(text)))
          children.(name) = text;
        end
        if(~isempty(attr))
          children.(name).('Attributes') = attr;
        end
      end
    else
      ptextflag = 'Text';
      if (strcmp(name, '#cdata_dash_section'))
        ptextflag = 'CDATA';
      elseif (strcmp(name, '#comment'))
        ptextflag = 'Comment';
      end
      
      %this is the text in an element (i.e., the parentNode)
      if ~isempty(regexprep(text.(textflag),'[\s]*',''))
        if ~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag))
          ptext.(ptextflag) = text.(textflag);
        else
          %what to do when element data is as follows:
          %<element>Text <!--Comment--> More text</element>
          
          %put the text in different cells:
          % if (~iscell(ptext)) ptext = {ptext}; end
          % ptext{length(ptext)+1} = text;
          
          %just append the text
          ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
        end
      end
    end
    
  end
end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr)))
  attr = [];
end

%parse child nodes
[childs,text,textflag] = parseChildNodes(theNode);

if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
  % Get the data of any childless nodes
  % faster than if any(strcmp(methods(theNode), 'getData'))
  % no need to try-catch (?)
  % faster than text = char(getData(theNode));
  text.(textflag) = toCharArray(getTextContent(theNode))';
end

end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = struct;
if hasAttributes(theNode)
  theAttributes = getAttributes(theNode);
  numAttributes = getLength(theAttributes);
  
  for count = 1:numAttributes
    %attrib = item(theAttributes,count-1);
    %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
    %attributes.(attr_name) = char(getValue(attrib));
    
    %Suggestion of Adrian Wanner
    str = toCharArray(toString(item(theAttributes,count-1)))';
    k = strfind(str,'=');
    attr_name = str(1:(k(1)-1));
    attr_name = strrep(attr_name, '-', '_dash_');
    attr_name = strrep(attr_name, ':', '_colon_');
    attr_name = strrep(attr_name, '.', '_dot_');
    attributes.(attr_name) = str((k(1)+2):(end-1));
  end
end
end