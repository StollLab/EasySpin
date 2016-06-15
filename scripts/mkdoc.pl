#!/usr/bin/perl

$noparams = ($#ARGV+1==0);

$docdir = '../documentation';

$templatefile = $docdir.'/'.'templates/headertemplate.html_';

if ($noparams) {
  opendir(dirr,$docdir) or die("Can't open $docdir: $!\n");
  @names = readdir(dirr) or die("Can't read $docdir: $!\n");
  closedir(dirr);
  @names = grep(/.html$/,sort(@names));
}
else {
  @names = @ARGV;
}
#print @names;

# read in header template file
open templatefile or die("Can't open $templatefile: $!\n");
$templatehtml = join('',<templatefile>);
close templatefile or die("Can't close $templatefile: $!\n");

# define match pattern and replacement string
$oldhtml = '\<header\>.*\<\/header\>';
$newhtml = $templatehtml;

print $oldhtml;
print $newhtml;

foreach $htmlfile (@names) {

  $file = $docdir.'/'.$htmlfile;
  print "processing ".$file."\n";

  # read in HTML file
  open file or die("Cannot open $file!");
  $html = join('',<file>);
  close file;
    
  # match and replace HTML code
  $html =~ s/$oldhtml/$newhtml/sg;
  
  # store modified HTML
  open(HANDLE,'>'.$file) or die("Cannot open $file!");
  print HANDLE $html;
  close(HANDLE) or die("Cannot close $file!");
}
