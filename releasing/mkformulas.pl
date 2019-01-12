$noparams = ($#ARGV+1==0);  

our ($repodir);

require './config.pl';

$HTMLdir = "$repodir/docs";
$pngdir = $HTMLdir.'/eqn';
$tempdir = "./latextemp/";

if (-e "$tempdir") {
    system("rm -R $tempdir");
}

system("mkdir $tempdir");

$templatefile = "$repodir/scripts/template.tex";

$latexoptions = '--interaction=nonstopmode --output-directory=./latextemp/';
# convert the first 100 pages of .dvi file to .svg and scale by a factor 1.3
$dvisvgmoptions = '--page=-10000 --transform=S1.3 --output=./latextemp/%f%1p';

$templateposition = '%formulaposition';

if ($noparams) {
  opendir(dirr,$HTMLdir) or die("Can't open $HTMLdir: $!\n");
  @names = readdir(dirr) or die("Can't read $HTMLdir: $!\n");
  closedir(dirr);
  @names = grep(/.html$/,sort(@names));
}
else {
  @names = @ARGV;
  #print @names;
}

foreach $htmlfile (@names) {

# read in LaTeX template file
#---------------------------------------------------------------
open templatefile or die("Can't open $templatefile: $!\n");
$templ = join('',<templatefile>);
close templatefile or die("Can't close $templatefile: $!\n");

# read in HTML file
#---------------------------------------------------------------
$file = $HTMLdir.'/'.$htmlfile;
open file or die("Cannot open $file!");
$html = join('',<file>);
close file;
#print $html;

# remove image tags from HTML file
#---------------------------------------------------------------
$imgkey = '<img[^\n]*?"\[eqn\]">';
$html =~ s/$imgkey//sg;

# extract latex snippets from HTML file and write to LaTeX file
#---------------------------------------------------------------
@latexcode = $html =~ /<\!--MATH(.*?)-->/sg;
$n_equations = scalar(@latexcode);


if ($n_equations>0) {
  print "=====================================================\n";
  print "  ".$htmlfile.":  ".$n_equations." latex expressions\n";
  print "=====================================================\n";
  
  $latex = join("\n\\newpage\n",@latexcode);

  $tmp = $templ;
  $tmp =~ s/$templateposition/$latex/;

  $tmpfile = $htmlfile;
  $tmpfile =~ s/\.html//;
  open(HANDLE,"> $tempdir$tmpfile.tex") or die ("Could not open $tmpfile!");
  print HANDLE $tmp;
  close(HANDLE) or die("Could not close $tmpfile!");

  # generate math graphics files
  #---------------------------------------------------------------
  system("latex $latexoptions $tempdir$tmpfile.tex");
  system("dvisvgm $dvisvgmoptions $tempdir$tmpfile.dvi");

  # move all image files
  system("mv  $tempdir$tmpfile*.svg $pngdir");
  # remove all temporary files
  system("rm $tempdir$tmpfile.*");
  
  # insert <img> tags into HTML file
  #---------------------------------------------------------------
  $html =~ s/<!--MATH/<!--IMG--><!--MATH/g;
  for ($k=1; $k<=$n_equations; $k++) {
    $imgtag = '<img src="eqn/'.$tmpfile.$k.'.svg" alt="[eqn]">';
    $html =~ s/<\!--IMG-->/$imgtag/s;
  }
  open(HANDLE,'>'.$file) or die("Cannot open $file!");
  print HANDLE $html;
  close(HANDLE) or die("Cannot close $file!");
}
}

if (-e "$tempdir") {
    system("rm -R $tempdir");
}