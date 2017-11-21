#!/usr/bin/perl

$noparams = ($#ARGV+1==0);

$HTMLdir = '../docs';
$pngdir = $HTMLdir.'/eqn';

$templatefile = "template.tex";

$latexoptions = '--interaction=nonstopmode';
$dvipsoptions = '-o -q';
$pstoimgoptions = '-antialias -scale=1.45 -multipage -crop=a -out=';

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
  open(HANDLE,">".$tmpfile.".tex") or die ("Could not open $tmpfile!");
  print HANDLE $tmp;
  close(HANDLE) or die("Could not close $tmpfile!");

  # generate math graphics files
  #---------------------------------------------------------------
  system('latex '.$latexoptions.' '.$tmpfile.'.tex');
  system('dvips '.$dvipsoptions.' '.$tmpfile);
  system('pstoimg '.$pstoimgoptions.$tmpfile.' '.$tmpfile.'.ps');
  
  # rename all image files to remove _ and leading zeros (e.g. myfile_01.png -> myfile1.png)
  system('rename s/_0/_/ *.png');
  system('rename s/_// *.png');
  
  # move all image files
  system('mv '.$tmpfile.'*.png '.$pngdir);
  # remove all temporary files
  system('\rm '.$tmpfile.'.*');
  
  # insert <img> tags into HTML file
  #---------------------------------------------------------------
  $html =~ s/<!--MATH/<!--IMG--><!--MATH/g;
  for ($k=1; $k<=$n_equations; $k++) {
    $imgtag = '<img src="eqn/'.$tmpfile.$k.'.png" alt="[eqn]">';
    $html =~ s/<\!--IMG-->/$imgtag/s;
  }
  open(HANDLE,'>'.$file) or die("Cannot open $file!");
  print HANDLE $html;
  close(HANDLE) or die("Cannot close $file!");
}
}
