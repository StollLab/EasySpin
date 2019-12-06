$noparams = ($#ARGV+1==0);  

use Cwd; # to get working directory
use File::Spec::Functions;

my $WorkingDir = getcwd; # get current working directory

$isWindows = "$^O" eq "MSWin32";

$TempRepoDir = "..";

# store all the paths that are needed to properly create the documentation
$sourcedir = File::Spec->catdir($TempRepoDir, "docsrc");
$targetdir = File::Spec->catdir($TempRepoDir, "documentation");
$pngdir    = File::Spec->catdir($targetdir, "eqn");
$scriptdir = File::Spec->catdir($TempRepoDir, "scripts");
$tempdir   = File::Spec->catdir(".", "latextemp");

# check if documentation folder already exists and delete if yes
if (-e "$targetdir") {
  if ($isWindows) {
    system("rmdir /s /q $targetdir ");
  }
  else{
    system("rm -R $targetdir");
  }
}

# delete the temporary latex directory if it exists
if (-e "$tempdir") {
  if ($isWindows) {
    system("rmdir /s /q $tempdir ");
  }
  else {
    system("rm -R $tempdir");
  }  
}

# copy doc source to documentation and delete html templates folder
$htmltemplatedir = File::Spec->catdir($targetdir, "templates");
if ($isWindows) {
    system("xcopy /E $sourcedir $targetdir\\ > NUL"); # double slash \\ is required so that Windows understand it is a directory
    system("rmdir /s /q $htmltemplatedir ");
  }
else{
    system("cp -r $sourcedir $targetdir");
    system("rm -R $htmltemplatedir");
}

system("mkdir $pngdir");
system("mkdir $tempdir");

if (-e $scriptdir) {
    # this catches the legacy version from before releasing and scripts folders where merged
    $templatefile = File::Spec->catfile($TempRepoDir, "scripts", "template.tex");
}
else {
    # this should be the default behavior, after merge of releasing and scripts folder
    $templatefile = File::Spec->catfile($TempRepoDir, "releasing", "template.tex");
}

# set up options for latex and related commands
$latexoptions = "--interaction=batchmode --output-directory=$tempdir";

$dvipsoptions = '-q';
$pstoimgoptions = '-antialias -scale=1.45 -quiet -multipage -crop=a';
$dvipngoptions = "-q* -T tight -D 150 -o";

$templateposition = '%formulaposition';

if ($noparams) {
  opendir(dirr,$sourcedir) or die("Can't open $sourcedir: $!\n");
  @names = readdir(dirr) or die("Can't read $sourcedir: $!\n");
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
    $file = File::Spec->catfile($sourcedir, $htmlfile);
    open file or die("Cannot open $file!");
    $html = join('',<file>);
    close file;
    # print $html;

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
      $htmlpng = $tmpfile;
      $tmpfile = File::Spec->catfile($tempdir, $tmpfile);

      $tmptex = join(".",$tmpfile,"tex");
      $tmpps  = join(".",$tmpfile,"ps");
      $tmpdvi = join(".",$tmpfile,"dvi");

      open(HANDLE,"> $tmptex") or die ("Could not open $tmptex!");
      print HANDLE $tmp;
      close(HANDLE) or die("Could not close $tmptex!");

      # generate math graphics files
      #---------------------------------------------------------------
      if ($isWindows) {
        system("latex $latexoptions $tmptex > NUL");
      }
      else {
        system("latex $latexoptions $tmptex >/dev/null");
      }
        
      if ($isWindows){
        $target = join("",$tmpfile,"%d.png");
        system("dvipng $dvipngoptions$target $tmpdvi > NUL");
      }
      else {
        system("dvips  $dvipsoptions -o $tmpps $tmpdvi");
        
        system("pstoimg $pstoimgoptions $tmpps");

        chdir($tempdir);
        # # rename all image files to remove _ and leading zeros (e.g. myfile_01.png -> myfile1.png)
        system('rename s/_0/_/ *.png');
        system('rename s/_// *.png');

        chdir($WorkingDir);
      }
      
      # move all images and remove temporary files
      if ($isWindows){
        system("move $tmpfile*.png $pngdir > NUL");
        system("del $tmpfile.*");

      }
      else {
        system("mv $tmpfile*.png $pngdir");
        system("rm $tmpfile.*");
      }

      # insert <img> tags into HTML file
      # ---------------------------------------------------------------
      $html =~ s/<!--MATH/<!--IMG--><!--MATH/g;
      for ($k=1; $k<=$n_equations; $k++) {
          $imgtag = '<img src="eqn/'.$htmlpng.$k.'.png" alt="[eqn]">';
          $html =~ s/<\!--IMG-->/$imgtag/s;
      }
      open(HANDLE,'>'.$file) or die("Cannot open $file!");
      print HANDLE $html;
      close(HANDLE) or die("Cannot close $file!");
    }
}

if (-e "$tempdir") {
  if ($isWindows) {
    system("rmdir /s /q $tempdir ");
  }
  else {
    system("rm -R $tempdir");
  }  
}