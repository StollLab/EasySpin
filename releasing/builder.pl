#!/usr/bin/perl

# lock for lock file, continue if file is too old
# wait 1 minute, check again
# create lock file 
# 

# gets arguments this script was called with
my $branch = @ARGV[0];

# directories
$topdir = './..'; # topdirectory 
$builddir = $topdir.'/easyspin-builds/'; # directory that builds are created in
$uploaddir =  $topdir.'/upload/'; # directory that will be used to upload files to easyspin.org
$serverdir = '~/public_html/easyspin/'; # parentfolder on the server, used to copy files from $uploaddir into

$indexhtml = 'index.html'; # name of the main html file
$downloadthml = 'download.html'; # name of the secondary html file

$pathtoindexhtml = $uploaddir.$indexhtml;
$pathtooldindexhtml = $uploaddir.$indexhtml.'.bak'; 

$pathtodownloadtml = $uploaddir.$downloadthml;
$patholddownloadhtml = $uploaddir.$downloadthml.'.bak'; 

# options
$allowedbranches = qw(stable,default,buildsys);

$matlaboptions = '-nosplash -nodesktop -nodisplay -r';
$matlabrunfile = qq( "run('./releasing/esbuild.m');exit;");

$login = 'easyspin@easyspin.org';

# -----------------------------------------------------------------
# updating and pulling from bitbucket.com 
print "Updating repository from bitbucket.com\n";
# system('hg pull -u');

print "Switching and updating to $branch\n";
# system('hg update '.$branch.);

my $tag = `hg log -r "." --template "{latesttag}"`;
print "Current tag is: $tag \n";

# -----------------------------------------------------------------
# building .svg versions of formulae
print "Converting LaTeX markup to svg \n";
# system('perl ./releasing/mkformulas.pl');

# -----------------------------------------------------------------
# Update tag in m file
# needs to be done
#
#---------------------
# calling Matlab and waiting until a new build was added to the build directory
# Open build dirrectory and get current number of files
opendir(dirr,$builddir) or die("Can't open $builddir: $!\n");
my $entriesinbuilddir = () = readdir(dirr);
my $currentnumberfiles = $entriesinbuilddir;
closedir(dirr);

print("Triggering Matlab build \n");
system('matlab.exe '.$matlaboptions.$matlabrunfile);

print("Waiting for Matlab to finish building and zipping Easyspin\n");

while ($currentnumberfiles == $entriesinbuilddir) {
    opendir(dirr,$builddir);
    $currentnumberfiles = () = readdir(dirr);
    sleep 5;
}
sleep 10;
# pgrep check if matlab is still running
closedir(dirr);

# renaming .zip file and moving to upload directory
my $filename = 'easyspin-'.$tag.'.zip';
print($filename." is moved to upload directory \n");
system('cp '.$builddir.$filename.' '.$uploaddir.$filename);

# -----------------------------------------------------------------
# getting html files from easyspin.org and updating them and delete old zip file

# get index.html from easyspin.org
print("Getting index.html\n");
system('scp '.$login.':'.$serverdir.$indexhtml.' '.$patholdindexhtml);

# regexp to find versions and corresponding zip files in the html files
$oldzipfile = '<!--'.$branch.'-->(.*?)<!--zip-->';
$newzipfile = '<!--'.$branch.'--><a href="easyspin-'.$tag.'.zip"><!--zip-->';

$oldversion = '<!--'.$branch.'-->(.*?)<!--version-->';
$newversion = '<!--'.$branch.'-->'.$tag.'<!--version-->';

# update the index.html file
print("Updating index.html\n");
open(input,'<'.$pathtooldindexhtml) or die("Cannot open $pathtooldindexhtml!");
open(output,'>'.$pathtoindexhtml) or die("Cannot open $pathtoindexhtml!");
while (<input>) {
    $_ =~ s/$oldzipfile/$newzipfile/g;
    $_ =~ s/$oldversion/$newversion/g;
    print output $_;    
}

close(input) or die("Cannot close $_!");
close(output) or die("Cannot close $_!");

# update the download.html file if stable branch is pushed
if ($branch eq 'stable') {
    print("Updating download.html \n");
    system('scp '.$login.':'.$serverdir.$downloadthml.' '.$patholddownloadhtml);
    open(input,'<'.$patholddownloadhtml) or die("Cannot open $_!");
    open(output,'>'.$pathtodownloadtml) or die("Cannot open $_!");
    while (<input>) {
        $_ =~ s/$oldzipfile/$newzipfile/g;
        $_ =~ s/$oldversion/$newversion/g;
        print output $_; 
    }       
    close(input) or die("Cannot close $_!");
    close(output) or die("Cannot close $_!");   
}

# removing backup files
print("deleting backup versions of html files \n");
system('rm '.$uploaddir.'*.bak');


# -----------------------------------------------------------------
# upload entire upload directory to easyspin org and then clean it
print("Uploading all new files to easyspin.org \n");
system('scp '.$uploaddir.'* '.$login.':'.$serverdir.'test/');

# clear upload directory
system("rm ".$uploaddir."*");

# -----------------------------------------------------------------
# ssh delete examples and documentation
# unzip new file on the server, copy documentation and examples, delete unzipped files


# delete lockfile
print "Finished! \n";