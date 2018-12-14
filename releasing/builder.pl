#!/usr/bin/perl

my $branch = @ARGV[0];

$topdir = './..';
$builddir = $topdir.'/easyspin-builds';

opendir(dirr,$builddir) or die("Can't open $builddir: $!\n");

my $entriesinbuilddir = () = readdir(dirr);
closedir(dirr);

$allowedbranches = qw(stable,default,buildsys);

print "Pulling from bitbucket.com\n";
system('hg pull -u');

print "Switching and updating to $branch\n";
# system('hg update '.$branch.);

print "Converting LaTeX markup to svg \n";
# system('perl ./releasing/mkformulas.pl');

my $output = `hg log -r "." --template "{latesttag}"`;

print "Adding tag $output \n";

$matlaboptions = '-nosplash -nodesktop -nodisplay -r';
$matlabrunfile = qq( "run('./releasing/esbuild.m');exit;");

print("Triggering Matlab build \n");
system('matlab.exe '.$matlaboptions.$matlabrunfile);

# Open build dirrectory and get current number of files
opendir(dirr,$builddir) or die("Can't open $builddir: $!\n");
my $currentnumberfiles = () = readdir(dirr);
closedir(dirr);

print("Waiting for Matlab to finish building and zipping Easyspin\n");
# wait until a file gets added
while ($currentnumberfiles == $entriesinbuilddir) {
    opendir(dirr,$builddir);
    $currentnumberfiles = () = readdir(dirr);
    sleep 5;
}

closedir(dirr);

my $filename = 'easyspin-'.$output.'.zip';
print($filename." is being uploaded to easyspin.org\n");
system('mv '.$builddir.'/'.$filename.' '.$builddir.'/test.zip');

$serverdir = '~/public_html/easyspin/';
$login = 'easyspin@easyspin.org';

system('scp '.$builddir.'/test.zip '.$login.':'.$serverdir.'test2.zip');

print "Finished! \n";