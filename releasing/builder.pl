#!/usr/bin/perl

my $branch = @ARGV[0];

print "Pulling from bitbucket.com\n";
# system('hg pull -u');

print "Switching and updating to '$branch'\n";
# system('hg update '.$branch.);

print "Converting LaTeX markup to svg \n";
system('perl ./releasing/mkformulas.pl');

my $output = `hg log -r "." --template "{latesttag}\n"`;

print "Adding tag $output \n";

$matlaboptions = '-nosplash -nodesktop -r';
$matlabrunfile = qq( "run('./releasing/esbuild.m');exit;");

# system('matlab.exe '.$matlaboptions.$matlabrunfile);

# wait or while to check while matlab is compiling

# login to server
# upload new files
# what about online documentation

print "Finished! \n";