use strict;
use warnings;

# other dependencies
use Fcntl ':flock'; # for locking on system level
use Net::SSH::Perl; # to use ssh
use Cwd;
my $workingdir = getcwd;

# values imported from config.pl
our ($parentdir, $builddir, $uploaddir, $serverdir, $repodir, $stableversion, $defaultversion,  @cutoffversion, $esbuild, $username, $hostname, @htmlfiles);

require './config.pl';

# settings ------------------------------------------------------------------
# options for file locking - to ensure only one instance is running:
our $numberallowedattempts = 3; # number of attempts to obtain a lock
our $waittime = 90; # time to wait between attempts in seconds

# ---------------------------------------------------------------------------------
# creating a lock file 
my $lockfilename = "build.lock";
open (my $lockfile,'>'.$parentdir.'/'.$lockfilename) or die $!;

my $attempts = 0;
my $lockobtained = 0;

while ($attempts < $numberallowedattempts) {
    $lockobtained = flock $lockfile, LOCK_EX|LOCK_NB;
    if ($lockobtained) {
        last;
    }
    ++$attempts;
    print("Another instance of build.pl appears to be running, trying again in $waittime seconds.\n");
    sleep($waittime);
}

if ($lockobtained) {
    print "Created lock file. \n";
}
else {
    print "Can not obtain lock, exiting. Maybe a build got stuck? \n";
    exit;
}

system("ssh-add ~/.ssh/no_phrase_rsa"); # private key to log into bitbucket, needs to be adapted specific user

if (-e "$repodir") {
    system("rm -R $repodir");
}

system("mkdir $repodir");
system(qq(hg clone ssh://hg\@bitbucket.org/sstoll/easyspin $repodir));


unless (-e "$builddir") {
    system("mkdir $builddir");
}

my $LinesToAdd = qq([extensions]\npurge = );

open(my $hgConf, '>>', "$repodir.hg/hgrc") or die "Could not open hg config file!";
say $hgConf $LinesToAdd;
close $hgConf;

# -----------------------------------------------------------------
my @tagstobuild = ();

my @newestversion = (0, 0, 0);  # (stable, default, dev)

my $callTags = qq(hg tags -R $repodir);
my @tagfile = `$callTags`;

my $NumericCutoff = 100000*$cutoffversion[0]+1000*$cutoffversion[1]+$cutoffversion[2];

unless ($ARGV[0]) {

    my @tags = ();

    # process the hg tags file and grab the version numbers, including the dev versions, but not 'tip'
    foreach (@tagfile) {
        my @buildID = ($_ =~ m/(\d+)\.(\d+)\.(\d+)(.*?)\s/);

        if (@buildID and $buildID[0]){
            my $version = 100000*$buildID[0]+1000*$buildID[1]+$buildID[2];

            if ($version >= $NumericCutoff) {
            my $ID = "$buildID[0].$buildID[1].$buildID[2]$buildID[3]";
            push @tags,$ID;
            }

        }
    }

    # get highest version number for the three branches

    # scan through all the tags, and compare them to the newestversions
    foreach (@tags) {
        my @buildID = ($_ =~ m/(\d+).(\d+).(\d+)(.*)/); # match major, minor, patch and everything that follows
        my $version = 100000*$buildID[0]+1000*$buildID[1]+$buildID[2];
        
        if ($buildID[3]){ # check if is an developer version
            my @devVersion = ($buildID[3] =~ m/-?([a-z]+).*(\d+)/);
            if ($devVersion[0] eq 'alpha') {
                $version = $version + 0.2 
            }
            elsif ($devVersion[0] eq 'beta') {
                $version = $version + 0.3
            }
            elsif ($devVersion[0] eq 'dev') {
                $version = $version + 0.1
            }

            if ($devVersion[1]) {
                $version = $version + 0.01*$devVersion[1];
            }

            $newestversion[2] = $version if $version > $newestversion[2];
        }
        elsif ($buildID[0] eq $defaultversion) {
            $newestversion[1] = $version if $version > $newestversion[1];
        }
        elsif ($buildID[0] eq $stableversion) {
            $newestversion[0] = $version if $version > $newestversion[0];
        }
    }

    print "The most recent versions are @newestversion \n";

    #  Get the currently available builds in the build directory
    opendir(my $buildfiles,$builddir);
    my @availablebuilds = readdir($buildfiles);

    my $zipfiles = '';
    for my $file (@availablebuilds) {
        my @buildID = ($file =~ m/(\d+).(\d+).(\d+)(.*).zip/);
        if (@buildID){
            my $ID = "$buildID[0].$buildID[1].$buildID[2]$buildID[3],";
            $zipfiles = join( "", $zipfiles, $ID);
        }
    }

    print "The following builds are already in the build directory: $zipfiles \n";

    # Identify which versions need to be built
    foreach (@tags) {
        unless ($zipfiles =~ m/$_,/) {
            push @tagstobuild, $_;
        }
    }
}
else {
    my $cmdlinearg = $ARGV[0];
    my $tagexists = 0;

    my @SemanticBuildID = ($cmdlinearg =~ m/(\d+)\.(\d+)\.(\d+)(.*?)/);
    print @SemanticBuildID;

    if (@SemanticBuildID){
        my $version = 100000*$SemanticBuildID[0]+1000*$SemanticBuildID[1]+$SemanticBuildID[2];

        if ($version < $NumericCutoff) {
            die "Only Easyspin versions starting from $cutoffversion[0].$cutoffversion[1].$cutoffversion[2] can be built using this script \n";
        }

    }

    foreach (@tagfile) {
        if ($_ =~ m/\b$cmdlinearg\b/) {
            $tagexists = 1  ;
        }
    }

    unless ($tagexists) {
        die "the tag '$cmdlinearg' was not found \n";
    }
    push @tagstobuild, $cmdlinearg;
}

print("The following versions will be built: @tagstobuild \n");

# ---------------------------------------------------------------------------------
# Processes all the tagstobuild
# Uploads happen only if the newest version of a release channel is being built
my $MatchPattern = '(\d+).(\d+).(\d+)-?([a-z]+)?[-.]?(\d+)?';

foreach (@tagstobuild) {
    my $thisBuild = $_;
    my @thisBuildID = ($thisBuild =~ m/$MatchPattern/);

    print "Building $thisBuild \n";

    system('hg purge -R '.$repodir);

    system("hg update $thisBuild -R $repodir -C");

    # ---------------------------------------------------------------------------------
    # Build .svg versions of formulae in the html documentation and replace the tags in html files
    chdir($repodir.'scripts/');

    system('perl mkexamples.pl');

    chdir($workingdir);

    # ---------------------------------------------------------------------------------
    # Build .svg versions of formulae in the html documentation and replace the tags in html files
    print "Converting LaTeX markup in documentation to svg and modifying html files. \n";
    system("perl mkformulas.pl");

    # ---------------------------------------------------------------------------------
    # Update ReleaseID (and betaVersion) in esbuild.m 
    my $findReleaseID = "ReleaseID(.*?); \% major.minor.patch"; # pattern to find ReleaseID 
    my $thisReleaseID = "ReleaseID = \'".$thisBuild."\'; \% major.minor.patch"; # Update ReleaseID
    
    my $findSourceDir = "SourceDir = (.*?);";
    my $replaceSourceDir = "SourceDir = ['$repodir'];";

    my $findZipDestDir = "ZipDestDir = (.*?);";
    my $replaceZipDestDir = "ZipDestDir = ['$builddir'];";

    my $esbuildNew = $esbuild.'new';

    print("Updating ReleaseID in esbuild.m\n");
    open(my $input,'<'.$esbuild) or die("Cannot open $esbuild!");
    open(my $output,'>'.$esbuildNew) or die("Cannot open $esbuildNew!");
    while (<$input>) {
        $_ =~ s/$findReleaseID/$thisReleaseID/g;
        $_ =~ s/$findSourceDir/$replaceSourceDir/g;
        $_ =~ s/$findZipDestDir/$replaceZipDestDir/g;
        print $output $_;    
    }

    close($input) or die("Cannot close $_!");
    close($output) or die("Cannot close $_!");

    system('mv '.$esbuildNew.' '.$esbuild);

    # ---------------------------------------------------------------------------------
    # Needs to be reworked - use pgrep call to Matlab does not work on my virtual machine linux
    # ---------------------------------------------------------------------------------
    # Call Matlab to run esbuild.m
    my $MatlabOptions = '-nosplash -nodesktop -nodisplay';
    my $MatlabTarget = qq(-r "run('$esbuild');exit;");

    # Open build dirrectory and get current number of files
    # opendir(my $tempdir,$builddir) or die("Can't open $builddir: $!\n");

    # my $EntriesInBuilddir = () = readdir($tempdir);
    # my $CurrentEntriesInBuilddir = $EntriesInBuilddir;

    # closedir($tempdir);

    print("Triggering Matlab build \n");
    system('matlab '.$MatlabOptions." ".$MatlabTarget);

    # print("Waiting for Matlab to finish building EasySpin $thisBuild \n");

    # while ($CurrentEntriesInBuilddir == $EntriesInBuilddir) {
    #     opendir(my $tempdir,$builddir);
    #     $CurrentEntriesInBuilddir = () = readdir($tempdir);
    #     closedir($tempdir);
    #     sleep 5;
    # }
    # sleep 10;

    # ---------------------------------------------------------------------------------
    # Translate semantic versioning and compare decived wether it needs to be uploaded

    my $version = 100000*$thisBuildID[0]+1000*$thisBuildID[1]+$thisBuildID[2];

    if ($thisBuildID[3]){ # check if is an developer version
        if ($thisBuildID[3] eq 'alpha') {
            $version = $version + 0.2 
        }
        elsif ($thisBuildID[3] eq 'beta') {
            $version = $version + 0.3
        }
        elsif ($thisBuildID[3] eq 'dev') {
            $version = $version + 0.1
        }
        if ($thisBuildID[4]) {
            $version = $version + 0.01*$thisBuildID[4];
        }
    }

    my $upload = 0;
    foreach my $VersionToCompare (@newestversion) {
        if ($version == $VersionToCompare) {
            $upload = 1;
        }
    }

    if ($upload) {

        system(qq(perl publish.pl $thisBuild));
       
    }
}


#-------------------------------
if (-e "$repodir") {
    system("rm -R $repodir");
}
# ---------------------------------------------------------------------------------
# Clean up lockfile and exit
print "removing lockfile \n";
close $lockfile;
system('rm '.$parentdir.'/'.$lockfilename);

print "All finished.\n";