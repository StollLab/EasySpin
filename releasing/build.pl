use strict;
use warnings;

# other dependencies
use Fcntl ':flock'; # for locking on system level
use Net::SSH::Perl; # to use ssh
use Cwd; # to get working directory
my $WorkingDir = getcwd; # get current working directory

# variables imported from config.pl
our ($SourceDir, $BuildsDir, $TempRepoDir, $StableMajorVersion, $DefaultMajorVersion, @VersionCutoff, $esbuild, $KeyBitBucket);

require './config.pl'; # load the configuration file

# settings ------------------------------------------------------------------
# options for File locking - to ensure only one instance is running:
our $NumberOfAttempts = 3; # number of Attempts to obtain a lock
our $WaitTime = 90; # time to wait between Attempts in seconds

# ---------------------------------------------------------------------------------
# creating a lock File 
my $LockFilename = "build.lock";
open (my $LockFile,'>'.$SourceDir.'/'.$LockFilename) or die $!;

my $Attempts = 0;
my $LockObtained = 0;

while ($Attempts < $NumberOfAttempts) {
    $LockObtained = flock $LockFile, LOCK_EX|LOCK_NB;
    if ($LockObtained) {
        last;
    }
    ++$Attempts;
    print("Another instance of build.pl appears to be running, trying again in $WaitTime seconds.\n");
    sleep($WaitTime);
}

if ($LockObtained) {
    print "Created lock File. \n";
}
else {
    print "Can not obtain lock, exiting. Maybe a build got stuck? \n";
    exit;
}

# ---------------------------------------------------------------------------------
# set up build environment
system("ssh-add $KeyBitBucket"); # private key to log into bitbucket, needs to be adapted specific user

# delete and reinitialize temporary directory of EasySpin if a previous build crashed
if (-e "$TempRepoDir") {
    system("rm -r $TempRepoDir");
}
system("mkdir $TempRepoDir");
system(qq(hg clone ssh://hg\@bitbucket.org/sstoll/easyspin $TempRepoDir));

# create the directory where builds are stored if not already available
unless (-e "$BuildsDir") {
    system("mkdir $BuildsDir");
}

# delete temporary build directory where p files are encoded in if a previous build crashed
system("rm -rf /tmp/easyspin*");



# update hg configuration file to include the purge extension, which is needed to cleanly update to a different commint
my $LinesToAdd = qq([extensions]\npurge = );

open(my $hgConf, '>>', "$TempRepoDir.hg/hgrc") or die "Could not open hg config File!";
say $hgConf $LinesToAdd;
close $hgConf;

# -----------------------------------------------------------------
# process tag
my @TagsToBuild = ();

my @NewestVersion = (0, 0, 0);  # (stable, default, dev)

my $callTags = qq(hg tags -R $TempRepoDir); # read tagfile
my @TagFile = `$callTags`;

# compute the numeric value of the cutoff version
my $NumericCutoff = 100000*$VersionCutoff[0]+1000*$VersionCutoff[1]+$VersionCutoff[2];

unless ($ARGV[0]) {
    #  if called without an argument, all missing versions are built

    my @AvailableTags = ();

    # process the hg AvailableTags File and grab the NumericVersion numbers, including the dev NumericVersions, but not 'tip'
    foreach (@TagFile) {
        my @BuildID = ($_ =~ m/(\d+)\.(\d+)\.(\d+)(.*?)\s/);

        if (@BuildID and $BuildID[0]){
            my $NumericVersion = 100000*$BuildID[0]+1000*$BuildID[1]+$BuildID[2];

            if ($NumericVersion >= $NumericCutoff) {
            my $ID = "$BuildID[0].$BuildID[1].$BuildID[2]$BuildID[3]";
            push @AvailableTags,$ID;
            }

        }
    }

    # get highest NumericVersion number for the three branches

    # scan through all the AvailableTags, and compare them to the NewestVersions
    foreach (@AvailableTags) {
        my @BuildID = ($_ =~ m/(\d+).(\d+).(\d+)(.*)/); # match major, minor, patch and everything that follows
        my $NumericVersion = 100000*$BuildID[0]+1000*$BuildID[1]+$BuildID[2];
        
        if ($BuildID[3]){ # check if the currently processed tag is a developer NumericVersion
            my @DevVersion = ($BuildID[3] =~ m/-?([a-z]+).*(\d+)/);
            # read whether tag is alpha beta or dev
            if ($DevVersion[0] eq 'alpha') {
                $NumericVersion = $NumericVersion + 0.2 
            }
            elsif ($DevVersion[0] eq 'beta') {
                $NumericVersion = $NumericVersion + 0.3
            }
            elsif ($DevVersion[0] eq 'dev') {
                $NumericVersion = $NumericVersion + 0.1
            }

            # in case the tag also contains a numeric value, eg. dev.3, alpha1
            if ($DevVersion[1]) {
                $NumericVersion = $NumericVersion + 0.01*$DevVersion[1];
            }

            # Update NewestVersion if necessary
            $NewestVersion[2] = $NumericVersion if $NumericVersion > $NewestVersion[2];
        }
        # if major version corresponds to the default branch:
        elsif ($BuildID[0] eq $DefaultMajorVersion) {
            $NewestVersion[1] = $NumericVersion if $NumericVersion > $NewestVersion[1];
        }
        # if major version corresponds to the stable branch:
        elsif ($BuildID[0] eq $StableMajorVersion) {
            $NewestVersion[0] = $NumericVersion if $NumericVersion > $NewestVersion[0];
        }
    }

    print "The most recent NumericVersions are @NewestVersion \n";

    #  Get the currently available builds in the build directory
    opendir(my $FilesInBuildDir,$BuildsDir);
    my @AvailableBuilds = readdir($FilesInBuildDir);

    my $zipFiles = '';
    for my $File (@AvailableBuilds) {
        my @BuildID = ($File =~ m/(\d+).(\d+).(\d+)(.*).zip/);
        if (@BuildID){
            my $ID = "$BuildID[0].$BuildID[1].$BuildID[2]$BuildID[3],";
            $zipFiles = join( "", $zipFiles, $ID);
        }
    }

    print "The following builds are already in the build directory: $zipFiles \n";

    # Identify which version need to be built
    foreach (@AvailableTags) {
        unless ($zipFiles =~ m/$_,/) {
            push @TagsToBuild, $_;
        }
    }
}
else {
    # if build.pl is called with a commandline argument
    my $cmdLineArgument = $ARGV[0];
    my $TagExists = 0;

    # try to match the commandline argument against the semantic versioning
    my @SemanticBuildID = ($cmdLineArgument =~ m/(\d+)\.(\d+)\.(\d+)(.*?)/);

    # if argument corresponds to semantic versioning, make sure version is newer than cutoff version
    if (@SemanticBuildID){
        my $NumericVersion = 100000*$SemanticBuildID[0]+1000*$SemanticBuildID[1]+$SemanticBuildID[2];

        if ($NumericVersion < $NumericCutoff) {
            die "Only Easyspin versions starting from $VersionCutoff[0].$VersionCutoff[1].$VersionCutoff[2] can be built using this script \n";
        }

    }

    # check if provided tag actually exists in tag file
    foreach (@TagFile) {
        if ($_ =~ m/\b$cmdLineArgument\b/) {
            $TagExists = 1  ;
        }
    }

    # error if tag is not existent
    unless ($TagExists) {
        die "the tag '$cmdLineArgument' was not found \n";
    }
    push @TagsToBuild, $cmdLineArgument;
}

print("The following versions will be built: @TagsToBuild \n");

# ---------------------------------------------------------------------------------
# Processes all the TagsToBuild
my $MatchPattern = '(\d+).(\d+).(\d+)-?([a-z]+)?[-.]?(\d+)?';

# loop over all tags that should be built
foreach (@TagsToBuild) {
    my $thisBuild = $_;
    my @thisBuildID = ($thisBuild =~ m/$MatchPattern/);

    print "Building $thisBuild \n";

    # clean cloned repo
    system('hg purge -R '.$TempRepoDir);

    # update to tag that is being built
    system("hg update $thisBuild -R $TempRepoDir -C");

    # ---------------------------------------------------------------------------------
    # Update html file that contains the examples
    chdir($TempRepoDir.'scripts/');

    system('perl mkexamples.pl');

    chdir($WorkingDir);

    # ---------------------------------------------------------------------------------
    # Build .svg NumericVersions of formulae in the html documentation and replace the AvailableTags in html Files
    print "Converting LaTeX markup in documentation to svg and modifying html Files. \n";
    system("perl mkformulas.pl");

    # ---------------------------------------------------------------------------------
    # Update variables in esbuild.m 
    my $matchReleaseID = "ReleaseID(.*?); \% major.minor.patch"; # pattern to find ReleaseID 
    my $replaceReleaseID = "ReleaseID = \'".$thisBuild."\'; \% major.minor.patch"; # Update ReleaseID
    
    my $matchSourceDir = "SourceDir = (.*?);";
    my $replaceSourceDir = "SourceDir = ['$TempRepoDir'];";

    my $matchZipDestDir = "ZipDestDir = (.*?);";
    my $replaceZipDestDir = "ZipDestDir = ['$BuildsDir'];";

    my $esbuildNew = $esbuild.'new';

    print("Updating ReleaseID in esbuild.m\n");
    open(my $Input,'<'.$esbuild) or die("Cannot open $esbuild!");
    open(my $Output,'>'.$esbuildNew) or die("Cannot open $esbuildNew!");
    while (<$Input>) {
        $_ =~ s/$matchReleaseID/$replaceReleaseID/g;
        $_ =~ s/$matchSourceDir/$replaceSourceDir/g;
        $_ =~ s/$matchZipDestDir/$replaceZipDestDir/g;
        print $Output $_;    
    }

    close($Input) or die("Cannot close $_!");
    close($Output) or die("Cannot close $_!");

    system('mv '.$esbuildNew.' '.$esbuild);

    # ---------------------------------------------------------------------------------
    # Call Matlab to run esbuild.m
    my $MatlabOptions = '-nosplash -nodesktop -nodisplay';
    my $MatlabTarget = qq(-r "run('$esbuild');exit;");

    print("Triggering Matlab build \n");
    system('matlab '.$MatlabOptions." ".$MatlabTarget);

    # ---------------------------------------------------------------------------------
    # Translate semantic version, compare to NewestVersions and decide wether it needs to be uploaded

    my $NumericVersion = 100000*$thisBuildID[0]+1000*$thisBuildID[1]+$thisBuildID[2];

    if ($thisBuildID[3]){ # check if is an developer NumericVersion
        if ($thisBuildID[3] eq 'alpha') {
            $NumericVersion = $NumericVersion + 0.2 
        }
        elsif ($thisBuildID[3] eq 'beta') {
            $NumericVersion = $NumericVersion + 0.3
        }
        elsif ($thisBuildID[3] eq 'dev') {
            $NumericVersion = $NumericVersion + 0.1
        }
        if ($thisBuildID[4]) {
            $NumericVersion = $NumericVersion + 0.01*$thisBuildID[4];
        }
    }

    my $UploadBuild = 0;
    foreach my $VersionToCompare (@NewestVersion) {
        if ($NumericVersion == $VersionToCompare) {
            $UploadBuild = 1;
        }
    }

    # upload the current build by calling publish.pl
    if ($UploadBuild) {

        system(qq(perl publish.pl $thisBuild));
       
    }
}

# ---------------------------------------------------------------------------------
# Clean up temporary EasySpin directories
if (-e "$TempRepoDir") {
    system("rm -r $TempRepoDir");
}

# ---------------------------------------------------------------------------------
# Clean up LockFile and exit
print "removing LockFile \n";
close $LockFile;
system('rm '.$SourceDir.'/'.$LockFilename);

# ---------------------------------------------------------------------------------
print "All finished.\n";