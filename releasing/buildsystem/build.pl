use strict;
use warnings;

# other dependencies
use Fcntl ':flock';  # for locking on system level
use Cwd;  # to get working directory
use IO::Handle;  # for $lockFile->autoflush, so lock-file diagnostics are visible immediately
use FindBin;  # to locate config.pl next to this script, regardless of cwd
my $workingDir = getcwd;  # get current working directory

# variables imported from config.pl
our ($SourceDir, $BuildsDir, $TempRepoDir);  # directories
our ($StableMajorVersion,  $KeyForStableVersion, $MonthsToExpireStable);  # settings for the stable versions
our ($DefaultMajorVersion, $KeyForDefaultVersion, $MonthsToExpireDefault);  # settings for the default versions
our ($KeyForDeveloperVersion, $KeyForExperimentalVersion, $MonthsToExpireDeveloper);  # settings for the experimental versions
our (@VersionCutoff, $KeyGitHub);  # some other settings

print "Reading config file\n";

require "$FindBin::Bin/config.pl";  # load configuration file (always next to this script)

# Config paths may use a leading "~" for the home directory. system() calls
# below are list-form (no shell involved), so nothing expands "~"
# automatically; expand it here with glob(), which resolves "~" without
# shelling out. Falls back to the original string if glob() finds nothing.
if (defined $KeyGitHub) {
    my ($expandedKey) = glob($KeyGitHub);
    $KeyGitHub = $expandedKey if defined $expandedKey;
}

# Settings
# ---------------------------------------------------------------------------------
# Options for file locking - to ensure only one instance is running
my $numberOfAttempts = 3; # number of Attempts to obtain a lock
my $waitTime = 90; # time to wait between Attempts in seconds

# Create a lock file
# ---------------------------------------------------------------------------------
print "Creating lock file \n";
my $lockFilename = "build.lock";
my $lockPath = $SourceDir.'/'.$lockFilename;
open (my $lockFile,'>'.$lockPath) or die "Cannot open lock file '$lockPath': $!\n";
$lockFile->autoflush(1);

my $lockObtained = 0;

for my $attempt (1..$numberOfAttempts) {
    $lockObtained = flock $lockFile, LOCK_EX|LOCK_NB;
    last if $lockObtained;
    print("Another instance of build.pl appears to be running, trying again in $waitTime seconds.\n");
    sleep($waitTime) if $attempt < $numberOfAttempts;
}

if ($lockObtained) {
    # Record who holds the lock and since when, for diagnostics if a future run has to wait on it.
    print $lockFile "Lock file created by pid $$ at ".scalar(localtime)."\n";
    print "Created lock file. \n";
}
else {
    print STDERR "Cannot obtain lock after $numberOfAttempts attempts (another instance appears to be running).\n";
    close $lockFile;
    exit 1;
}

# ---------------------------------------------------------------------------------
# Main build sequence, wrapped in eval so that any failure (die) still falls
# through to the cleanup below instead of leaving the lock file held and the
# temporary EasySpin clone behind.
# ---------------------------------------------------------------------------------
my $anyFailures = 0;  # for tracking non-fatal per-tag failures
my $success = eval {

    # Set up build environment
    # ---------------------------------------------------------------------------------
    my $sshadd_status = system('ssh-add', $KeyGitHub); # private key to log into GitHub, needs to be adapted to specific user
    die "ssh-add failed: exit code  $?\n" if $sshadd_status!=0;

    # Delete and reinitialize temporary directory of EasySpin if a previous build crashed
    if (-e "$TempRepoDir") {
        my $rm_status = system('rm', '-rf', $TempRepoDir);
        die "Could not remove stale temporary directory '$TempRepoDir': exit code  $?\n" if $rm_status!=0;
    }
    my $mkdir_status = system('mkdir', $TempRepoDir);
    die "Could not create temporary directory '$TempRepoDir': exit code  $?\n" if $mkdir_status!=0;
    my $clone_status = system('git', 'clone', 'git@github.com:StollLab/EasySpin.git', $TempRepoDir);
    die "git clone failed: exit code  $?\n" if $clone_status!=0;

    # Create the directory where builds are stored if not already available
    unless (-e "$BuildsDir") {
        my $mkdir_builds_status = system('mkdir', $BuildsDir);
        die "Could not create builds directory '$BuildsDir': exit code  $?\n" if $mkdir_builds_status!=0;
    }

    # Delete temporary build directory where p files are encoded in if a previous build crashed
    # (best-effort cleanup of another run's leftovers; not fatal if it fails)
    my @stalePaths = glob('/tmp/easyspin*');
    if (@stalePaths) {
        system('rm', '-rf', @stalePaths) == 0
            or warn "Could not remove stale /tmp/easyspin* files: exit code $?\n";
    }


    # Process tag
    # -----------------------------------------------------------------
    my @tagsToBuild = ();

    # If a specific tag is requested on the command line, it is always
    # published once built -- @newestVersion is only ever computed in the
    # "no argument" branch below, so the usual newest-in-channel comparison
    # has nothing to compare against for an explicitly requested tag.
    my $forceUpload = 0;

    my @newestVersion = (0, 0, 0);  # (stable, default, dev)

    my $callTags = qq(git --git-dir=$TempRepoDir/.git tag);  # read tagfile
    my @tagFile = `$callTags`;
    die "git tag failed: exit code  $?\n" if $?!=0;
    die "git tag returned no tags; the repository may not have cloned correctly\n" unless @tagFile;

    # Compute the numeric value of the cutoff version
    my $numericCutoff = 100000*$VersionCutoff[0]+1000*$VersionCutoff[1]+$VersionCutoff[2];

    unless ($ARGV[0]) {
        # If called without an argument, all missing versions are built

        my @availableTags = ();

        # Process the availableTags File and grab the numericVersion numbers, including the dev numericVersions, but not 'tip'
        foreach (@tagFile) {
            my @vBuildID = ($_ =~ m/(v?)(\d+)\.(\d+)\.(\d+)(.*?)\s/);

            if (@vBuildID and $vBuildID[1]){
                my $numericVersion = 100000*$vBuildID[1]+1000*$vBuildID[2]+$vBuildID[3];

                if ($numericVersion >= $numericCutoff) {
                    my $id = "$vBuildID[0]$vBuildID[1].$vBuildID[2].$vBuildID[3]$vBuildID[4]";
                    push @availableTags,$id;
                }

            }
        }

        # Get highest numericVersion number for the three branches

        # Scan through all the availableTags, and compare them to the newestVersions
        foreach (@availableTags) {
            my @buildID = ($_ =~ m/v?(\d+).(\d+).(\d+)(.*)/); # match major, minor, patch and everything that follows
            my $numericVersion = 100000*$buildID[0]+1000*$buildID[1]+$buildID[2];

            if ($buildID[3]){ # check if the currently processed tag is a developer numericVersion
                my @devVersion = ($buildID[3] =~ m/-?([a-z]+).*?(\d{1,3})/);
                # read whether tag is alpha beta or dev
                if ($devVersion[0] eq 'alpha') {
                    $numericVersion = $numericVersion + 0.2
                }
                elsif ($devVersion[0] eq 'beta') {
                    $numericVersion = $numericVersion + 0.3
                }
                elsif ($devVersion[0] eq 'dev') {
                    $numericVersion = $numericVersion + 0.1
                }

                # in case the tag also contains a numeric value, eg. dev.3, alpha1
                if ($devVersion[1]) {
                    $numericVersion = $numericVersion + 0.0001*$devVersion[1];
                }

                # Update newestVersion if necessary
                $newestVersion[2] = $numericVersion if $numericVersion > $newestVersion[2];
            }
            # if major version corresponds to the default branch:
            elsif ($buildID[0] eq $DefaultMajorVersion) {
                $newestVersion[1] = $numericVersion if $numericVersion > $newestVersion[1];
            }
            # if major version corresponds to the stable branch:
            elsif ($buildID[0] eq $StableMajorVersion) {
                $newestVersion[0] = $numericVersion if $numericVersion > $newestVersion[0];
            }
        }

        print "The most recent NumericVersions are @newestVersion \n";

        # Get the currently available builds in the build directory
        opendir(my $filesInBuildDir,$BuildsDir);
        my @availableBuilds = readdir($filesInBuildDir);

        my $zipFiles = '';
        for my $file (@availableBuilds) {
            my @buildID = ($file =~ m/(v?)(\d+).(\d+).(\d+)(.*).zip/);
            if (@buildID){
                my $id = "$buildID[0]$buildID[1].$buildID[2].$buildID[3]$buildID[4],";
                $zipFiles = join( "", $zipFiles, $id);
            }
        }

        print "The following builds are already in the build directory: $zipFiles \n";

        # Identify which version need to be built. $zipFiles is built from the
        # actual filenames in $BuildsDir, which never carry a leading "v" (e.g.
        # "easyspin-6.0.14.zip"), while @availableTags keeps whatever prefix the
        # git tag itself uses (e.g. "v6.0.14"); strip it for this comparison only
        # so an existing build isn't mistaken for a missing one. The original,
        # unstripped tag is still what gets pushed to @tagsToBuild and later used
        # for "git checkout", since that's the real tag name known to git.
        foreach (@availableTags) {
            (my $bareTag = $_) =~ s/^v//;
            unless ($zipFiles =~ m/\Q$bareTag\E,/) {
                push @tagsToBuild, $_;
            }
        }
    }
    else {  # If build.pl is called with a commandline argument

        my $cmdLineArgument = $ARGV[0];

        # Reject anything that isn't a plausible tag identifier
        die("Invalid tag identifier '$cmdLineArgument': only letters, digits, dots, dashes and underscores are allowed.\n")
            unless $cmdLineArgument =~ /^[A-Za-z0-9_.-]+$/;

        my $tagExists = 0;

        # Try to match the commandline argument against the semantic versioning
        my @semanticBuildID = ($cmdLineArgument =~ m/v?(\d+)\.(\d+)\.(\d+)(.*?)/);

        # If argument corresponds to semantic versioning, make sure version is newer than cutoff version
        if (@semanticBuildID){
            my $numericVersion = 100000*$semanticBuildID[0]+1000*$semanticBuildID[1]+$semanticBuildID[2];

            if ($numericVersion < $numericCutoff) {
                die "Only Easyspin versions starting from $VersionCutoff[0].$VersionCutoff[1].$VersionCutoff[2] can be built using this script \n";
            }

        }

        # Check if provided tag actually exists in tag file
        foreach (@tagFile) {
            if ($_ =~ m/\b\Q$cmdLineArgument\E\b/) {
                $tagExists = 1;
            }
        }

        # Error if tag is not existent
        unless ($tagExists) {
            die "The tag '$cmdLineArgument' was not found \n";
        }
        push @tagsToBuild, $cmdLineArgument;
        $forceUpload = 1;
    }

    if (@tagsToBuild==0) {
        print("No new versions to build. \n");
    }
    else {
        print("The following versions will now be built: @tagsToBuild \n");
    }


    # Loop over all tags that should be built
    # ---------------------------------------------------------------------------------
    foreach (@tagsToBuild) {
        my $thisTag = $_;

        print "Building tag $thisTag \n";

        # Update git repository
        # ---------------------------------------------------------------------------------
        # Clean cloned repo
        my $gitclean_status = system('git', "--git-dir=$TempRepoDir/.git", "--work-tree=$TempRepoDir", 'clean', '-f');
        die "git clean failed for $thisTag: exit code  $?\n" if $gitclean_status!=0;

        # Update to tag that is being built
        my $checkout_status = system('git', "--git-dir=$TempRepoDir/.git", "--work-tree=$TempRepoDir", '-c', 'advise.detachedHead=false', 'checkout', $thisTag);
        die "git checkout failed for $thisTag: exit code  $?\n" if $checkout_status!=0;

        # Generate HTML file that contains list of examples
        # ---------------------------------------------------------------------------------
        chdir($TempRepoDir.'releasing/');
        print "Creating list with examples. \n";
        my $mkexamples_status = system('perl', 'mkexamples.pl');
        chdir($workingDir);
        die "mkexamples.pl failed for $thisTag: exit code  $?\n" if $mkexamples_status!=0;

        # Build documentation (compile math with LaTeX etc.)
        # ---------------------------------------------------------------------------------
        print "Building documentation.\n";
        my $docbuilder_status = system('perl', 'docbuilder.pl');
        die "docbuilder.pl failed for $thisTag: exit code  $?\n" if $docbuilder_status!=0;

        # Write esbuild_config.m
        # ---------------------------------------------------------------------------------
        my $vMatchPattern = '(v)(\d+.\d+.*)';
        my $releaseID;

        my @shortTag = ($thisTag =~ m/$vMatchPattern/);

        if ($shortTag[1]) {
            $releaseID = $shortTag[1];
        }
        else {
            $releaseID = $thisTag;
        }

        my $releaseChannel;
        my $monthsToExpire;
        my $versionSlot;  # index into @newestVersion (0=stable, 1=default, 2=dev) this tag belongs to

        my $matchVersionPattern = '(v?)(\d+).(\d+).(\d+)-?([a-z]+)?[-.]?(\d+)?';
        my @thisBuildID = ($thisTag =~ m/$matchVersionPattern/);
        if ($thisBuildID[1]) {
            if ($thisBuildID[4]) {
                $releaseChannel = $KeyForDeveloperVersion;
                $monthsToExpire = $MonthsToExpireDeveloper;
                $versionSlot = 2;
            }
            elsif ($thisBuildID[1] eq $StableMajorVersion) {
                $releaseChannel = $KeyForStableVersion;
                $monthsToExpire = $MonthsToExpireStable;
                $versionSlot = 0;
            }
            elsif ($thisBuildID[1] eq $DefaultMajorVersion) {
                $releaseChannel = $KeyForDefaultVersion;
                $monthsToExpire = $MonthsToExpireDefault;
                $versionSlot = 1;
            }
        }
        else {  # if tag does not follow semantic versioning, e.g. easyspin-evolve.zip
            $releaseChannel = $KeyForExperimentalVersion;
            $monthsToExpire = $MonthsToExpireDeveloper;
        }

        print("Writing esbuild config file\n");
        print "  ReleaseID:       $releaseID\n";
        print "  Release channel: $releaseChannel\n";
        print "  Months:          $monthsToExpire\n";
        print "  Source dir:      $TempRepoDir\n";
        print "  Zip dest. dir:   $BuildsDir\n";

        my $esbuildconfigfile = "esbuild_config.m";
        open(my $output,'>'.$esbuildconfigfile) or die("Cannot open $esbuildconfigfile!");
        print $output "releaseID = '$releaseID';\n";
        print $output "releaseChannel = '$releaseChannel';\n";
        print $output "monthsToExpiry = $monthsToExpire;\n";
        print $output "sourceDir = '$TempRepoDir';\n";
        print $output "zipDestDir = '$BuildsDir';\n";
        close($output) or die("Cannot close $esbuildconfigfile!");


        # Call Matlab to run esbuild.m to generates a packaged zip file
        # ---------------------------------------------------------------------------------
        print("Running MATLAB build\n");
        my @matlabOptions = ('-nosplash', '-nodesktop', '-nodisplay');
        my $matlabCommand = "run('esbuild.m');exit;";
        my $matlab_status = system('matlab', @matlabOptions, '-r', $matlabCommand);
        system('rm', $esbuildconfigfile)==0
            or warn "Could not remove '$esbuildconfigfile': exit code $?\n";
        if ($matlab_status!=0) {
            die "MATLAB build failed for $thisTag: exit code  $matlab_status\n" ;
        }
        print("MATLAB build completed\n");


        # Publish zip file
        # ---------------------------------------------------------------------------------
        # Decide wether to upload version. A tag explicitly requested on the
        # command line is always uploaded; otherwise, the first conditional
        # checks if the build version follows semantic versioning or not
        # (e.g. easyspin-evolve), and only the newest version in its channel
        # gets uploaded.
        my $uploadBuild = 0;
        if ($forceUpload) {
            $uploadBuild = 1;
        }
        elsif ($thisBuildID[1] and defined $versionSlot) {
            # Translate semantic version and compare to the newest version already
            # found for this tag's own channel (not all three), to decide wether it
            # needs to be uploaded.
            my $numericVersion = 100000*$thisBuildID[1]+1000*$thisBuildID[2]+$thisBuildID[3];

            if ($thisBuildID[4]){ # check if is an developer NumericVersion
                if ($thisBuildID[4] eq 'alpha') {
                    $numericVersion = $numericVersion + 0.2
                }
                elsif ($thisBuildID[4] eq 'beta') {
                    $numericVersion = $numericVersion + 0.3
                }
                elsif ($thisBuildID[4] eq 'dev') {
                    $numericVersion = $numericVersion + 0.1
                }
                if ($thisBuildID[5]) {
                    $numericVersion = $numericVersion + 0.0001*$thisBuildID[5];
                }
            }
            $uploadBuild = 1 if $numericVersion == $newestVersion[$versionSlot];
        }

        # Upload the current build by calling publish.pl. A failure here does
        # not abort the rest of the batch -- other tags in @tagsToBuild are
        # independent of this one -- but it is tracked in $anyFailures so the
        # final exit code still reflects that something went wrong.
        if ($uploadBuild) {
            print "Calling publish script for upload.\n";
            my $publish_status = system('perl', 'publish.pl', $thisTag);
            if ($publish_status != 0) {
                warn "publish.pl failed for $thisTag: exit code $publish_status\n";
                $anyFailures = 1;
            }
        }
        else {
            print "No upload.\n";
        }
    }

    1;  # signal success to eval
};
my $error = $@;

if (!$success) {
    print STDERR "ERROR: $error";
}

# Clean up temporary EasySpin directory whether or not the build succeeded
# ---------------------------------------------------------------------------------
if (-e "$TempRepoDir") {
    system('rm', '-rf', $TempRepoDir) == 0
        or warn "Could not remove temporary directory '$TempRepoDir': exit code $?\n";
}

# Remove lock file and exit
# ---------------------------------------------------------------------------------
print "removing lock file \n";
close $lockFile;
unlink($lockPath) or warn "Could not remove lock file: $!\n";

# ---------------------------------------------------------------------------------
if ($success && !$anyFailures) {
    print "Finished.\n";
    exit 0;
}
elsif ($success) {
    print "Finished, but one or more publish steps failed; see warnings above.\n";
    exit 1;
}
else {
    print "Finished with errors.\n";
    exit 1;
}
