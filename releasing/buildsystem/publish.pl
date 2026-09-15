use strict;
use warnings;

# other dependencies
use Fcntl ':flock';  # for locking on system level
use FindBin;         # to locate config.pl next to this script, regardless of cwd
use IO::Handle;      # for $lockFile->autoflush, so lock-file diagnostics are visible immediately

my $build;
if ($ARGV[0]) {
    $build = $ARGV[0];
}
else {
    die("publish.pl must be called with an argument that specifies the version to upload (e.g. v6.0.14). \n");
}

# Reject anything that isn't a plausible build/tag identifier before it is used
# to build local shell commands and the remote SSH command string.
if ($build !~ /^[A-Za-z0-9_.-]+$/) {
    die("Invalid build identifier '$build': only letters, digits, dots, dashes and underscores are allowed.\n");
}

print "Publishing\n";
print "------------------------------------------------------\n";
print "Build: $build.\n";

# Variables imported from config.pl
our ($SourceDir, $BuildsDir, $UploadDir, $ServerDir, $StableMajorVersion, $DefaultMajorVersion, $KeyForStableVersion, $KeyForDefaultVersion, $KeyForDeveloperVersion, $KeyForExperimentalVersion, $ChannelForDocumentation, $username, $hostname, @HTMLfiles, $KeyWebserver);

print "Loading config file...\n";
require "$FindBin::Bin/config.pl";  # load configuration file (always next to this script)

# Config paths may use a leading "~" for the home directory. system() is now
# called in list form (no shell involved, see below), so nothing expands "~"
# automatically; expand it here with glob(), which resolves "~" without
# shelling out. Falls back to the original string if glob() finds nothing.
if (defined $KeyWebserver) {
    my ($expandedKey) = glob($KeyWebserver);
    $KeyWebserver = $expandedKey if defined $expandedKey;
}

# Creating a lock file to prevent another instance of this script from running
# ---------------------------------------------------------------------------------
print "Creating lock file...\n";
my $numberOfAttempts = 3;  # number of attempts to obtain a lock
my $waitTime = 90;  # time to wait between attempts in seconds
my $lockFilename = "upload.lock";
my $lockPath = $SourceDir.'/'.$lockFilename;
open (my $lockFile,'>'.$lockPath) or die "Cannot open lock file '$lockPath': $!\n";
$lockFile->autoflush(1);

my $lockObtained = 0;
for my $attempt (1..$numberOfAttempts) {
    $lockObtained = flock $lockFile, LOCK_EX|LOCK_NB;
    last if $lockObtained;
    print("  Another instance of publish.pl appears to be running, trying again in $waitTime seconds.\n");
    sleep($waitTime) if $attempt < $numberOfAttempts;
}

if ($lockObtained) {
    # Record who holds the lock and since when, for diagnostics if a future run has to wait on it.
    print $lockFile "Lock file created by pid $$ at ".scalar(localtime)."\n";
    print "  Lock file created. \n";
}
else {
    print STDERR "  Cannot obtain lock after $numberOfAttempts attempts (another instance appears to be running), exiting.\n";
    close $lockFile;
    exit 1;
}

# Share one SSH connection across every ssh/scp call below instead of opening
# a new TCP+SSH handshake each time. Without this, the webserver's sshd can
# refuse or drop connections ("kex_exchange_identification: Connection closed
# by remote host") when too many separate connections arrive in the short
# burst this script produces (one scp per HTML file, plus the final upload
# and doc-extraction ssh call).
# ---------------------------------------------------------------------------------
my $webServerLogin = $username."@".$hostname;
my $sshControlPath = "/tmp/.easyspin-publish-ssh-$$";
my @sshMuxOpts = ('-o', 'ControlMaster=auto', '-o', "ControlPath=$sshControlPath", '-o', 'ControlPersist=10m');

# ---------------------------------------------------------------------------------
# Main publishing sequence, wrapped in eval so that any failure (die) still
# falls through to the cleanup below instead of leaving the lock file and the
# upload directory behind.
# ---------------------------------------------------------------------------------
my $success = eval {

    # Add key to webserver to keychain
    # ---------------------------------------------------------------------------------
    print "Adding SSH key to keychain...\n";
    my $sshadd_status = system('ssh-add', $KeyWebserver);  # private key to log into webserver
    die "Adding SSH key failed: exit code  $?\n" if $sshadd_status!=0;

    # Set up environment
    # ---------------------------------------------------------------------------------
    if (-e $UploadDir) {
        my $rm_status = system('rm', '-r', $UploadDir);
        die "Could not remove upload directory: exit code  $?\n" if $rm_status!=0;
    }

    my $mkdir_status = system('mkdir', $UploadDir);
    die "Could not create upload directory: exit code  $?\n" if $mkdir_status!=0;

    # Determine the release channel
    # ---------------------------------------------------------------------------------
    print "Determining release channel from build ID.\n";
    my $releaseChannel;
    # Anchored and dot-escaped so a version-looking substring embedded in an
    # unrelated tag (e.g. "nightly-2024.01.15") isn't mistaken for a semantic
    # version; an optional leading "v" is allowed since git tags commonly use it.
    my $matchPattern = '^v?(\d+)\.(\d+)\.(\d+)-?([a-z]+)?[-.]?(\d+)?$';

    my @buildID = ($build =~ m/$matchPattern/);

    if ($buildID[0]) {
        if ($buildID[3]) {
            $releaseChannel = $KeyForDeveloperVersion;
        }
        elsif ($buildID[0] eq $StableMajorVersion) {
            $releaseChannel = $KeyForStableVersion;
        }
        elsif ($buildID[0] eq $DefaultMajorVersion) {
            $releaseChannel = $KeyForDefaultVersion;
        }
    }
    else {
        # if tag does not follow semantic versioning, e.g. easyspin-evolve.zip
        $releaseChannel = $KeyForExperimentalVersion;
    }
    print "  Release channel: $releaseChannel.\n";

    # Git tags conventionally include a leading "v" (e.g. "v6.0.14"), but the
    # build artifact on disk, the HTML version strings, and the remote
    # directory names all use the bare version number, so strip it here.
    if ($buildID[0]) {
        $build =~ s/^v//;
    }

    # Look for zip file with the provided tag
    # ---------------------------------------------------------------------------------
    my $zipFileName = 'easyspin-'.$build.'.zip';

    print("Copying $zipFileName to upload directory. \n");
    die("$zipFileName does not exist in $BuildsDir \n") unless (-e "$BuildsDir$zipFileName");
    my $cp_status = system('cp', $BuildsDir.$zipFileName, $UploadDir.$zipFileName);
    die "cp command failed: exit code  $?\n" if $cp_status!=0;

    # Regexp to find versions and links to zip files in the html files
    # ---------------------------------------------------------------------------------
    # Escaped once for use inside the regex patterns below (the replacement-side
    # strings don't need this, since s/// treats the replacement as a plain
    # string, not a pattern).
    my $quotedChannel = quotemeta($releaseChannel);

    my $findLinkTozipFile = '<!--'.$quotedChannel.'zip-->';
    my $matchzipFile = "easyspin-(.*?).zip";
    my $replacezipFile = "easyspin-$build.zip";
    my $replaceLinkTozipFile = '<!--'.$releaseChannel.'--><a href="easyspin-'.$build.'.zip"><!--zip-->';

    my $matchOldVersion = '<!--'.$quotedChannel.'-->(.*?)<!--version-->';
    my $replaceOldVersion = '<!--'.$releaseChannel.'-->'.$build.'<!--version-->';

    my $matchInVersionsFile = "$quotedChannel:.*";
    my $replaceInVersionsFile = "$releaseChannel:$build";

    # Get html files from webserver and update them with the new version tags and zipfile names
    # ---------------------------------------------------------------------------------
    print "Downloading and updating HTML files...\n";
    foreach (@HTMLfiles) {
        my $currentFile = $_;
        print("Getting $currentFile...\n");

        # Download the current html file
        my $scp_dl_status = system('scp', @sshMuxOpts, $webServerLogin.':'.$ServerDir.$currentFile, $UploadDir.$currentFile.'.bak');
        die "scp command failed for $currentFile: exit code  $?\n" if $scp_dl_status!=0;

        # Scan through the html file and replace strings
        print("Updating $currentFile...\n");
        open(my $inputHTML,'<'.$UploadDir.$currentFile.'.bak') or die("Cannot open $currentFile.bak!");
        open(my $outputHTML,'>'.$UploadDir.$currentFile) or die("Cannot open $currentFile!");
        while (<$inputHTML>) {

            if ($_ =~ m/$findLinkTozipFile/) {
                $_ =~ s/$matchzipFile/$replacezipFile/g;
            }
            $_ =~ s/$matchOldVersion/$replaceOldVersion/g;
            $_ =~ s/$matchInVersionsFile/$replaceInVersionsFile/g;
            print $outputHTML $_;
        }

        close($inputHTML) or die("Cannot close $inputHTML!");
        close($outputHTML) or die("Cannot close $outputHTML!");
    }

    print("Deleting backup versions of html files...\n");
    my @bakFiles = glob($UploadDir.'*.bak');
    if (@bakFiles) {
        system('rm', @bakFiles) == 0
            or warn "Could not remove backup HTML files: exit code $?\n";
    }

    # Upload entire upload directory to webserver
    # ---------------------------------------------------------------------------------
    print("Uploading all new files to webserver via SCP...\n");

    my @uploadFiles = glob($UploadDir.'*');
    my $scp_status = @uploadFiles ? system('scp', @sshMuxOpts, @uploadFiles, $webServerLogin.':'.$ServerDir) : 0;
    die "SCP command failed: exit code  $?\n" if $scp_status!=0;

    # Unzip the build on the server and extract documentation
    # ---------------------------------------------------------------------------------
    # only happens for the release channel that is specified with $ChannelForDocumentation in the config file

    if ($releaseChannel eq $ChannelForDocumentation) {

        # Compose server-side command
        my $changeDir = "cd ".$ServerDir;
        my $rmFolders = "rm -rf ./documentation ./examples";
        my $unzipDoc = "unzip -qq ".$zipFileName." 'easyspin-".$build."/documentation/*' -d ./tmp/";
        my $unzipExamples = "unzip -qq ".$zipFileName." 'easyspin-$build/examples/*' -d ./tmp/";
        my $moveFiles = "cp -r ./tmp/easyspin-".$build."/* ./";
        my $rmTempDir = "rm -r ./tmp";
        my $issueCmd = join('; ', $changeDir, $rmFolders, $unzipExamples, $unzipDoc, $moveFiles, $rmTempDir);

        print("Unzipping new stable version and updating documentation and examples...\n");

        # Execute command via SSH
        my $ssh_status = system('ssh', @sshMuxOpts, '-o','IdentitiesOnly=yes','-i',$KeyWebserver,$webServerLogin,$issueCmd);
        die "SSH command failed: exit code  $?\n" if $ssh_status!=0;
    }

    1; # signal success to eval
};
my $error = $@;  # immediately capture error

if (!$success) {
    print STDERR "ERROR: $error";
}

# Remove the upload directory whether or not publishing succeeded
# ---------------------------------------------------------------------------------
if (-e $UploadDir) {
    print "Removing upload directory...\n";
    system('rm', '-r', $UploadDir) == 0
        or warn "Could not remove upload directory '$UploadDir': exit code $?\n";
}

# Close the shared SSH connection (if one was ever opened), so no background
# master process or control socket is left behind.
# ---------------------------------------------------------------------------------
if (-S $sshControlPath) {
    system('ssh', '-o', "ControlPath=$sshControlPath", '-O', 'exit', $webServerLogin);
}

# Remove lock file and exit
# ---------------------------------------------------------------------------------
print "Removing lock file...\n";
close $lockFile;
unlink($lockPath) or warn "Could not remove lock file: $!\n";

# ---------------------------------------------------------------------------------
if ($success) {
    print "Finished.\n";
    exit 0;
}
else {
    print "Finished with errors.\n";
    exit 1;
}
