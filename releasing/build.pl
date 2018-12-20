use strict;
use warnings;

# other dependencies
use Fcntl ':flock'; # for locking on system level
use Net::SSH::Perl; # to use ssh

# directories
my $parentdir = './..'; # topdirectory 
my $builddir = $parentdir.'/easyspin-builds/'; # directory that builds are created in
my $uploaddir =  $parentdir.'/upload/'; # directory that will be used to upload files to easyspin.org
my $serverdir = '~/public_html/easyspin/test/'; # parentfolder on the server, used to copy files from $uploaddir into


# esbuild.m location
my $esbuild = './releasing/esbuild.m';

my $username = "easyspin";
my $hostname = "easyspin.org";

my $WebServerLogin = $username."@".$hostname;

my @htmlfiles = ("index.html","download.html","version.html");

# options for file locking - to ensure only one instance is running:
my $numberallowedattempts = 3; # number of attempts to obtain a lock
my $waittime = 90; # time to wait between attempts in seconds

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

# ---------------------------------------------------------------------------------
# Need to figure out release channel here
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

my $releasechannel = 'development';
my $branch = 'default';

# ---------------------------------------------------------------------------------
# Ensure ssh client is running and add keys to keychain
system('eval "$(ssh-agent -s)"');
system("ssh-add ~/.ssh/hostmonster_rsa"); # private key to log into hostmonster.com
system("ssh-add ~/.ssh/no_phrase_rsa"); # private key to log into bitbucket, needs to be adapted specific user

# ---------------------------------------------------------------------------------
# This section needs to be reworked
# ---------------------------------------------------------------------------------
# updating and pulling from bitbucket.com  
system('hg pull -u');

print "Switching and updating to $branch\n";
# system('hg update '.$branch.);

my $tag = `hg log -r "." --template "{latesttag}"`;
 
print "Current tag is: $tag \n";

# ---------------------------------------------------------------------------------
# Build .svg versions of formulae in the html documentation and replace the tags in html files
print "Converting LaTeX markup in documentation to svg and modifying html files. \n";
# system('perl ./releasing/mkformulas.pl');

# ---------------------------------------------------------------------------------
# Update ReleaseID (and betaVersion) in esbuild.m 
my $findReleaseID = "ReleaseID(.*?); \% major.minor.patch"; # pattern to find ReleaseID 
my $thisReleaseID = "ReleaseID = \'".$tag."\'; \% major.minor.patch"; # Update ReleaseID

# Extract semantic versioning from $tag
my ($major, $minor, $patch, $dev) = ($tag =~ m/(\d+).(\d+).(\d+)(.*)/);

if ($dev) {
   print "Non-numeric characters found in version info, making a developer build. \n";
}

my $esbuildNew = $esbuild.'new';

print("Updating ReleaseID in esbuild.m\n");
open(my $input,'<'.$esbuild) or die("Cannot open $esbuild!");
open(my $output,'>'.$esbuildNew) or die("Cannot open $esbuildNew!");
while (<$input>) {
    $_ =~ s/$findReleaseID/$thisReleaseID/g;
    # do we need this? --------------------------------------
    if ($dev){
        $_ =~ s/betaVersion(.*?);/betaVersion = false;/g;
        }
    else {
        $_ =~ s/betaVersion(.*?);/betaVersion = false;/g;
    }
    # -------------------------------------------------------
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
opendir(my $tempdir,$builddir) or die("Can't open $builddir: $!\n");

my $EntriesInBuilddir = () = readdir($tempdir);
my $CurrentEntriesInBuilddir = $EntriesInBuilddir;

closedir($tempdir);

print("Triggering Matlab build \n");
system('matlab.exe '.$MatlabOptions." ".$MatlabTarget);

print("Waiting for Matlab to finish building EasySpin $tag. \n");

while ($CurrentEntriesInBuilddir == $EntriesInBuilddir) {
    opendir(my $tempdir,$builddir);
    $CurrentEntriesInBuilddir = () = readdir($tempdir);
    closedir($tempdir);
    sleep 5;
}
sleep 10;

# copy easyspin.zip file to upload directory
my $zipfilename = 'easyspin-'.$tag.'.zip';
print("Copying $zipfilename to upload directory. \n");
system('cp '.$builddir.$zipfilename.' '.$uploaddir.$zipfilename);

# ---------------------------------------------------------------------------------
# get html files from easyspin.org and update them with the new version tags and zipfile names

# regexp to find versions and links to zip files in the html files
my $oldzipfile = '<!--'.$releasechannel.'-->(.*?)<!--zip-->';
my $newzipfile = '<!--'.$releasechannel.'--><a href="easyspin-'.$tag.'.zip"><!--zip-->';

my $oldversion = '<!--'.$releasechannel.'-->(.*?)<!--version-->';
my $newversion = '<!--'.$releasechannel.'-->'.$tag.'<!--version-->';

# loop over all the html files
foreach (@htmlfiles)
{
    my $currentfile = $_;
    print("Getting $currentfile \n");
    # download the current zipfile from easyspin.org
    system('scp '.$WebServerLogin.':'.$serverdir.$currentfile.' '.$uploaddir.$currentfile.'.bak');
    
    print("Updating index.html \n");
    open(my $inputhtml,'<'.$uploaddir.$currentfile.'.bak') or die("Cannot open $currentfile.bak!");
    open(my $outputhtml,'>'.$uploaddir.$currentfile) or die("Cannot open $currentfile!");
    while (<$inputhtml>) {
        $_ =~ s/$oldzipfile/$newzipfile/g;
        $_ =~ s/$oldversion/$newversion/g;
        print $outputhtml $_;    
    }

    close($inputhtml) or die("Cannot close $inputhtml!");
    close($outputhtml) or die("Cannot close $outputhtml!");  
}

print("Deleting backup versions of html files \n");
system('rm '.$uploaddir.'*.bak');

# upload entire upload directory to easyspin org and then clean it
print("Uploading all new files to easyspin.org \n");
system('scp '.$uploaddir.'* '.$WebServerLogin.':'.$serverdir);

# clear upload directory
print("Clean upload directory \n");
system("rm ".$uploaddir."*");

# ---------------------------------------------------------------------------------
# SSH into the server, unzip the build and extract documentation
# only happens for the 'stable' release channel!
if ($releasechannel eq 'stable') {
    print("Logging into easyspin.org via SSH \n");
    my $ssh = Net::SSH::Perl->new($hostname);
    $ssh -> login("$username");

    my $changedirectory = "cd ".$serverdir." \n";
    my $removefolders = "rm -r ./documentation ./examples \n";
    my $unzipdocumentation = qq(unzip -qq $zipfilename 'easyspin-$tag/documentation/*' -d ./tmp/ \n);
    my $unzipexamples = qq(unzip -qq $zipfilename 'easyspin-$tag/examples/*' -d ./tmp/ \n);
    my $movefiles = qq(mv ./tmp/easyspin-$tag/* ./ \n);
    my $rmtemp = qq(rm -r ./tmp \n);

    print("Unzipping new stable version and updating documentation and examples \n");
    my $cmd = $changedirectory.$removefolders.$unzipexamples.$unzipdocumentation.$movefiles.$rmtemp;

    my ($stdout,$stderr,$exit) = $ssh->cmd("$cmd");
}
# ---------------------------------------------------------------------------------
# Clean up lockfile and exit
close $lockfile;
system('rm '.$parentdir.'/'.$lockfilename);

print "Lock released.\n";
print "All finished.\n";