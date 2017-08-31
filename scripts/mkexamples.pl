#!/usr/bin/perl

$ExamplesMain = '../examples';

chdir($ExamplesMain);

# Get list of example subdirectories
opendir(ED,$ExamplesMain);
@allcategories = sort grep {-d && !/^\.\.?\z/} readdir(ED);
closedir(ED);

foreach $subdir (@allcategories) {
  $Description{$subdir} = $subdir;
}

$Description{'solidstate'} = 'Solid-state cw EPR simulations';
$Description{'isotropic'} = 'Isotropic cw EPR simulations';
$Description{'fastmotion'} = 'Fast-motion cw EPR simulations';
$Description{'slowmotion'} = 'Slow-motion cw EPR simulations';
$Description{'endor'} = 'ENDOR simulations';
$Description{'pulse'} = 'Pulse EPR simulations';
$Description{'shaping'} = 'Pulse shaping etc.';
$Description{'analysis'} = 'Data analysis';
$Description{'varia'} = 'Other examples';
$Description{'fitting'} = 'Least-squares fitting';
$Description{'magnetometry'} = 'Magnetometry';

open(D,'>../documentation/examplesmain.html');

print D  <<QQQ;

<!DOCTYPE html>
<html>
<head>
   <meta charset="utf-8">
   <link rel="icon" href="img/eslogo196.png">
   <link rel="stylesheet" type="text/css" href="style.css">
   <link rel="stylesheet" href="highlight/matlab.css">
   <script src="highlight/highlight.pack.js"></script>
   <script>hljs.initHighlightingOnLoad();</script>
   <title>Examples</title>
</head>

<body>

<header>
<ul>
<li><img src="img/eslogo42.png">
<li class="header-title">EasySpin
<li><a href="index.html">Documentation</a>
<li><a href="references.html">Publications</a>
<li><a href="http://easyspin.org" target="_blank">Website</a>
<li><a href="http://easyspin.org/forum" target="_blank">Forum</a>
</ul>
</header>

<section>


<!-- ============================================================= -->
<div class="subtitle" style="margin-top:0ex;">Examples</div>

<p>
Here is a collection of examples for the use of EasySpin, divided into
several categories. They are an excellent starting point for writing
your own code.
</p>

QQQ

foreach $category (@allcategories) {

  chdir($category);
  opendir(DIR,'.');
  @files = readdir(DIR);
  closedir(DIR);
  
  @files = sort(grep(/\.m$/,@files));

  print D '<b>'.$Description{$category}.'</b>';
  print D "<table width=100%>\n";

  @col = ('ffffff', 'f3f3f3');
  $icol = 1;
  foreach $ex (@files) {
    $thiscol = $col[$icol];
    $name = $ex;
    open EX, $name;
    $firstline = <EX>;
    close EX;
    $descr = $firstline;
    $descr =~ s/%\s*//;
    #$descr =~ s/%(.*)%\s*//;
    #$category = $1;
    #$category =~ s/\s*//g;
    $descr = ucfirst($descr);
    print D qq(<tr bgcolor="$thiscol"><td width=150><a href="../examples/$category/$name">$name</a></td><td>$descr</td></tr>\n);
    if ($icol==1) { $icol=2; } else { $icol=1; }
  }
  print D "</table><p></p>\n";
  chdir('..');
}

print D <<QQQ;

<hr>
</section>

<footer></footer>

</body>
</html>

QQQ
