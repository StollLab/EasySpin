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

$Description{'analysis'} = 'Data analysis';
$Description{'endor'} = 'ENDOR simulations';
$Description{'fitting'} = 'Least-squares fitting';
$Description{'liquids'} = 'Isotropic cw EPR simulations';
$Description{'magnetometry'} = 'Magnetometry';
$Description{'photoexcitation'} = 'Spin-polarized systems';
$Description{'pulse evolve'} = 'Pulse EPR simulations (evolve)';
$Description{'pulse saffron'} = 'Pulse EPR simulations (saffron)';
$Description{'pulse spidyan'} = 'Pulse EPR simulations (spidyan)';
$Description{'pulse shaping'} = 'Pulse shaping etc.';
$Description{'slowmotion'} = 'Slow-motion cw EPR simulations';
$Description{'solidstate'} = 'Solid-state cw EPR simulations';
$Description{'trajectories'} = 'EPR spectra from MD trajectories';
$Description{'varia'} = 'Other examples';

open(D,'>../docsrc/examplesmain.html');

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
<li><a href="http://easyspin.org/academy" target="_blank">Academy</a>
<li><a href="http://easyspin.org/forum" target="_blank">Forum</a>
</ul>
</header>

<section>


<!-- ============================================================= -->
<div class="subtitle" style="margin-top:0ex;">Examples overview</div>

<p>
Here is a collection of examples for the use of EasySpin, divided into
several categories. They are an excellent starting point for writing
your own code.
</p>

<p>
The examples are organized in the following groups:
</p>

<ul>
<li><a href="#analysis">Data analysis</a></li>
<li><a href="#endor">ENDOR simulations</a></li>
<li><a href="#exchange">Chemical exchange</a></li>
<li><a href="#fastmotion">Fast-motion cw EPR simulations</a></li>
<li><a href="#fitting">Least-squares fitting</a></li>
<li><a href="#isotropic">Isotropic cw EPR simulations</a></li>
<li><a href="#magnetometry">Magnetometry</a></li>
<li><a href="#pulse evolve">Pulse EPR (evolve)</a></li>
<li><a href="#pulse saffron">Pulse EPR (saffron)</a></li>
<li><a href="#pulse shaping">Pulse shaping</a></li>
<li><a href="#pulse spidyan">Pulse EPR (spidyan)</a></li>
<li><a href="#slowmotion">Slow-motion cw EPR simulations</a></li>
<li><a href="#trajectories">Slow-motion cw EPR simulations from molecular-dynamics trajectories</a></li>
<li><a href="#solidstate">Solid-state cw EPR simulations</a></li>
<li><a href="#photoexcitation">EPR simulations for photoexcited systems</a></li>
<li><a href="#varia">Other examples</a></li>
</ul>

<p>
</p>

<!-- ============================================================= -->
<div class="subtitle">Examples list</div>

QQQ

foreach $category (@allcategories) {

  chdir($category);
  opendir(DIR,'.');
  @files = readdir(DIR);
  closedir(DIR);
  
  @files = sort(grep(/\.m$/,@files));

  print D "<!-- ===================================================================== -->\n";
  print D '<a name="'.$category.'"><b>'.$Description{$category}."</b></a>\n\n";
  print D "<table width=100%>\n";

  @col = ('ffffff', 'f3f3f3');
  $icol = 1;
  foreach $ex (@files) {
    $thiscol = $col[$icol];
    $exfilename = $ex;
    open EX, $exfilename;
    $firstline = <EX>;
    close EX;
    $description = $firstline;
    $description =~ s/%\s*//;
    $description = ucfirst($description);
    $exname = substr $exfilename, 0, -2;  # remove .m
    print D qq(<tr bgcolor="$thiscol">\n<td width=200><a href="../examples/$category/$exfilename">$exname</a></td>\n<td>$description</td></tr>\n);
    if ($icol==1) { $icol=2; } else { $icol=1; }
  }
  print D "</table>\n<p></p>\n\n";
  chdir('..');
}

print D <<QQQ;

<hr>
</section>

<footer></footer>

</body>
</html>
QQQ
