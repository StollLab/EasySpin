function computefreq() {
 var Field = document.cform.field.value;
 var gValue = document.cform.g.value;
 var Frequency = 0.01399624604872811*Field*gValue;
 document.cform.freq.value = Frequency;
}
function computefield() {
 var Frequency = document.cform.freq.value;
 var gValue = document.cform.g.value;
 var Field = Frequency/0.01399624604872811/gValue;
 document.cform.field.value = Field;
}
function computeg() {
 var Frequency = document.cform.freq.value;
 var Field = document.cform.field.value;
 var gValue = Frequency/0.01399624604872811/Field;
 document.cform.g.value = gValue;
}
function setfreq(frq) {
  document.cform.freq.value = frq;
}
function setg(g) {
  document.cform.g.value = g;
}

function dd_computedistance() {
var nu = document.ddform.dd_freq.value;
var r3 = 51.92052556862238/nu;
var r = Math.pow(r3,1/3);
document.ddform.dd_dist.value = r;
}

function dd_computefrequency() {
var r = document.ddform.dd_dist.value;
var nu = 51.92052556862238/r/r/r;
document.ddform.dd_freq.value = nu;
}

function dd_setdistance(r) {
  document.ddform.dd_dist.value = r;
}


function dh_computedistance() {
var nu = document.dhform.dh_freq.value;
var r3 = 0.078972741853244/nu;
var r = Math.pow(r3,1/3);
document.dhform.dh_dist.value = r;
}

function dh_computefrequency() {
var r = document.dhform.dh_dist.value;
var nu = 0.078972741853244/r/r/r;
document.dhform.dh_freq.value = nu;
}

function table(N) {
document.isotopelist.selectedelement.value = N;
}
