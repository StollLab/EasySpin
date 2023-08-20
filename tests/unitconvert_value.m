function ok = test()

%                'cm^-1->eV'  'cm^-1->K'  'cm^-1->mT'  'cm^-1->MHz'       
%   'eV->cm^-1'               'eV->K'     'eV->mT'     'eV->MHz'
%   'K->cm^-1'   'K->eV'                  'K->mT'      'K->MHz'
%   'mT->cm^-1'  'mT->eV'     'mT->K'                  'mT->MHz'
%   'MHz->cm^-1' 'MHz->eV'    'MHz->K'    'MHz->mT'
%
%    https://physics.nist.gov/cuu/Constants/energy.html
ok = zeros(16,1);
v = unitconvert(1,'cm^-1->eV');
ok(1) = areequal(v,1.239841e-4,1e-6,'rel');
v = unitconvert(1,'cm^-1->K');
ok(2) = areequal(v,1.438776,1e-6,'rel');
v = unitconvert(1,'cm^-1->mT');
ok(3) = areequal(v,1.069734e3,1e-6,'rel');
v = unitconvert(1,'cm^-1->MHz' );
ok(4) = areequal(v,2.997924e4,1e-6,'rel');

v = unitconvert(1,'eV->cm^-1' );
ok(5) = areequal(v,8.065543e3,1e-6,'rel');
v = unitconvert(1,'eV->K');
ok(6) = areequal(v,1.160451e4,1e-6,'rel');
v = unitconvert(1,'eV->mT');
ok(7) = areequal(v,8.627987e6,1e-6,'rel');
v = unitconvert(1,'eV->MHz');
ok(8) = areequal(v,2.417989e8,1e-6,'rel');

v = unitconvert(1,'K->cm^-1');
ok(9) = areequal(v,69.503480e-2,1e-6,'rel');
v = unitconvert(1,'K->eV');
ok(10) = areequal(v,8.617333e-5,1e-6,'rel');
v = unitconvert(1,'K->mT');
ok(11) = areequal(v,7.435024e2,1e-6,'rel');
v = unitconvert(1,'K->MHz');
ok(12) = areequal(v,2.083661e4,1e-6,'rel');

v = unitconvert(1,'MHz->cm^-1');
ok(13) = areequal(v,3.335640e-5,1e-6,'rel');
v = unitconvert(1,'MHz->eV');
ok(14) = areequal(v,4.135667e-9,1e-6,'rel');
v = unitconvert(1,'MHz->K');
ok(15) = areequal(v,4.799243e-5,1e-6,'rel');
v = unitconvert(1,'MHz->mT');
ok(16) = areequal(v,3.568248e-2,1e-6,'rel');