function ok = test()

%                'cm^-1->eV'  'cm^-1->K'  'cm^-1->mT'  'cm^-1->MHz'       
%   'eV->cm^-1'               'eV->K'     'eV->mT'     'eV->MHz'
%   'K->cm^-1'   'K->eV'                  'K->mT'      'K->MHz'
%   'mT->cm^-1'  'mT->eV'     'mT->K'                  'mT->MHz'
%   'MHz->cm^-1' 'MHz->eV'    'MHz->K'    'MHz->mT'
ok = zeros(10,1);
v = unitconvert(unitconvert(1,'eV->cm^-1' ),'cm^-1->eV');
ok(1) = areequal(v,1,1e-10,'rel');
v = unitconvert(unitconvert(1,'K->cm^-1'),'cm^-1->K');
ok(2) = areequal(v,1,1e-10,'rel');
v = unitconvert(unitconvert(1,'mT->cm^-1'),'cm^-1->mT');
ok(3) = areequal(v,1,1e-10,'rel');
v = unitconvert(unitconvert(1,'MHz->cm^-1'),'cm^-1->MHz');
ok(4) = areequal(v,1,1e-10,'rel');

v = unitconvert(unitconvert(1,'K->eV'),'eV->K');
ok(5) = areequal(v,1,1e-10,'rel');
v = unitconvert(unitconvert(1,'mT->eV'),'eV->mT');
ok(6) = areequal(v,1,1e-10,'rel');
v = unitconvert(unitconvert(1,'MHz->eV'),'eV->MHz');
ok(7) = areequal(v,1,1e-10,'rel');

v = unitconvert(unitconvert(1,'mT->K'),'K->mT');
ok(8) = areequal(v,1,1e-10,'rel');
v = unitconvert(unitconvert(1,'MHz->K'),'K->MHz');
ok(9) = areequal(v,1,1e-10,'rel');

v = unitconvert(unitconvert(1,'MHz->mT'),'mT->MHz');
ok(10) = areequal(v,1,1e-10,'rel');


