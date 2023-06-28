function ok = test()

%                'cm^-1->eV'  'cm^-1->K'  'cm^-1->mT'  'cm^-1->MHz'       
%   'eV->cm^-1'               'eV->K'     'eV->mT'     'eV->MHz'
%   'K->cm^-1'   'K->eV'                  'K->mT'      'K->MHz'
%   'mT->cm^-1'  'mT->eV'     'mT->K'                  'mT->MHz'
%   'MHz->cm^-1' 'MHz->eV'    'MHz->K'    'MHz->mT'

v = unitconvert(rand,'cm^-1->eV');
v = unitconvert(rand,'cm^-1->K');
v = unitconvert(rand,'cm^-1->mT');
v = unitconvert(rand,'cm^-1->MHz' );

v = unitconvert(rand,'eV->cm^-1' );
v = unitconvert(rand,'eV->K');
v = unitconvert(rand,'eV->mT');
v = unitconvert(rand,'eV->MHz');

v = unitconvert(rand,'K->cm^-1');
v = unitconvert(rand,'K->eV');
v = unitconvert(rand,'K->mT');
v = unitconvert(rand,'K->MHz');

v = unitconvert(rand,'MHz->cm^-1');
v = unitconvert(rand,'MHz->eV');
v = unitconvert(rand,'MHz->K');
v = unitconvert(rand,'MHz->mT');

ok = true;
