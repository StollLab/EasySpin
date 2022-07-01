function [ok,data] = test(opt,olddata)

% regression test for pepper with photoselection
% results of script have been verified against an independent program
% with high fidelity

Sys = struct( 'S',1,'g',2,'lw',1.5,'D',[800,-120] ) ;
Sys.tdm = [ pi/2 0 ] ;
Sys.Pop = [ 1 0 0 ] ;

Exp = struct( 'mwFreq',9.75,'Range',[ 300 400 ],'Harmonic',0 ) ;
Exp.lightBeam = 'parallel' ;
Exp.lightScatter = 0.2 ;

[ x,spc_1 ] = pepper( Sys,Exp ) ;

Exp.lightBeam = 'perpendicular' ;

[ ~,spc_2 ] = pepper( Sys,Exp ) ;

if opt.Display
  if isempty(olddata)
    subplot( 1,2,1 )
    plot( x,spc_1 )
    title( 'E_{light} || B_0')
    xlabel( 'magnetic field / mT' )
    ylabel( 'intensity / a.u.' )
    subplot( 1,2,2 )
    plot( x,spc_2 )
    title( 'E_{light} \perp B_0')
    xlabel( 'magnetic field / mT' )
    ylabel( 'intensity / a.u.' )
  else
    subplot( 4,2,[ 1 3 5 ] )
    plot( x,spc_1,'k',x,olddata.spc(1,:),'r' )
    legend( 'old','new' )
    title( 'E_{light} || B_0')
    ylabel( 'intensity / a.u.' )
    subplot( 4,2,7 )
    plot( x,spc_1 - olddata.spc(1,:) )
    xlabel( 'magnetic field / mT' )
    subplot( 4,2,[ 2 4 6 ] )
    plot( x,spc_2,'k',x,olddata.spc(2,:),'r' )
    legend( 'old','new' )
    title( 'E_{light} \perp B_0')
    ylabel( 'intensity / a.u.' )
    subplot( 4,2,8 )
    plot( x,spc_2 - olddata.spc(2,:) )
    xlabel( 'magnetic field / mT' )
  end
end

spc = [ spc_1; spc_2 ] ;
data.spc = spc ;

if ~isempty(olddata)
  ok = areequal( spc,olddata.spc,1e-4,'rel' ) ;
else
  ok = [] ;
end