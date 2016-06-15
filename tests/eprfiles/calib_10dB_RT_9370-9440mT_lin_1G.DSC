#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX,CPLX
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D,D
IIFMT	D,D
*
*	Data Ranges and Resolutions:
*
XPTS	2801
XMIN	93700.000000
XWID	700.000000
*
*	Documentational Text:
*
TITL	'calib_10dB_RT_9370-9440mT_lin_1G'
IRNAM	'Intensity',''
IINAM	'',''
XNAM	'Field'
IRUNI	'',''
IIUNI	'',''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    07/15/13
TIME    10:21:08
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    9.405
A1RS    2801
A1SW    0.07
MWFQ    2.634649e+11
MWPW    0.0015
AVGS    1
SPTP    0.04096
RCAG    40
RCHM    1
B0MA    0.0001
B0MF    100000
RCPH    225.0
RCOF    0.0
RCTC    0
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     fieldCtrl, 1.0

CenterField        94050.000 G
CryoFreeMagnet     Default
CurrWaitCurTol     0.030 V
CurrWaitVltTol     0.030 V
Delay              0.0 s
FieldWait          At rest (1st)
HighFieldVal       94400.000 G
HystPs5120         True
Linearize          True
LowFieldVal        93700.000 G
Resolution         2801
SweepDirection     Up
SweepProfile       Fast flyback
SweepWidth         700.000 G

.DVC     freqCounter, 1.0

FrequencyMon       263.464865 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              1.500 mW
PowerAtten         10.0 dB

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       0

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           0
ConvTime           40.96 ms
DModDetectSCT      First
DualDetect         0/90 Degrees
EliDelay           1.0 us
Enable1stHarm      True
Enable1stHarm90    True
Enable2ndHarm      False
Enable2ndHarm90    False
Enable3rdHarm      False
Enable3rdHarm90    False
Enable4thHarm      False
Enable4thHarm90    False
Enable5thHarm      False
Enable5thHarm90    False
EnableDisp         True
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               40 dB
GainB              39 dB
Harmonic           1
HighPass           True
InputPlugA         AC3
InputPlugB         AC3
Integrator         False
IsCalibExp         False
ModAmp             1.000 G
ModFreq            100.00 kHz
ModPhase           225.0
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadMode           True
ResCenter          0 mm
ResLength          0 mm
Resonator          1
SctNorm            False
SctRevision        Allegro
SetAllOrd          False
SetOrdDef          False
SpuExtension       True
SpuRevision        MultiHarmonic
SweepTime          114.72896 s
TimeConst          0
TimeExp            False
TuneCaps           247

.DVC     trigger, 1.0


*
************************************************************
