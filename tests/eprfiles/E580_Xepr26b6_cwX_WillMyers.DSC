#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1024
XMIN	3280.000000
XWID	200.000000
*
*	Documentational Text:
*
TITL	'003_100K40dB124nsLeveled'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    04/13/11
TIME    17:43:58
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.338
A1SW    0.02
MWFQ    9.487074e+09
MWPW    1.988e-05
AVGS    124
RESO    Standard_6K
SPTP    0.04096
RCAG    60
RCHM    1
B0MA    5e-05
B0MF    100000
RCPH    353.0
RCOF    0.0
A1RS    1024
RCTC    0.08192
STMP    50
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     fieldCtrl, 1.0

AllegroMode        False
CenterField        3380.00 G
Delay              0.0 s
FieldFlyback       On
FieldResol         2.56000
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     12140.000 G
SweepDirection     Up
SweepWidth         200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.487074 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              0.01988 mW
PowerAtten         40 dB

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         124
NbScansDone        124
NbScansToDo        4096
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       1

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        False
Calibrated         True
ConvTime           40.96 ms
DModAFCTrap        True
DModAmp            1.00 G
DModCalibrated     True
DModDetectSCT      First
DModEliDelay       1.0 us
DModExtLockIn      False
DModExtTrigger     False
DModFieldMod       First
DModGain           60 dB
DModHighPass       True
DModIntegrator     True
DModModOutput      Internal
DModSignalInput    Internal
DModTimeConst      1.28 ms
DoubleModFreq      5.00 kHz
DoubleModPhase     0.0
DoubleMode         False
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               60 dB
Harmonic           1
HighPass           True
Integrator         False
ModAmp             0.50 G
ModFreq            100.00 kHz
ModInput           Internal
ModOutput          Internal
ModPhase           353.0
Offset             0.0 %
QuadMode           False
QuadPhase          90.0
Resolution         1024
Resonator          1
SamplingTime       40.96 ms
SctNorm            True
SctRevision        DSC2
SignalInput        Internal
SweepTime          41.94 s
TimeConst          81.92 ms
TimeExp            False
TuneCaps           32

*
************************************************************
