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
XMIN	3375.000000
XWID	50.000000
*
*	Documentational Text:
*
TITL	'cw:strong-pitch'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    raba
DATE    11/06/98
TIME    09:21:01
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.34
A1SW    0.005
MWFQ    9.58667e+09
MWPW    0.0200195
AVGS    1
SPTP    0.04096
RCAG    47
RCHM    1
B0MA    0.0001
B0MF    100000
RCPH    0.0
RCOF    0.7
A1RS    1024
RCTC    0.01024
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     fieldCtrl, 1.0

CenterField        3400.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
SweepDirection     Up
SweepWidth         50.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.586670 GHz

.DVC     mwBridge, 1.0

AcqFineTuning      Never
Power              20.02 mW
PowerAtten         10.0 dB

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
Calibrated         True
ConvTime           40.96 ms
DModAFCTrap        True
DModAmp            1.00 G
DModCalibrated     True
DModDetectSCT      First
DModEliDelay       1.0
DModExtLockIn      False
DModExtTrigger     False
DModFieldMod       First
DModGain           60 dB
DModHighPass       True
DModModOutput      Internal
DModSignalInput    Internal
DModTimeConst      1.28 ms
DoubleModFreq      5.00 kHz
DoubleModPhase     0.0
DoubleMode         False
EliDelay           1.0
ExtLockIn          False
ExtTrigger         False
Gain               47 dB
Harmonic           1
HighPass           True
ModAmp             1.00 G
ModFreq            100.00 kHz
ModInput           Internal
ModOutput          Internal
ModPhase           0.0
Offset             0.7 %
QuadMode           False
QuadPhase          90.0
Resolution         1024
Resonator          1
SignalInput        Internal
SweepTime          41.94 s
TimeConst          10.24 ms
TuneCaps           32

*
************************************************************
