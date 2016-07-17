#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	MAN
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	I
*
*	Data Ranges and Resolutions:
*
XPTS	1024
XMIN	200
XWID	14000
*
*	Documentational Text:
*
TITL	'Q-Band oxygen ambient pressure'
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
DATE    10/02/98
TIME    13:47:58
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.72
A1SW    1.4
MWFQ    3.42773e+10
MWPW    0.11377
AVGS    1
RESO    qt9802
SPTP    0.16384
RCAG    50
RCHM    1
B0MA    0.0006
B0MF    100000
RCPH    4.0
RCOF    0.0
A1RS    1024
RCTC    0.08192
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*
.DVC     acqStart, 1.0
.DVC     fieldCtrl, 1.0
CenterField        7200.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
SweepDirection     Up
SweepWidth         14000.0 G
.DVC     fieldSweep, 1.0
.DVC     freqCounter, 1.0
FrequencyMon       34.277 GHz
.DVC     mwBridge, 1.0
AcqFineTuning      Never
Power              113.8 mW
PowerAtten         0.0 dB
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
ConvTime           163.84 ms
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
Gain               50 dB
Harmonic           1
HighPass           True
ModAmp             6.00 G
ModFreq            100.00 kHz
ModInput           Internal
ModOutput          Internal
ModPhase           4.0
Offset             0.0 %
QuadMode           False
QuadPhase          90.0
Resolution         1024
Resonator          1
SignalInput        Internal
SweepTime          167.77 s
TimeConst          81.92 ms
TuneCaps           46
*
************************************************************
*
#MHL	1.0 * MANIPULATION HISTORY LAYER by BRUKER
*
SOURCE_PRIM
'/usr/people/xuser/xeprFiles/Data/PEH/Qband_O2/oxi_1'
END_SOURCE_PRIM
SOURCE_SCND
  SOURCE
  '/usr/people/xuser/xeprFiles/Data/PEH/Qband_O2/oxi_1'
  END_SOURCE
  SELECT  'qualiRegions'
  'qualiRegions'
  END_SELECT
  PROCESS 'prLinRegr'
  REPORT
                  a = -1.710 +- 0.002
                  b = 620    +- 4
  reduced chi-square = 939.7
  END_REPORT
END_SOURCE_SCND
SELECT  'qualiRegions'
'qualiRegions'
END_SELECT
PROCESS 'prDiff'
PAR_VAL Gain(Sec.) =  1.000e+00
PAR_VAL x-Shift(Sec.) =  0.000e+00
PAR_VAL x-Stretch(Sec.) =  1.000e+00
MDATE   10/29/98 08:54:57
CHG_FMT 'integer'
MDATE   12/14/01 14:23:26
*
************************************************************
