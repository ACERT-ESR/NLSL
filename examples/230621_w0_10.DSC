#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D,D,D,D,D,D,D,D,D,D
*
*	Data Ranges and Resolutions:
*
XPTS	3413
XMIN	3430.200000
XWID	149.956050
*
*	Documentational Text:
*
TITL	'230621_w0_10'
IRNAM	'Intensity','','','','','','','','',''
XNAM	'Field'
IRUNI	'','','','','','','','','',''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    06/21/23
TIME    15:26:31
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.35052
A1SW    0.015
MWFQ    9.835726e+09
MWPW    0.002
AVGS    15
SPTP    0.004
RCAG    81
RCHM    1
B0MA    0.000223
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    3413
RCTC    0
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        3505.20 G
Delay              0.5 s
FieldFlyback       On
FieldWait          Wait stable
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     3430.200 G
SweepDirection     Up
SweepWidth         150.0 G
WidthTM            200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.835726 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Each Slice Scan
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        58.9
Power              2.000 mW
PowerAtten         20.0 dB
QValue             8700

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         15
NbScansDone        15
NbScansToDo        15
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       1

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           9.271e-09
ConvTime           4.00 ms
DModDetectSCT      First
DoubleModAcc       1
DoubleModFreq      5.000 kHz
DoubleMode         False
DualDetect         OFF
EliDelay           1.0 us
Enable1stHarm      True
Enable1stHarm90    True
Enable2ndHarm      True
Enable2ndHarm90    True
Enable3rdHarm      True
Enable3rdHarm90    True
Enable4thHarm      True
Enable4thHarm90    True
Enable5thHarm      True
Enable5thHarm90    True
EnableDisp         False
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               81 dB
GainB              81 dB
Harmonic           1
HighPass           True
InputPlugA         AC3
InputPlugB         AC3
Integrator         False
IsCalibExp         False
ModAmp             2.230 G
ModFreq            100.00 kHz
ModPhase           0.0
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0.99652,0,0,0.00737177,0,0,-0.00559614,0,0,-2.88221e-05,0,0,1.00404e-05,0,0,3.43695e-08,0,0,-5.0404e-09,0,0,-1.4783e-11,0,0,-1.29132e-12,0
QuadMode           False
ResCenter          62.5 mm
ResLength          40 mm
Resolution         3413
Resonator          1
SctNorm            True
SctRevision        Allegro
SetAllOrd          False
SetOrdDef          False
SpuExtension       True
SpuRevision        MultiHarmonic
SweepTime          13.65200 s
TimeConst          0
TimeExp            False
TuneCaps           51
dModSW             True

*
************************************************************
