#include <math.h>
#include <mex.h>

bool Display = false;

__inline int isodd(int k) { return (k % 2); }
__inline int parity(int k) { return (isodd(k) ? -1 : +1); }
__inline int mini(int a,int b) { return ((a<b) ? a : b); }

/* Lanczos approximation to the logarithm of the Gamma function
 for integer arguments.
 (taken from "Numerical Recipes", error 2e-10) */
double facln(const int N)
{
const double values[201] =
{0, 0,
0.69314718055994530942, 1.7917594692280550008, 3.1780538303479456196,
4.7874917427820459942, 6.5792512120101009951, 8.5251613610654143002, 
10.604602902745250228, 12.801827480081469611, 15.104412573075515295, 
17.502307845873885839, 19.987214495661886150, 22.552163853123422886, 
25.191221182738681500, 27.899271383840891566, 30.671860106080672804, 
33.505073450136888884, 36.395445208033053576, 39.339884187199494036, 
42.335616460753485030, 45.380138898476908026, 48.471181351835223880, 
51.606675567764373570, 54.784729398112319190, 58.003605222980519939, 
61.261701761002001985, 64.557538627006331059, 67.889743137181534983, 
71.257038967168009010, 74.658236348830164385, 78.092223553315310631, 
81.557959456115037179, 85.054467017581517414, 88.580827542197678804, 
92.136175603687092483, 95.719694542143202485, 99.330612454787426929, 
102.96819861451381270, 106.63176026064345913, 110.32063971475739543, 
114.03421178146170323, 117.77188139974507154, 121.53308151543863396, 
125.31727114935689513, 129.12393363912721488, 132.95257503561630988, 
136.80272263732636847, 140.67392364823425940, 144.56574394634488601, 
148.47776695177303207, 152.40959258449735784, 156.36083630307878519, 
160.33112821663090703, 164.32011226319518141, 168.32744544842765233, 
172.35279713916280156, 176.39584840699735172, 180.45629141754377105, 
184.53382886144949050, 188.62817342367159119, 192.73904728784490244, 
196.86618167288999399, 201.00931639928152668, 205.16819948264119854, 
209.34258675253683565, 213.53224149456326119, 217.73693411395422725, 
221.95644181913033395, 226.19054832372759333, 230.43904356577695232, 
234.70172344281826774, 238.97838956183432305, 243.26884900298271418, 
247.57291409618688394, 251.89040220972319438, 256.22113555000952546, 
260.56494097186320931, 264.92164979855280104, 269.29109765101982254, 
273.67312428569370415, 278.06757344036614291, 282.47429268763039603, 
286.89313329542699395, 291.32395009427030757, 295.76660135076062402, 
300.22094864701413175, 304.68685676566871547, 309.16419358014692194, 
313.65282994987906178, 318.15263962020932685, 322.66349912672617689, 
327.18528770377521720, 331.71788719692847314, 336.26118197919847703, 
340.81505887079901787, 345.37940706226685411, 349.95411804077023693, 
354.53908551944080885, 359.13420536957539878, 363.73937555556349014, 
368.35449607240474959, 372.97946888568902068, 377.61419787391865645, 
382.25858877306002911, 386.91254912321755248, 391.57598821732961963, 
396.24881705179152580, 400.93094827891574549, 405.62229616114488919, 
410.32277652693730542, 415.03230672824963956, 419.75080559954473410, 
424.47819341825707467, 429.21439186665157013, 433.95932399501482019, 
438.71291418612118484, 443.47508812091894096, 448.24577274538460572, 
453.02489623849613510, 457.81238798127818110, 462.60817852687492219, 
467.41219957160817874, 472.22438392698059624, 477.04466549258563310, 
481.87297922988793423, 486.70926113683941223, 491.55344822329800350, 
496.40547848721762066, 501.26529089157929278, 506.13282534203487520, 
511.00802266523602674, 515.89082458782239760, 520.78117371604415136, 
525.67901351599506273, 530.58428829443349218, 535.49694318016954419, 
540.41692410599766910, 545.34417779115487380, 550.27865172428556555, 
555.22029414689486985, 560.16905403727303813, 565.12488109487429886, 
570.08772572513420614, 575.05753902471020676, 580.03427276713078116, 
585.01787938883911760, 590.00831197561785390, 595.00552424938196897, 
600.00947055532742811, 605.02010584942368386, 610.03738568623860819, 
615.06126620708488458, 620.09170412847732004, 625.12865673089094920, 
630.17208184781019582, 635.22193785505973286, 640.27818366040804092, 
645.34077869343500772, 650.40968289565523925, 655.48485671088906617, 
660.56626107587352917, 665.65385741110591324, 670.74760761191267558, 
675.84747403973687400, 680.95341951363745461, 686.06540730199399784, 
691.18340111441075295, 696.30736509381401187, 701.43726380873708535, 
706.57306224578734711, 711.71472580229000695, 716.86222027910346000, 
722.01551187360123894, 727.17456717281576797, 732.33935314673928203, 
737.50983714177743381, 742.68598687435126295, 747.86777042464334810, 
753.05515623048410309, 758.24811308137431347, 763.44661011264013922, 
768.65061679971693457, 773.86010295255835551, 779.07503871016734113, 
784.29539453524566594, 789.52114120895886719, 794.75224982581345382, 
799.98869178864340302, 805.23043880370304540, 810.47746287586353154, 
815.72973630391016142, 820.98723167593794297, 826.24992186484282852, 
831.51778002390615665, 836.79077958246990345, 842.06889424170042068, 
847.35209797043840919, 852.64036500113294442, 857.93366982585743682, 
863.23198719240547350};

const double c[6] = {76.18009172947146,-86.50532032941677,
  24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
  -0.5395239384953e-5};

double x, y, tmp, Series, lnGamma;
int j;

if (N<=200)
  return values[N];

x = N + 1;
tmp = x + 5.5;
tmp -= (x+0.5)*log(tmp);

y = N + 1;
Series = 1 + 1.90015e-10;
for (j=0;j<6;j++) Series +=c[j]/++y;

lnGamma = -tmp + log(2.5066282746310005024*Series/x);

return lnGamma;
}

/*=========================================================================*/
double jjj(int j1, int j2, int j3, int m1, int m2, int m3)
{

double w3j;


int k;
const int prty = parity(j1+j2+j3);
int phase = 1;

/* Guarantee j1>=j2>=j3 */
if (j1<j3) {
  k = j1; j1 = j3; j3 = k;
  k = m1; m1 = m3; m3 = k;
  phase *= prty;
}
if (j1<j2) {
  k = j1; j1 = j2; j2 = k;
  k = m1; m1 = m2; m2 = k;
  phase *= prty;
}
if (j2<j3) {
  k = j2; j2 = j3; j3 = k;
  k = m2; m2 = m3; m3 = k;
  phase *= prty;
}

/* Guarantee m3>=0 */
if (m3<0) {
  m1 = -m1; m2 = -m2; m3 = -m3;
  phase *= prty;
}


/* Check triangle condition */
if (j1-j2>j3) {
  w3j = 0;
  return w3j;
}

if (m3<0)
  mexPrintf("  (%d %d %d | %d %d %d)\n",j1,j2,j3,m1,m2,m3);


if (j3<=2) {
  /* j3<=2 (covers almost all cases) */
  /* Explicit formulae from Edmonds' book (Table 2, p.125-127) */
  
  double tmp;
    
  int jd = j1 - j2;
  int J = j2;
  int M = m1;

  int x = 2*J - 1 + jd;
  int y = J - M;
  int z = J + M;
  
  /* Prefactor (-1)^(J-M) */
  phase *= parity(J-M);

  switch (j3) {
  case 0:
  
    w3j = 1/sqrt(2.0*J+1);
    break;
  
  case 1:
  
    if (jd==0) {
      if (m3==0) 
        w3j = M/sqrt(J*(J+1)*(2.0*J+1));
      else
        w3j = sqrt(2.0*y*(z+1)/((x+3)*(x+2)*(x+1)));
    }
    else {
      tmp = (x+1)*(x+2)*(x+3);
      if (m3==0)
        w3j = -sqrt(2*(y+1)*(z+1)/tmp);
      else
        w3j = -sqrt(y*(y+1)/tmp);
    }
    break;
  
  case 2:
  
    tmp = x*(x+1)*(x+2)*(x+3)*(x+4);
    switch (m3) {
    case 0:
      switch (jd) {
      case 0:
        w3j = 2*(3*M*M-J*(J+1))/sqrt(tmp);
        break;
      case 1:
        w3j = -2*M*sqrt(6*(z+1)*(y+1)/tmp); break;
      default:
        w3j = sqrt(6*(z+2)*(z+1)*(y+2)*(y+1)/tmp); break;
      }
      break;
    case 1:
      switch (jd) {
      case 0:
        w3j = (2*M+1)*sqrt(6*(z+1)*y/tmp); break;
      case 1:
        w3j = -(2*J+4*M+4)*sqrt(y*(y+1)/tmp); break;
      default:
        w3j = 2*sqrt((z+2)*y*(y+1)*(y+2)/tmp); break;
      }
      break;
    case 2:
      switch (jd) {
      case 0:
        w3j = sqrt(6*(y-1)*y*(z+1)*(z+2)/tmp); break;
      case 1:
        w3j = -2*sqrt((y-1)*y*(y+1)*(z+2)/tmp); break;
      default:
        w3j = sqrt((y-1)*y*(y+1)*(y+2)/tmp); break;
      }
      break;
    }
    break;
  }
}
else {
  /* General algorithm */

  double delta_ln = facln(j1+j2-j3) + facln(j1-j2+j3) + facln(-j1+j2+j3) - facln(j1+j2+j3+1);
  double fact_ln = facln(j1+m1) + facln(j1-m1) + facln(j2+m2) + facln(j2-m2) + facln(j3+m3) + facln(j3-m3);

  int t1 = j2 - m1 - j3;
  int t2 = j1 + m2 - j3;
  int t3 = j1 + j2 - j3;
  int t4 = j1 - m1;
  int t5 = j2 + m2;
  int tmin, tmax, t;
  double Sum, denom_ln;
  
  tmin = (t1>t2) ? t1 : t2; tmin = (tmin>0) ? tmin : 0;
  tmax = (t3<t4) ? t3 : t4; tmax = (tmax<t5) ? tmax : t5;

  /*t = max([t1,t2,0]) : min([t3,t4,t5]); */
  
  Sum = 0;
  for (t = tmin;t<=tmax;t++) {
    denom_ln = facln(t) + facln(t-t1) + facln(t-t2) + facln(t3-t) + facln(t4-t) + facln(t5-t);
    Sum += parity(t) * exp(delta_ln/2 + fact_ln/2 - denom_ln);
  }
  
  w3j = parity(j1-j2-m3) * Sum;
  /*value = (-1)^(j1-j2-m3) * sum((-1).^t.*exp(delta_ln/2 + fact_ln/2 - denom_ln)); */
  
}

return phase*w3j;

}

long maxElements = 1000000;
int maxRows = 50000;

const double sqrt12 = 0.70710678118655; /* sqrt(1.0/2.0); */
const double sqrt13 = 0.57735026918963; /* sqrt(1.0/3.0); */
const double sqrt23 = 0.81649658092773; /* sqrt(2.0/3.0); */

double *ridx, *cidx;
double *MatrixRe, *MatrixIm;

int nRows, nElements;

int Lemax, Lomax, Kmax, Mmax;
int jKmin, pSmin, pImax, deltaL, deltaK;

struct SystemStruct {
  double I;
  double EZI0, NZI0, HFI0;
  double *ReEZI2, *ImEZI2, *ReHFI2, *ImHFI2;
};

struct DiffusionStruct {
  double* xlk, *d2psi;
  double psi, Rxx, Ryy, Rzz;
  double Exchange;
  int maxL;
};

struct SystemStruct Sys;
struct DiffusionStruct Diff;

/*============================================================================ */
/*============================================================================ */
/*============================================================================ */
/*============================================================================ */
void makematrix()
{

const double I = Sys.I;
const double NZI0 = Sys.NZI0;
const double HFI0 = Sys.HFI0;
const double EZI0 = Sys.EZI0;
const double *ReEZI2 = Sys.ReEZI2;
const double *ImEZI2 = Sys.ImEZI2;
const double *ReHFI2 = Sys.ReHFI2;
const double *ImHFI2 = Sys.ImHFI2;

/* Diffusion parameters */
/*--------------------------------------------- */
const double psi = Diff.psi;
const double *d2psi = Diff.d2psi;
const double Rxx = Diff.Rxx;
const double Ryy = Diff.Ryy;
const double Rzz = Diff.Rzz;

const double ExchangeFreq = Diff.Exchange;
const double *xlk = Diff.xlk;

/* Lband = (lptmx>=2) ? lptmx*2 : 2; */
const int Lband = 8;
/* Kband = kptmx*2; */
const int Kband = 8;
/*========================================================== */

/* Loop variables */
int L1, jK1, K1, M1, pS1, qS1, pI1, qI1;
int L2, jK2, K2, M2, pS2, qS2, pI2, qI2;
int Ld, Ls, jKd, Kd, Ks, KK;

/* Min and max values for loop variables */
int K1max, M1max, qS1max, qI1max;
int L2max, jK2min, K2min, K2max, M2min, M2max;
int pS2min, qS2min, qS2max, pI2min, qI2min, qI2max;

int iRow = 0, iCol = 0, iElement = 0;

bool Potential = Diff.maxL>=0;
bool Exchange = (ExchangeFreq!=0);

double IsoDiffKdiag, IsoDiffKm2, IsoDiffKp2, PotDiff;
bool diagRC, Ld2, diagLK, diagI;
double N_L, N_K, NormFactor, R_EZI2, R_HFI2, a1, g1, a2, g2;
int parityLK2, L, pId, qId;
double Term1, Term2, Xd, Xs;
int pd;
bool includeRank0;
double LiouvilleElement, GammaElement, t;
double d2jjj;

const bool RhombicDiff = (Rxx!=Ryy);

if (Display) {
mexPrintf("   Potential: %d\n",Potential);
mexPrintf("   Exchange:  %d\n",Exchange);
mexPrintf("   Lemax Lomax deltaL:   %d  %d  %d\n",Lemax,Lomax,deltaL);
mexPrintf("   jKmin Kmax Mmax:   %d  %d  %d\n",jKmin,Kmax,Mmax);
mexPrintf("   pSmin  pImax:   %d  %d\n",pSmin, pImax);
}


/* All equation numbers refer to Meirovich et al, J.Chem.Phys. 77 (1982) */

for (L1=0;L1<=Lemax;L1+=deltaL) {
  if (isodd(L1) && (L1>Lomax)) continue;
  for (jK1=jKmin;jK1<=1;jK1+=2) {
    K1max = mini(Kmax,L1);
    for (K1=0;K1<=K1max;K1+=deltaK) {
      if ((K1==0)&&(parity(L1)!=jK1)) continue;

      /* Potential-independent part of diffusion operator */
      /*-------------------------------------------------------- */
      IsoDiffKdiag = (Rxx+Ryy)/2*(L1*(L1+1))+K1*K1*(Rzz-(Rxx+Ryy)/2);
      if (RhombicDiff) {
        KK = K1-2;
        IsoDiffKm2 = (Rxx-Ryy)/4*sqrt((L1-KK-1)*(L1-KK)*(L1+KK+1)*(L1+KK+2));
        KK = K1+2;
        IsoDiffKp2 = (Rxx-Ryy)/4*sqrt((L1+KK-1)*(L1+KK)*(L1-KK+1)*(L1-KK+2));
      }
      else {
        IsoDiffKp2 = IsoDiffKm2 = 0;
      }
/*
      else {
        IsoDiffKdiag = Rxx*L1*(L1+1)*pow(1+pmxl*L1*(L1+1),-pl) +
                        K1*K1*(Rzz*pow(1+ pmz*K1*K1,-pkzz) -
                               Rxx*pow(1+pmxk*K1*K1,-pkxy));
      }
*/
      /*-------------------------------------------------------- */
      
      
      M1max = mini(Mmax,L1);
      for (M1=-M1max;M1<=M1max;M1++) {
        
        for (pS1=pSmin;pS1<=1;pS1++) {
          qS1max = 1 - abs(pS1);
          for (qS1=-qS1max;qS1<=qS1max;qS1+=2) {
            for (pI1 = -pImax;pI1<=pImax;pI1++) {
              if ((psi==0) && ((pI1+pS1-1)!=M1)) continue;
              qI1max = ((int)(2*I)) - abs(pI1);
              for (qI1=-qI1max;qI1<=qI1max;qI1+=2) {
                /*mexPrintf("%3d %3d %3d %3d    %2d %2d %2d %2d\n",L1,jK1,K1,M1,pS1,qS1,pI1,qI1); */
                iCol = iRow;
                diagRC = true;
                L2max = mini(Lemax,L1+Lband);
                for (L2=L1;L2<=L2max;L2+=deltaL) {
                  if (isodd(L2)&&(L2>Lomax)) continue;
                  Ld = L1 - L2; Ls = L1 + L2;
                  Ld2 = abs(Ld)<=2;

                  /* N_L normalisation factor, see after Eq. (A11) */
                  N_L = sqrt((double)((2.0*L1+1.0)*(2.0*L2+1.0)));

                  jK2min  = (diagRC) ?  jK1 : jKmin;
                  for (jK2=jK2min;jK2<=1;jK2+=2) {
                    jKd = jK1 - jK2;
                    K2max = mini(Kmax,L2);
                    K2min = (diagRC) ? K1 : 0;
                    for (K2=K2min;K2<=K2max;K2+=deltaK) {
                      if ((K2==0)&&(parity(L2)!=jK2)) continue;
                      diagLK = (L1==L2) && (K1==K2);
                      Kd = K1 - K2;
                      Ks = K1 + K2;
                      parityLK2 = parity(L2+K2);
                      
                      /*---------------------------------------------- */
                      /* R(mu=EZI,HFI;l=2), see Eq. (A42) and (A44) */
                      /*---------------------------------------------- */
                      R_EZI2 = 0; R_HFI2 = 0;
                      if (Ld2) {
                        a1 = 0; g1 = 0;
                        if (abs(Kd)<=2) {
                          const double coeff = jjj(L1,2,L2,K1,-Kd,-K2);
                          /*mexPrintf("[wigner3j(%d,%d,%d,%d,%d,%d) %g]\n",L1,2,L2,K1,-Kd,-K2,coeff); */
                          if (jK1==jK2) {
                            a1 = coeff*ReHFI2[Kd+2];
                            g1 = coeff*ReEZI2[Kd+2];
                          }
                          else {
                            if (ImHFI2) a1 = coeff*ImHFI2[Kd+2]*jK1;
                            if (ImEZI2) g1 = coeff*ImEZI2[Kd+2]*jK1;
                          }
                        }
                        a2 = 0; g2 = 0;
                        if (abs(Ks)<=2) {
                          const double coeff = jjj(L1,2,L2,K1,-Ks,K2);
                          if (jK1==jK2) {
                            a2 = coeff*ReHFI2[Ks+2];
                            g2 = coeff*ReEZI2[Ks+2];
                          }
                          else {
                            if (ImHFI2) a2 = coeff*ImHFI2[Ks+2]*jK1;
                            if (ImEZI2) g2 = coeff*ImEZI2[Ks+2]*jK1;
                          }
                        }
                        R_EZI2 = g1 + jK2*parityLK2*g2;
                        R_HFI2 = a1 + jK2*parityLK2*a2;
                      }
                      /*---------------------------------------------- */

                      /* N_K(K_1,K_2)  normalization factor, Eq. (A43) */
                      N_K = 1.0;
                      if (K1==0) N_K /= sqrt(2.0);
                      if (K2==0) N_K /= sqrt(2.0);
                      
                      /* Normalisation prefactor in Eq.(A40) and Eq.(A41) */
                      NormFactor = N_L*N_K*parity(M1+K1);
                      
                      /*---------------------------------------------------------------------- */
                      /* Potential-dependent term of diffusion operator, Eq. (A40) */
                      /*---------------------------------------------------------------------- */
                      PotDiff = 0;
                      if (Potential) {
                        if ((abs(Ld)<=Lband)&&(parity(Ks)==1)&&
                            (jKd==0)&&(abs(Kd)<=Kband)&&(abs(Ks)<=Kband)) {
                          for (L=0; L<=Lband; L+=2) {
                            Term1 = Term2 = 0;
                            if (Kd>=-L) {
                              Xd = xlk[(Kd+L)*(Diff.maxL+1) + L];
                              if (Xd!=0)
                                Term1 = Xd * jjj(L1,L,L2,K1,-Kd,-K2);
                            }
                            if (Ks<=L) {
                              Xs = xlk[(Ks+L)*(Diff.maxL+1) + L];
                              if (Xs!=0)
                                Term2 = parityLK2*jK2* Xs *jjj(L1,L,L2,K1,-Ks,K2);
                            }
                            if (Term1 || Term2)
                              PotDiff += (Term1+Term2) * jjj(L1,L,L2,M1,0,-M1);
                          }
                          
                          PotDiff *= N_L*N_K*parity(M1+K1);
                          
                          if (Display)
                            if (abs(PotDiff)>1e-10)
                              mexPrintf(" pot (%d,%d;%d,%d): %e\n",L1,L2,K1,K2,PotDiff);
                        }
                      }
                      /*---------------------------------------------------------------------- */
                      
                      M2max = mini(Mmax,L2);
                      M2min = (diagRC) ? M1 : -M2max;
                      for (M2=M2min;M2<=M2max;M2++) {
                        int Md = M1 - M2;
                        
                        bool diagLKM = diagLK && (jKd==0) && (Md==0);

                        /* Pre-compute 3j symbol in Eq. (A41) for l = 2 */
                        const double Liou3j = (Ld2) ? jjj(L1,2,L2,M1,-Md,-M2) : 0;
                        
                        pS2min = (diagRC) ? pS1 : pSmin;
                        for (pS2=pS2min;pS2<=1;pS2++) {
                          int pSd = pS1 - pS2;
                          qS2max = 1 - abs(pS2);
                          qS2min = (diagRC) ? qS1 : -qS2max;
                          for (qS2=qS2min;qS2<=qS2max;qS2+=2) {
                            bool diagS = (pS1==pS2) && (qS1==qS2);
                            int qSd = qS1 - qS2;
                            pI2min = (diagRC) ? pI1 : -pImax;
                            for (pI2=pI2min;pI2<=pImax;pI2++) {
                              if ((psi==0) && ((pI2+pS2-1)!=M2)) continue;
                              pId = pI1 - pI2;
                              qI2max = ((int)(2*I)) - abs(pI2);
                              qI2min = (diagRC) ? qI1 : -qI2max;
                              for (qI2=qI2min;qI2<=qI2max;qI2+=2) {
                                qId = qI1 - qI2;
                                diagI = (pId==0) && (qId==0);
                                

                                /*-------------------------------------- */
                                /* Matrix element of Liouville operator */
                                /*-------------------------------------- */
                                /* The isotropic terms are not normalized by any factor. */
                                /* - The N_K factor is not needed because */
                                /*   the l=0 ISTO components have not been transformed by */
                                /*   the K-symmetrization. (following the formula, N_K is */
                                /*   cancelled by R_0) */
                                /* - The factor N_L*parityMK1 is canceled by the product of */
                                /*   wigner3j(L,M;0,0;L,-M) and wigner3j(L,K,0,0,L,-K), */
                                /*   which therefore do not need to be calculated. */
                                pd = pSd + pId;
                                includeRank0 = diagLKM && (pd==0);
                                
                                LiouvilleElement = 0;
                                
                                if (Ld2 && (abs(pSd)<=1) && (abs(pId)<=1) && 
                                   /*((psi!=0) || (pd==Md)) &&*/
                                   (abs(Md)<=2) &&
                                   (abs(pSd)==abs(qSd)) && (abs(pId)==abs(qId))) {
                                  
                                  d2jjj = d2psi[(pd+2)+(Md+2)*5]*Liou3j;
                                  
                                  /* Hyperfine interaction */
                                  /*----------------------------------------- */
                                  if ((I>0) && (pSd*pId==qSd*qId)) {
                                    /* Compute Clebsch-Gordan coeffs and S_A from Eq. (B7) */
                                    double C0, C2, S_A;
                                    if (pId==0) {
                                      if (pSd==0) {
                                        S_A = (pS1*qI1+pI1*qS1)/2.0;
                                        C0 = -sqrt13; /* (110|000) */
                                        C2 = +sqrt23; /* (112|000) */
                                      }
                                      else {
                                        S_A = -(pI1*pSd+qI1*qSd)/sqrt(8.0);
                                        C0 = 0; /* no rank 0 term */
                                        C2 = +sqrt12; /* (112|101), (112|-10-1) */
                                      }
                                    }
                                    else {
                                      int t = qI1*qId + pI1*pId;
                                      double KI = sqrt(I*(I+1.0)-t*(t-2.0)/4.0);
                                      if (pSd==0) {
                                        S_A = -(pS1*pId+qS1*qId)*KI/sqrt(8.0);
                                        C0 = 0; /* no rank 0 term */
                                        C2 = +sqrt12; /* (112|011), (112|0-1-1) */
                                      }
                                      else {
                                        S_A = pSd*qId*KI/2.0;
                                        C0 = +sqrt13; /* (110|1-10), (110|-110) */
                                        if (pd==0)
                                          C2 = +sqrt(1.0/6.0); /* (112|1-10), (112|-110) */
                                        else
                                          C2 = +1.0; /* (112|112), (112|-1-1-2) */
                                      }
                                    }
                                    /* Rank-2 term, Eq. (A41) */
                                    LiouvilleElement += NormFactor*d2jjj*R_HFI2 * C2*S_A;
                                    /* Rank-0 term */
                                    if (includeRank0)
                                      LiouvilleElement += HFI0*C0*S_A;
                                  }
                                  
                                  /* Electronic Zeeman interaction */
                                  /*----------------------------------------- */
                                  if (diagI) { /* i.e. pId==0 and qId==0 */
                                    /* Rank-2 term, Eq. (B7) */
                                    /* Compute Clebsch-Gordan coeffs and S_g Eq. (B8) */
                                    double C2, S_g;
                                    if (pSd==0) {
                                      C2 = +sqrt23; /* (112|000) */
                                      S_g = C2*pS1;
                                    }
                                    else {
                                      C2 = +sqrt12;  /* (112|-10-1), (112|101) */
                                      S_g = -C2*qSd/sqrt(2.0);
                                    }
                                    LiouvilleElement += NormFactor*d2jjj*R_EZI2*S_g;
                                    /* Rank-0 term */
                                    if (includeRank0) {
                                      const double C0 = -sqrt13; /* (110|000) */
                                      LiouvilleElement += C0*EZI0*pS1;
                                    }
                                  }
                                  
                                  /* Nuclear Zeeman interaction */
                                  /*----------------------------------------- */
                                  if (diagS && diagI && includeRank0) {
                                    const double C0 = -sqrt13; /* (110|000) */
                                    LiouvilleElement += pI1*C0*NZI0;
                                  }
                                  
                                }
                                /*------------------------------------------- */
                                
                                
                                /*------------------------------------------- */
                                /* Matrix element of diffusion operator       */
                                /*------------------------------------------- */
                                GammaElement = 0;
                                
                                /* Potential-independent terms */
                                if (diagS && diagI && (Ld==0) && (Md==0) && (jKd==0)) {
                                  if (Kd==0) GammaElement += IsoDiffKdiag;
                                  else if (Kd==+2) GammaElement += IsoDiffKm2/N_K;
                                  else if (Kd==-2) GammaElement += IsoDiffKp2/N_K;
                                }
                                
                                /* Potential-dependent terms */
                                if (Potential) {
                                  if (diagS && diagI && (Md==0) && (jKd==0))
                                    GammaElement += PotDiff;
                                }
                                
                                /* Exchange term */
                                if (Exchange) {
                                  if ((pSd==0) && (pId==0) && diagLKM) {
                                    t = 0;
                                    if ((qId==0) && (qSd==0)) t += 1.0;
                                    if ((qId==0) && (pS1==0)) t -= 0.5;
                                    if ((pI1==0) && (qSd==0)) t -= 1.0/(2.0*I+1);
                                    GammaElement += t*ExchangeFreq;
                                  }
                                }
                                /*------------------------------------------- */

                                
                                /*------------------------------------------- */
                                /* Store element values and indices           */
                                /*------------------------------------------- */
                                /*mexPrintf("%d/%d %d/%d %d/%d %d/%d\n",L1,L2,jK1,jK2,K1,K2,M1,M2); */
                                if ((GammaElement!=0) || (LiouvilleElement!=0)) {
                                  /*if (Display) mexPrintf("     saving element %5d: %f  %f\n",iElement+1,GammaElement,LiouvilleElement);*/
                                  MatrixRe[iElement] = GammaElement;
                                  MatrixIm[iElement] = -LiouvilleElement;
                                  ridx[iElement] = iRow;
                                  cidx[iElement] = iCol;
                                  iElement++;
                                  if (!diagRC) {
                                    MatrixRe[iElement] = GammaElement;
                                    MatrixIm[iElement] = -LiouvilleElement;
                                    ridx[iElement] = iCol;
                                    cidx[iElement] = iRow;
                                    iElement++;
                                  }
                                  if (iElement>=maxElements)
                                    mexErrMsgTxt("Number of non-zero elements too large. Increase Opt.Allocation(1).");
                                }
                                iCol++;
                                diagRC = false;
                                if (iCol>=maxRows)
                                  mexErrMsgTxt("Dimension too large. Increase Opt.Allocation(2).");
                                
                              }
                            } 
                          }
                        }
                      }
                    }
                  }
                } /* all column index loops */

                iRow++;
              }
            }
          }
        }
      }
    }
  }
} /* all row index loops */

nRows = iRow;
nElements = iElement;

}


/*======================================================================= */
/* Determine basis size */
int BasisSize()
/*======================================================================= */
{

/* Loop variables */
int L1, jK1, K1, M1, pS1, qS1, pI1, qI1;

/* Min and max values for loop variables */
int K1max, M1max, qS1max, qI1max;

int iRow = 0;
double I = Sys.I;

for (L1=0;L1<=Lemax;L1+=deltaL) {
  if (isodd(L1) && (L1>Lomax)) continue;
  for (jK1=jKmin;jK1<=1;jK1+=2) {
    K1max = mini(Kmax,L1);
    for (K1=0;K1<=K1max;K1+=deltaK) {
      if ((K1==0)&(parity(L1)!=jK1)) continue;
      M1max = mini(Mmax,L1);
      for (M1=-M1max;M1<=M1max;M1++) {
        for (pS1=pSmin;pS1<=1;pS1++) {
          qS1max = 1 - abs(pS1);
          for (qS1=-qS1max;qS1<=qS1max;qS1+=2) {
            for (pI1 = -pImax;pI1<=pImax;pI1++) {
              if ((Diff.psi==0)&&((pI1+pS1-1)!=M1)) continue;
              qI1max = (int)(2*I) - abs(pI1);
              for (qI1=-qI1max;qI1<=qI1max;qI1+=2) {
                /*mexPrintf("%3d %3d %3d %3d    %2d %2d %2d %2d\n",L1,jK1,K1,M1,pS1,qS1,pI1,qI1);*/
                iRow++;
              }
            }
          }
        }
		
      }
    }
  }
}

return iRow; /* number of rows */

}

/*============================================================================ */
/*============================================================================ */
/*============================================================================ */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  mxArray *T, *D;
  double *R;
  double *basisopts, *allocationOptions;
  int N, idxS, idxB, idxD;
  if (nrhs==1) {
    double *jm = mxGetPr(prhs[0]);
    /*mexPrintf("%0.10f\n",jjj(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5])); */
    return;
  }

  if (nrhs!=4) mexErrMsgTxt("4 input arguments expected!");
  if (nlhs!=5) mexErrMsgTxt("3 output arguments expected!");

  /*mexPrintf("Parsing system...\n"); */
  idxS = 0;
  Sys.I = mxGetScalar(mxGetField(prhs[idxS],0,"I"));
  Sys.EZI0 = mxGetScalar(mxGetField(prhs[idxS],0,"EZ0"));
  Sys.NZI0 = mxGetScalar(mxGetField(prhs[idxS],0,"NZ0"));
  Sys.HFI0 = mxGetScalar(mxGetField(prhs[idxS],0,"HF0"));

  T = mxGetField(prhs[idxS],0,"EZ2");
  Sys.ReEZI2 = mxGetPr(T);
  Sys.ImEZI2 = mxGetPi(T);
  T = mxGetField(prhs[idxS],0,"HF2");
  Sys.ReHFI2 = mxGetPr(T);
  Sys.ImHFI2 = mxGetPi(T);
  
  /*mexPrintf("Parsing basis structure...\n"); */
  idxB = 1;
  /*Lemax = (int)mxGetScalar(mxGetField(prhs[idxB],0,"evenLmax"));
  Lomax = (int)mxGetScalar(mxGetField(prhs[idxB],0,"oddLmax"));
  Kmax = (int)mxGetScalar(mxGetField(prhs[idxB],0,"Kmax"));
  Mmax = (int)mxGetScalar(mxGetField(prhs[idxB],0,"Mmax"));
  jKmin = (int)mxGetScalar(mxGetField(prhs[idxB],0,"jKmin"));
  pSmin = (int)mxGetScalar(mxGetField(prhs[idxB],0,"pSmin"));
  pImax = (int)mxGetScalar(mxGetField(prhs[idxB],0,"pImax"));
  deltaL = (int)mxGetScalar(mxGetField(prhs[idxB],0,"deltaL"));
  deltaK = (int)mxGetScalar(mxGetField(prhs[idxB],0,"deltaK"));
  */
  basisopts = mxGetPr(prhs[1]);
  
  Lemax = (int)basisopts[0];
  Lomax = (int)basisopts[1];
  Kmax = (int)basisopts[2];
  Mmax = (int)basisopts[3];
  jKmin = (int)basisopts[4];
  pSmin = (int)basisopts[5];
  pImax = (int)basisopts[6];
  deltaL = (int)basisopts[7];
  deltaK = (int)basisopts[8];

  /*
  mexPrintf("(%d %d %d %d) jKmin %d, pSmin %d, pImax %d, deltaL %d, deltaK %d\n",
    Lemax,Lomax,Kmax,Mmax,jKmin,pSmin,pImax,deltaL,deltaK);
  */
  
  /*mexPrintf("Parsing diffusion structure...\n"); */
  idxD = 2;
  Diff.Exchange = mxGetScalar(mxGetField(prhs[idxD],0,"Exchange"));
  Diff.psi = mxGetScalar(mxGetField(prhs[idxD],0,"psi"));
  Diff.xlk = mxGetPr(mxGetField(prhs[idxD],0,"xlk"));
  Diff.maxL = mxGetScalar(mxGetField(prhs[idxD],0,"maxL"));
  Diff.d2psi = mxGetPr(mxGetField(prhs[idxD],0,"d2psi"));
  
  R = mxGetPr(mxGetField(prhs[idxD],0,"Diff"));
  Diff.Rxx = R[0];
  Diff.Ryy = R[1];
  Diff.Rzz = R[2];

  allocationOptions = mxGetPr(prhs[3]);
  maxElements = (long)allocationOptions[0];
  if (Display) mexPrintf("  max elements: %d\n",maxElements);
  
  N = BasisSize();
  maxRows = N+1;

  if (Display) mexPrintf("  basis function count: %d\n",N);
  
  plhs[0] = mxCreateDoubleMatrix(maxElements,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(maxElements,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(maxElements,1,mxCOMPLEX);
  ridx = mxGetPr(plhs[0]);
  cidx = mxGetPr(plhs[1]);
  MatrixRe = mxGetPr(plhs[2]);
  MatrixIm = mxGetPi(plhs[2]);
  
  if (Display) mexPrintf("  starting matrix loops...\n");

  makematrix();

  if (Display) mexPrintf("  finishing matrix loops...\n");

  plhs[3] = mxCreateDoubleScalar(nRows);
  plhs[4] = mxCreateDoubleScalar(nElements);
  
  return;
}
