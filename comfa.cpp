
#include "openeye.h"

#include <fstream>
#include <iostream>
#include <vector>
#include "oeplatform.h"
#include "oesystem.h"
#include "oechem.h"
#include "oedepict.h"

using namespace std;
using namespace OEPlatform;
using namespace OEChem;
using namespace OESystem;
using namespace OEDepict;
using namespace OEMath;


const char *InterfaceData =
    "!PARAMETER -in\n"
    "  !TYPE string\n"
    "  !DEFAULT in.sdf\n"
    "  !BRIEF 3D from TCFA to model\n"
    "!END\n"
    "!PARAMETER -pred\n"
    "  !TYPE string\n"
    "  !DEFAULT in.sdf\n"
    "  !BRIEF 3D from TC to predict\n"
    "!END\n"
    "!PARAMETER -biolbl\n"
    "  !TYPE string\n"
    "  !DEFAULT LOGBIO\n"
    "  !BRIEF tag for measured affinity\n"
    "!END\n"
        "!PARAMETER -rbatt\n"
        "  !TYPE double\n"
        "  !DEFAULT 0.85\n"
        "  !BRIEF rotor attenuation factor\n"
        "!END\n"
        "!PARAMETER -ncomp\n"
        "  !TYPE int\n"
        "  !DEFAULT -1\n"
        "  !BRIEF # PLS components (-1 for SAMPLS)\n"
        "!END\n"
        "!PARAMETER -mkmodel\n"
        "  !TYPE bool\n"
        "  !DEFAULT true\n"
        "  !BRIEF making new model\n"
        "!END\n"
        "!PARAMETER -predict\n"
        "  !TYPE bool\n"
        "  !DEFAULT true\n"
        "  !BRIEF predicting (with new or existing) model\n"
        "!END\n"
        "!PARAMETER -model\n"
        "  !TYPE string\n"
        "  !DEFAULT model.tcfa\n"
        "  !BRIEF model file name (new or exixting)\n"
        "!END\n"
;


#include <fstream>
#include <iostream>
#include <cmath>

double vsd(int n, double *a ) {
// ***********************************
    double t = 0.0;
    for (int i = 0; i < n; ++i) t += a[i];
    t /= ((double)n);
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += (a[i] - t) * (a[i] - t );
    return sqrt(s/((double)(n-1)));
}

void ccfit( int maxnc, int nMol, double *cc, double *beta,
	double *vcoef, double *work, double *hwork, double *ycmpnt, 
	double *bio, double *y1 ) {
// **********************************************
    double pac = 1.0e-08;
    double pacd = sqrt((double)(nMol-1))*pac;
    double pacn = vsd(nMol,bio) * pacd;

	double *wy = new double[nMol];
	for (int i = 0; i < nMol*maxnc; i++) vcoef[i] = 0.0;

    double a = 0.0;
    for (int i = 0; i < nMol; ++i) a += bio[i];
    *y1 = a/((double)nMol);
	double wt = (double) (nMol - 1);

/// cross validation LOO loop = omitting each structure individually
  for (int lxv = 0; lxv < nMol; lxv++){

	  for (int i = 0; i < nMol; ++i) wy[i] = 1.0;
	  wy[lxv] = 0.0;

	  double a2 = 0.0;
	  for (int i = 0; i < nMol; ++i) {
		  work[i*3] = bio[i] - *y1;
		  a2 += wy[i] * work[i*3] * work[i*3];
	  }

  for (int nc = 0; nc < maxnc; ++nc) {
	  for (int i = 0; i < nMol; ++i) {
		work[i * 3 + 2] = work[i * 3];
		vcoef[i * maxnc + nc] = work[i * 3 + 2];
		}
	  double yy = 0.0;
	  for (int i = 0; i < nMol; ++i) yy += work[i*3+2] * work[i*3+2] * wy[i];
	  for (int j = 0; j < nMol; ++j) {
		a = 0.0;
		for (int i = 0; i < nMol; ++i)
			a += cc[i*nMol + j] * work[i*3+2] * wy[i];
		work[j*3 + 1] = a;
	  }
// Center "s"=w2, i.e. orthogonalize to any constant-vector
	  a = 0.0;
	  for (int i = 0; i < nMol; ++i) a += work[i*3 + 1] * wy[i];
	  a /= wt;
	  for (int i = 0; i < nMol; ++i) work[i*3 + 1] -= a;
// Orthogonalize "s"=w2 to all previous t(jh) 
	  for (int ic = 0; ic < nc-1; ic++) {
		  a = 0.0;
		  for (int i = 0; i < nMol; ++i) a += work[i*3 + 1] * ycmpnt[ i * maxnc + ic ] * wy[i];
		  //(use hwork(jh),the stored value of ycmpnt**2 for prev ic)
		  double f = a / hwork[ic];
		  for (int i = 0; i < nMol; ++i) work[ i*3 + 1] -= f * ycmpnt[ i * maxnc + ic ];
		  for (int i = 0; i < nMol; ++i) vcoef[ i*maxnc + nc ] -= f * vcoef[ i*maxnc + ic] ;
	    }
	  // Set (t)=(s)
	double tt = 0.0;
	double ty = 0.0;
	for (int i = 0; i < nMol; ++i) {
	    double w = work[ i*3 + 1 ];
	    tt += w*w * wy[i];
	    ty += w * work[ i*3 ] * wy[i];
	}
	double r = (ty + pacn) / (tt+pacd);
	beta[nc] = r; 
	for (int j = 0; j < nMol; j++) 
		ycmpnt[j*maxnc + nc] = r * work[j*3 + 1];
	for (int i = 0; i < nMol; ++i) 
		vcoef[ i*maxnc + nc ] = r * vcoef[ i*maxnc + nc ];
	a2 = 0.0;
	for (int i = 0; i < nMol; ++i) {
		a2 += ycmpnt[i*maxnc + nc] * ycmpnt[i*maxnc + nc] * wy[i];
		work[i*3] -= ycmpnt[i*maxnc + nc];
	}
	hwork[ nc ] = a2; 
    }
  }
}


void getCovar( const int nMol, const int nfld, const double *fld, double *Covar ) {
// **********************************************
/*
ofstream f;
	f.open("SAMPLS.tc");
	f.precision(6);
	f << "21 rows in this table\n";
*/
	double *dists = new double[ nMol * nMol ];
    for (int i = 0; i<nMol; i++) {
	 dists[ i*nMol + i ] = 0.0;  //diagonal
	 for (int j = 0; j < i; j++) {
		double sum = 0.0;
		for (int f = 0; f < nfld; f++) {
			double d1 = fld[i * nfld + f];
			double d2 = fld[j * nfld + f];
			sum += (d1 - d2) * (d1 - d2);
		}
		dists[i * nMol + j] = dists[j * nMol + i] = sum;
//		 f << i+1 << ' ' << j+1 << ' ' << sqrt(sum) << '\n';
	 }
    }
//	f.close();
    double s = 0.0;
	for (int i = 1; i < nMol; i++) for (int j = 0; j < i; j++) s+= dists[i*nMol + j];
	double s2 = s;
	double s1 = s2/((double)nMol);
	for (int i = 0; i < nMol; i++) {
		double t = -s1;
		for (int j = 0; j < nMol; j++) t += dists[i * nMol + j];
		Covar[i * nMol + i] = t / ((double) nMol);
	}
	for (int i = 1; i < nMol; i++) for (int j = 0; j < i; j++) {
			double t = 0.5 * (Covar[i * nMol + i] + Covar[j * nMol + j] - dists[i * nMol + j]);
			Covar[i * nMol + j] = Covar[j * nMol + i] = t;
		}
}

int sampls(int maxnc, int nMol, double *bio, 
	int nfld, double *fld, 
	double *sdep, double *q2, double *res) {
// **************************************************

    double *Covar = new double[ nMol * nMol ];
    getCovar( nMol, nfld, fld, Covar);

    double *ycmpnt = new double[nMol*maxnc];
    double wy[nMol];
    double *yrwork = new double[nMol*maxnc];
    double *beta = new double[maxnc];
    double *zwork = new double[nMol*maxnc];
    double *work =  new double[3*nMol];
    double *hwork = new double[maxnc];
	double *vcoeff = new double[nMol*maxnc];
	for (int i = 0; i < nMol*maxnc; i++) vcoeff[i] = 0.0;
	double y1work = 0.0;

    double avgBio = 0.0;
    for (int nm = 0; nm < nMol; nm++) avgBio += bio[nm];
    avgBio /= (double) nMol;

    for (int lxv=0; lxv < nMol; ++lxv) {
		for (int i = 0; i < nMol; ++i) wy[i] = 1.0;
		wy[lxv] = 0.0;

		ccfit( maxnc, nMol, Covar, beta, zwork, work, hwork, yrwork, bio, &y1work );

		for (int ic=0 ; ic < maxnc; ++ic) {
	    	ycmpnt[ lxv*maxnc + ic ] = yrwork[ lxv*maxnc + ic ];
	    	if (ic == 0) ycmpnt[lxv*maxnc] += (y1work - avgBio);
			for (int i = 0; i < nMol; i++) if (i != lxv)
					vcoeff[ i*maxnc + ic] += zwork[i*maxnc+ic]/((double) nMol);
		}
    }

    double a = 0.0;
	for (int nm = 0; nm < nMol; nm++) a += bio[nm];
	double y1 = a / ((double) nMol);
	double y2 = 0.0;
	for (int nm = 0; nm < nMol; nm++) y2 += (bio[nm] - y1) * (bio[nm] - y1);

	double *pwork = new double[ 2*maxnc ];
	for (int nc = 0; nc < maxnc; nc++ ) pwork[ maxnc + nc ] = 0.0;
	for (int nm = 0; nm < nMol; nm++) {
		double dy = bio[nm] - y1;
		for (int nc = 0; nc < maxnc; ++nc) {
			dy -= ycmpnt[nm * maxnc + nc];
			pwork[nc] = dy;
			pwork[ maxnc + nc ] += dy * dy;
		}
	}
    int optnc = 0;
    double optsd = sqrt( y2 / ((double)(nMol - 1)));
	double varBio = optsd;
    for (int nc = 0; nc < maxnc; ++nc) {
		int nDegF = nMol - nc - 1;
		nDegF = nDegF > 0 ? nDegF : 1;
		double nowsd = sqrt(pwork[maxnc + nc]/((double)(nDegF)) );
		if (nowsd > optsd) break;
		if (nowsd < 0.3) break;
		optsd = nowsd;
		optnc = nc + 1;
    }
    *sdep = optsd;
    *q2 = 1.0 - optsd*optsd/varBio;
    
    for (int nm = 0; nm < nMol; nm++) {
		double biodiff = bio[nm] - avgBio;
		for (int nc = 0; nc < optnc; ++nc) {
	    	biodiff -= ycmpnt[nm*maxnc+nc];
	    	pwork[nc] = biodiff;
		}
		res[nm] = pwork[optnc-1];
    }
    return optnc;
}
/*
int main() {

    double bio[14] = {-2.15,-1.28,-1.19,-1.0,-.75,-.63,-.42,-.40,-.15,
	-.05,.02,.04,.35,.40};
    double fld[42] = {-.61,4.4,9.63, -.61,6.0,11.0, 0,3.39,9.67, 0,3.3,8.34,
	-.05,2.243,10.321, .42,3.197,8.836, -.17,3.95,8.8, -.24,4.233,9.417,
	-.4,5.849,11.216, -.2,5.119,12.085, .23,4.21,9.113, .23,4.06,8.828,
	.03,4.655,10.07, .54,4.27,8.847};
    
    double sdep;
    double q2;
    double *res = new double[14];
    std::cout << sampls(3, 14, bio, 3, fld, &sdep, &q2, res );
}*/
#include <fstream>
#include <iostream>
#include <cmath>

#define FORARR(p,dim1size,dim1,dim2) ( p[(dim1size)*(dim2)+(dim1)] )

void matrix( int *pnvar1, int *pnvar2, int *picomp, double *beta, 
	double *w1, double *w2, double *b, double *ro, double *z ) {
// ************************************************************

 int nvar1 = *pnvar1;
 int nvar2 = *pnvar2;
 int icomp = *picomp + 1;

 int temp= nvar1* nvar2;
  /* Initialize z and beta  */
 for (int i=0;i<temp;i++) beta[i]=0. ;

 for (int ii=0;ii< nvar1; ii++)
  { for (int i=0;i<nvar1;i++) z[i] = 0.;
    z[ii]= 1.0 ;
     /* start the summation over icomp */
    for (int ic=0;ic<icomp;ic++)
     {double ss = 0.0 ;
        /* Dummy latent variable */
       for (int i=0;i<nvar1;i++)  ss += z[i] * FORARR(w1,nvar1, i,ic) ;
       ss *= ro[ic];
       for (int k=0;k<nvar2;k++)
        FORARR(beta,nvar1, ii,k) += ss * FORARR(w2,nvar2, k,ic);
        /* Update z matrix */
       ss=0;
       for (int j=0;j<nvar1;j++) ss   += z[j] * FORARR(w1,nvar1, j,ic)  ;
       for (int j=0;j<nvar1;j++) z[j] -= ss * FORARR(b,nvar1, j,ic) * ro[ic] ;
     } /* ic loop */
   }   /* ii loop */
}

double r1unif( int *pcseed, int *ptseed, int ibyte[4], int *pfcn) {
// ******************************************************
 double r1unif;
 int i,k,icarry,i1,i2,j1,j2,it1,it2;
 static int first = true;
 static int cseed[6] ={0},
            tseed[32]={0},
            iscr[5],
            t[29] = {1,2,3,3,2,1,4,5,6,7,5,4,7,6,1,6,7,4,5,2,3,7,6,5,4,3,2,1,0},
            jcseed = 12345,
            jtseed = 1073   ;
 r1unif = 0.0 ;

 if ( first || *pfcn>0 )
  { if (*pfcn <= 0) { jcseed = abs(*pcseed); jtseed = abs(*ptseed); }
    first = false;
    cseed[0] = jcseed;
    for (i=0;i<5;i++)
     { cseed[i+1] = cseed[i]  /64;
       cseed[i]  -= cseed[i+1]*64;  }
    cseed[5] = cseed[5] % 4;
    if (jcseed && (cseed[0] % 2)==0 ) cseed[0] ++;
    tseed[0] = jtseed;
    for (i=0;i<11;i++)
     { tseed[i+1] = tseed[i]  /2;
       tseed[i]  -= tseed[i+1]*2; }
    for (i=11;i<32;i++) tseed[0] = 0;
    if (jtseed) tseed[0] = 1;
    if ( ! *pfcn ) return r1unif;
  }

 for (i= 0;i<17;i++) tseed[i] = abs ( tseed[i] - tseed[i+15] );
 for (i=17;i<32;i++) tseed[i] = abs ( tseed[i] - tseed[i-17] );
 cseed[5] = 13* cseed[5] + 55* cseed[4] + 16* cseed[3];
 cseed[4] = 13* cseed[4] + 55* cseed[3] + 16* cseed[2];
 cseed[3] = 13* cseed[3] + 55* cseed[2] + 16* cseed[1];
 cseed[2] = 13* cseed[2] + 55* cseed[1] + 16* cseed[0];
 cseed[1] = 13* cseed[1] + 55* cseed[0] ;
 cseed[0] = 13* cseed[0] ;

 k= -6;
 icarry =0;
 for (i=0;i<5;i++)
  { k += 6;
    cseed[i] += icarry;
    icarry    = cseed[i]/64;
    cseed[i] -= 64*icarry;
    i2 = cseed[i] / 8;
    i1 = cseed[i] - 8*i2;
    j1 = 4*tseed[k+2]+tseed[k+1]+tseed[k+1]+tseed[k  ];
    j2 = 4*tseed[k+5]+tseed[k+4]+tseed[k+4]+tseed[k+3];
    it1 = 28;
    if (i1>j1) it1 = (i1*i1-i1)/2+j1;
    if (i1<j1) it1 = (j1*j1-j1)/2+i1;
    it2 = 28;
    if (i2>j2) it2 = (i2*i2-i2)/2+j2;
    if (i2<j2) it2 = (j2*j2-j2)/2+i2;
    iscr[i] = 8*t[it2]+t[it1];
    r1unif = (r1unif + iscr[i])  / 64.0 ;
  } /* i loop */
 cseed[5] = (cseed[5]+icarry) % 4;
 j1 = tseed[30]+tseed[31]+tseed[31];
 it1 = abs(cseed[5]-j1);
 if (it1 == 1 && (cseed[5]+j1)==3 ) it1=3;
 r1unif = (r1unif + it1) / 4.0 ;
 if ( *pfcn ==1 ) return r1unif;
 ibyte[3] = iscr[0] + (iscr[1] % 4) * 64 ;
 ibyte[2] = iscr[1]/4+(iscr[2] % 16)* 16 ;
 ibyte[1] = iscr[2]/16+iscr[3]*4;
 ibyte[0] = iscr[4]+it1*64;
 return r1unif;
}

void ranums( double *x, int *pn) {
// *******************************************************
 int i;
 static int ibyte[4],
            icseed   = 0,
            itseed   = 0,
            ifcn     = 1 ;
 for (i=0;i< *pn; i++) x[i] = r1unif( &icseed, &itseed, ibyte, &ifcn);
}

void resid( int *ppat, int *pvar1, int *pvar2, double *x, double *y, 
	double *weyt, double *xbar, double *xscal, double *off, double *beta, 
	double *res, double *ypred, double *ss, int *perr ) {
// *****************************************************

 *perr=0;
  int nvar1 = *pvar1;
  int nvar2 = *pvar2;

 for (int j=0;j< nvar2;j++) {
    ss[j] = 0. ;
    int jj = j+nvar1;
    off[j] = xbar[jj];
    for (int k=0;k<nvar1;k++) {
        FORARR(beta,nvar1,k,j) *= (xscal[jj]/xscal[k]) ;
       off[j] -= FORARR(beta,nvar1,k,j) *xbar[k];
    }
  } /* j loop */
 for (int i=0;i< *ppat;i++) {
     if (weyt[i]>=0.0) continue;
     for (int j=0;j<nvar2;j++) {
            double s  = off[j];
            for (int k=0;k<nvar1;k++)
                s += FORARR(x,nvar1,k,i) * FORARR(beta,nvar1,k,j);
            FORARR(ypred,nvar2,j,i) = s;
            FORARR(res,nvar2+1,j,i) = FORARR(y    ,nvar2,j,i)
                               - FORARR(ypred,nvar2,j,i)   ;
            s = FORARR(res,nvar2+1,j,i);
            ss[j] += s*s * fabs(weyt[i]) ;
     } /* j loop */
   }   /* i loop */
}

void sort(double *v, int *a, int ii, int *jj) {
// ********************************************************
int i, ij, j, k, l, m, t, tt;
int il[20], iu[20];
double vt, vtt;

/* following two lines are only to convince lint that we aren't using
   uninitialized variables:                                             */
  tt = 123456789;
 vtt = 123456789.0;

/*     PUTS INTO A THE PERMUTATION VECTOR WHICH SORTS V INTO */
/*     INCREASING ORDER.  ONLY ELEMENTS FROM II TO JJ ARE CONSIDERED. */
/*     ARRAYS IU(K) AND IL(K) PERMIT SORTING UP TO 2**(K+1)-1 ELEMENTS */

/*     THIS IS A MODIFICATION OF CACM ALGORITHM #347 BY R. C. SINGLETON, */
/*     WHICH IS A MODIFIED HOARE QUICKSORT. */
// remains in F77=>C syntax
m = 1;
i = (ii);
j = (*jj);
L_10:
	if( i >= j )
		goto L_80;
L_20:
	k = i;
ij = (j+i)/2;
t = a[ij-1];
vt = v[ij-1];
if( v[i-1] <= vt )
	goto L_30;
a[ij-1] = a[i-1];
a[i-1] = t;
t = a[ij-1];
v[ij-1] = v[i-1];
v[i-1] = vt;
vt = v[ij-1];
L_30:
	l = j;
if( v[j-1] >= vt )
	goto L_50;
a[ij-1] = a[j-1];
a[j-1] = t;
t = a[ij-1];
v[ij-1] = v[j-1];
v[j-1] = vt;
vt = v[ij-1];
if( v[i-1] <= vt )
	goto L_50;
a[ij-1] = a[i-1];
a[i-1] = t;
t = a[ij-1];
v[ij-1] = v[i-1];
v[i-1] = vt;
vt = v[ij-1];
goto L_50;
L_40:
	a[l-1] = a[k-1];
a[k-1] = tt;
v[l-1] = v[k-1];
v[k-1] = vtt;
L_50:
	l = l - 1;
if( v[l-1] > vt )
	goto L_50;
tt = a[l-1];
vtt = v[l-1];
L_60:
	k = k + 1;
if( v[k-1] < vt )
	goto L_60;
if( k <= l )
	goto L_40;
if( l - i <= j - k )
	goto L_70;
il[m-1] = i;
iu[m-1] = l;
i = k;
m = m + 1;
goto L_90;
L_70:
	il[m-1] = k;
iu[m-1] = j;
j = l;
m = m + 1;
goto L_90;
L_80:
	m = m - 1;
if( m == 0 )
	return;
i = il[m-1];
j = iu[m-1];
L_90:
	if( j - i > 10 )
		goto L_20;
if( i == (ii) )
	goto L_10;
i = i - 1;
L_100:
	i = i + 1;
if( i == j )
	goto L_80;
t = a[i+1-1];
vt = v[i+1-1];
if( v[i-1] <= vt )
	goto L_100;
k = i;
L_110:
	a[k+1-1] = a[k-1];
v[k+1-1] = v[k-1];
k = k - 1;
if( vt < v[k-1] )
	goto L_110;
a[k+1-1] = t;
v[k+1-1] = vt;
goto L_100;
} /* end of func */

void qq( int *ppat, double *w, double *q1, int *perr) {
// ***************************************************

 *perr = 0;
 int npat = *ppat;
 double xnpat = 0. ;
 for (int i=0;i<npat;i++) xnpat += w[i];
 if (xnpat <= 0.) { *perr=1; return; }
 double u = -0.5 * w[0]/xnpat;
 for (int i=0;i<npat;i++) {
    u += w[i]/xnpat;
    double t = (u>0.5) ? sqrt( -2.0 * log(1.0-u)) : sqrt( -2.0 * log(u)) ;
    q1[i] = t - (2.515517 + 0.802853*t + 0.010328*t*t )
               /(1. + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t) ;
    if (u<=0.5) q1[i] =  -q1[i];
  }
}

void see(double *ss, double *press, int *pnss, int *pnpress, int *pnpat, int *pncomp) {
// *******************************************************

    double RAT = *pnpat - *pncomp - 1;
    if (RAT <= 0.) RAT = 1.0;
    RAT = *pnpat / RAT;
    for (int i = 0; ss && i < *pnss; i++) ss[i] = sqrt(ss[i] * RAT);
    for (int i = 0; press && i < *pnpress; i++) press[i] = sqrt(press[i] * RAT);
 }

double powi( double fptv, int intv)
// *********************************************************
{
    double t= 1.0;
 double a = abs( intv);
 for (int i=0; i<a; i++) t *= fptv ;
 if (intv<0) t = 1.0 / t;		/* zero divide! */
 return t;
}

/* f77c.h - definitions to support FORTRAN translations to C
			made by RTC PLUS from COBALT BLUE (Feb, 1987) */

#define mod(a,b)	((a) % (b))
#define lmod(a,b)	((a) % (b))

	/* externals to use in F77 std lib macros TO AVOID side effects */
	/* static int		f77_i_, f77_j_; */
	/* static long		f77_l_, f77_m_; */
	/* static double   	f77_d_, f77_e_; */


int nobootp(double *x, double *y, int *npat) {
// *****************************************************************

for(int i=0; i <= *npat; i+=1 ) x[i] = y[i];
return 1;
} 

int crossp(double *x, double *y, int *ix, int *npat, int *iout, int *iex, 
	int *ic) {
// **********************************************************
int ii;
int i, i1, i2;	/* changed from long for prototyping 7/13/95 */

/*  SET WEIGHTS FOR CROSS-VALIDATION */

int lmin2 = *iex < (*ic)  ? *iex : (*ic);
i1 = ((*ic)-1)*(*iout) + lmin2;
lmin2 = *iex < *ic ? *iex : *ic;
i2 = (*ic)*(*iout) + lmin2;

for(i=0; i < *npat ; i+=1 ){
	ii = ix[i];
	x[ii] = y[ii];
	if( i >= i1 && i <= i2 ) x[ii] =  - x[ii];
}

return 1;
} /* end of func */

int pcpls( int *nvar1, int *nvar2, int *npat, int *icomp, 
	int *nitmax, double *eps, double *x, double *y, double *weyt,
	double *w1, double *w2, double *b, double *ro, double *u,
	double *v, double *vt, int *ierr ) {
// *********************************************************
#define X(I_,J_)        (x+(I_)*( *nvar1)+(J_))
#define Y(I_,J_)        (y+(I_)*( *nvar2)+(J_))
#define W1(I_,J_)       (w1+(I_)*( *nvar1)+(J_))
#define W2(I_,J_)       (w2+(I_)*( *nvar2)+(J_))
#define B(I_,J_)        (b+(I_)*( *nvar1)+(J_))
#define U(I_,J_)        (u+(I_)*( *npat)+(J_))
#define V(I_,J_)        (v+(I_)*( *npat)+(J_))

/*     CALCULATES THE ICOMP-TH PLS COMPONENT */

(*ierr) = 0;

/*     Initialize V */

for(int i=0; i < *npat; i+=1 )
                 *V((*icomp),i) =  *Y(i,0);
int iter = 0;

    L_100:
	;
iter = iter + 1;

/*     Weights of first block */

for(int j=0; j < *nvar1; j+=1 ){
    double ss = 0.;
    for(int i=0; i < *npat; i+=1 ) if( weyt[i] >  0. )
	    ss = ss +  *X(i,j) * *V( *icomp, i) * weyt[i];
    *W1( *icomp, j ) = ss;
}
double ss = 0.;
for(int k=0; k < *nvar1; k+=1 )
    ss = ss + powi((double) *W1( *icomp, k ),2);
if( ss <= 0. ) { (*ierr) =  - 1; return 0; }

ss = sqrt(ss);
for(int k=0; k < *nvar1; k+=1 ) *W1( *icomp, k ) /=  ss;

/*     Latent variable of first block */
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. ) {
    ss = 0.;
    for(int k=0; k < *nvar1; k+=1 )
    	ss += *X(i,k) * *W1(*icomp,k);
    *U( *icomp, i) = ss;
}

/*     Norm of U */
double sumu = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
	sumu += powi((double) *U(*icomp,i),2);

/*     Weights of second block */
for(int j=0; j < *nvar2; j+=1 ){
    ss = 0.;
    for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
		ss += *Y(i,j) * *U(*icomp,i) *weyt[i];
    *W2( *icomp, j ) = ss;
}

ss = 0.;
for(int k=0; k < *nvar2; k+=1 )
    ss += powi((double) *W2(*icomp,k),2);

if( ss <= 0. ){ (*ierr) =  - 2; return 0; }
ss = sqrt(ss);
for(int k=0; k < *nvar2; k+=1 ) *W2(*icomp, k) /= ss;

/*     Latent variable of second block */
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. ) {
	ss = 0.;
	for(int k=0; k < *nvar2; k+=1 )
        ss += *Y(i, k) * *W2(*icomp, k);
    vt[i] = ss;
}
/*     convergence? */
ss = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
    ss += powi((double)(vt[i] - *V(*icomp,i)),2);
for(int i=0; i < *npat; i+=1 ) *V(*icomp,i) = vt[i];

if( ss > (*eps) && iter < (*nitmax) ) goto L_100;

/*     Inner relationship */
if( sumu <= 0. ){ (*ierr) =  - 3; return 0; }
ss = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
    ss += *U( *icomp, i ) * *V( *icomp, i ) * weyt[i];

ro[(*icomp)] = ss/sumu;

/*     Norm of the model vector */
sumu = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
    sumu += powi((double)( *U(*icomp,i) * ro[(*icomp)]),2) * weyt[i];

/*     Loadings */
if( sumu <= 0. ){ (*ierr) =  - 4; return 0; }
for(int j=0; j < *nvar1; j+=1 ){
    ss = 0.;
    for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
	    ss += *X( i, j ) * *U(*icomp, i) * ro[(*icomp)] * weyt[i]/sumu;
    *B(*icomp,j) = ss;
}

/*     Residuals */
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. ){
    ss = ro[(*icomp)] * *U(*icomp,i);
    for(int j=0; j < *nvar1; j+=1 )
	    *X(i,j) -= *B(*icomp,j) * ss;
    for(int j=0; j < *nvar2; j+=1 )
	    *Y(i,j) -= *W2(*icomp,j) * ss;
}

return 1;
#undef V
#undef U
#undef B
#undef W2
#undef W1
#undef Y
#undef X

} /* end of func */

int randomp(double *x, double *y, int *ix, int *npat) {
// ********************************************************
 int i, ii;
/*  RANDOM NUMBER GENERATOR FOR DRAWING BOOTSTRAP SAMPLE */

ranums(x,npat);
for(i=0; i < *npat; i+=1 ) ix[i] = (*npat)*((int)x[i]) + 1;

for(i=0; i < *npat; i+=1 ) x[i] = 0.;

for(i=0; i < *npat; i+=1 ){
	ii = ix[i];
	x[ii] = x[ii] + y[ii];
}
return 1;
} /* end of func */

int bootpls( double *ss, double *r2, double *press, double *cr2, double *off, 
	double *beta, int *nopt, int *ncomp, int *iboot, int *nvar2,
	int *nvar1 )
// ********************************************************
{

/*     MEAN AND AVERAGE OF THE BOOTSTRAPED QUANTITIES */

/*  ZERO THE STORAGES */
for(int j=0; j < *nvar2; j+=1 ) {
    if (ss) { 
	FORARR( ss, *nvar2, *iboot, j ) = 0.0;
	FORARR( r2, *nvar2, *iboot, j ) = 0.0;
	FORARR( ss, *nvar2, (*iboot)+1, j ) = 0.0;
	FORARR( r2, *nvar2, (*iboot)+1, j ) = 0.0;
    }
    for(int k=0; press && k < *ncomp; k+=1 ){
	press[((*iboot))*(*ncomp)*(*nvar2+1) + k*(*nvar2+1) + j] = 0.;
	press[((*iboot)+1)*(*ncomp)*(*nvar2+1) + k*(*nvar2+1) + j] = 0.;
	cr2[((*iboot))*(*ncomp)*(*nvar2) + k*(*nvar2+1) + j] = 0.;
	cr2[((*iboot)+1)*(*ncomp)*(*nvar2) + k*(*nvar2+1) + j] = 0.;
    }
}

nopt[(*iboot)] = 0;
nopt[(*iboot)+1] = 0.;
for(int j=0; j < *nvar2; j+=1 ){
	FORARR( off, *nvar2, *iboot, j ) = 0.0;
	FORARR( off, *nvar2, (*iboot)+1, j ) = 0.0;
	for(int i=0; i < *nvar1; i+=1 ) {
	    beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i] = 0.;
	    beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1) + i] = 0.;
	}
}

/*  SUM THE ELEMENTS OF ALL VECTORS INTO THE IBOOT+1-TH ELEMENTS */
for(int ib=0; ib < *iboot; ib+=1 ) for(int k=0; k < *nvar2; k+=1 ){
    if (ss) { 
	FORARR( ss, *nvar2, *iboot, k ) = FORARR( ss, *nvar2, *iboot, k ) 
		+ FORARR( ss, *nvar2, ib, k );
	FORARR( r2, *nvar2, *iboot, k ) = FORARR( r2, *nvar2, *iboot, k ) 
		+ FORARR( r2, *nvar2, ib, k );
    }
    for(int j=0; press && j < *ncomp; j+=1 ){
	cr2[((*iboot))*(*ncomp)*(*nvar2) + j*(*nvar2)+k] =  
	    cr2[((*iboot))*(*ncomp)*(*nvar2)+j*(*nvar2)+k] + 
	    cr2[(ib)*(*ncomp)*(*nvar2)*(*nvar2)+k];
	press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] = 
	    press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] +  
	    press[ib*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k];
    }	

    nopt[(*iboot)] = nopt[(*iboot)] + nopt[ib];
    for(int j=0; j < *nvar2; j+=1 ){
	FORARR( off, *nvar2, *iboot, j ) =
	    FORARR( off, *nvar2, *iboot, j ) + FORARR( off, *nvar2, ib, j);
	for(int i=0; i < *nvar1; i+=1 )
	    beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i] =  
		beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i] +  
		beta[ib*(*nvar2)*(*nvar1) + j*(*nvar1) +i];
    }
}

/*  CALCULATE AVERAGE */
double boo = (double) *iboot;
for(int k=0; k < *nvar2; k+=1 ){
    if (ss) {
	FORARR(ss, *nvar2, *iboot, k) /= boo;
	FORARR(r2, *nvar2, *iboot, k) /= boo;
    }
    for(int j=0; press && j < *ncomp; j+=1 ){
	press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] /= boo;
	cr2[((*iboot))*(*ncomp)*(*nvar2)+j*(*nvar2)+k] /= boo;
    }
}

nopt[(*iboot)] /= *iboot;
for(int j=0; j < *nvar2; j+=1 ){
    FORARR( off, *nvar2, *iboot, j ) /= boo;
    for(int i=0; i < *nvar1; i+=1 )
	beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) +i] /= boo;
}

/*  SUM THE SQUARED DEVIATION FROM THE AVERAGES IN THE IBOOT+2-TH ELEMENT */
for(int ib=0; ib < *iboot; ib+=1 ) {
  for(int k=0; k < *nvar2; k+=1 ){
    if (ss) { 
	FORARR( ss, *nvar2, *iboot+1, k ) =
 	    FORARR( ss, *nvar2, *iboot+1, k ) +
powi((double)( FORARR( ss, *nvar2, ib, k ) - FORARR( ss, *nvar2, *iboot, k)), 2);
	FORARR( ss, *nvar2, *iboot+1, k ) =
 	    FORARR( ss, *nvar2, *iboot+1, k ) +
powi((double)( FORARR( ss, *nvar2, ib, k ) - FORARR( ss, *nvar2, *iboot, k)), 2);
    }
    for(int j=0; press && j < *ncomp; j+=1 ){
	cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k] =  
	    cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k] + 
	    powi((double)( cr2[ib*(*ncomp)*(*nvar2)+j*(*nvar2)+k] -
	       cr2[((*iboot))*(*ncomp)*(*nvar2)+j*(*nvar2)+k]),2);
	press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] =  
	    press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] + 
	    powi((double)( press[ib*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] - 
		press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k]),2);
    }
  }
  nopt[(*iboot)+1] = nopt[(*iboot)+1] + (int) powi((double)(nopt[ib] -nopt[(*iboot)]),2);
  for(int j=0; j < *nvar2; j+=1 ){
	FORARR(off, *nvar2, *iboot+1, j) =
	    FORARR( off, *nvar2, *iboot+1, j ) +
	    powi((double)( FORARR(off,*nvar2,ib,j) - FORARR(off,*nvar2,*iboot,j)),2);
	for(int i=0; i < *nvar1; i+=1 )
	    beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1) + i] = 
    		beta[((*iboot)+1)*(*nvar1)*(*nvar2) + j*(*nvar1) + i] + 
    		powi((double)( beta[ib*(*nvar2)*(*nvar1) + j*(*nvar1) + i] - 
		    beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i]),2);
  }
}

/*  CALCULATE VARIANCE */
for(int k=0; k < *nvar2; k+=1 ){
   if (ss) {
	FORARR( ss, *nvar2, *iboot+1, k )
	     = sqrt( FORARR( ss, *nvar2, *iboot+1, k )/boo);
	FORARR( r2, *nvar2, *iboot+1, k ) =
	    sqrt( FORARR( r2, *nvar2, *iboot+1, k )/boo );
   }
   for(int j=0; press && j < *ncomp; j+=1 ){
	press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] = 
	    sqrt( press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k]/boo);
	cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k] = 
	    sqrt( cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k]/boo);
   }
}

nopt[(*iboot)+1] = (int)(sqrt(float(nopt[(*iboot)+1])/boo));
for(int j=0; j < *nvar2; j+=1 ){
   FORARR( off, *nvar2, *iboot+1, j ) = 
	sqrt( FORARR( off, *nvar2, *iboot+1, j )/boo);
   for(int i=0; i <= *nvar1; i+=1 )
	beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1)+i] = 
	    sqrt( beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1) + i]/boo);
}

return 1;
} /* end of func */


int plsjer( int *pnpat, int *pnvar1, int *pnvar2, int *pncomp, int *nitmax, 
	double *eps, int *iboot, int *icros, int *icent, double *x, double *y, 
	double *weyt, double *weytb, double *xbar, double *xscal, 
	double *off, double *beta, double *varnc, double *w1, double *w2, 
	double *b, double *u, double *v, double *ro, double *ss, double *r2, 
	double *press, double *cr2, int *nopt, double *ypred, double *res, 
	double *sss, double *ssy, int *ierr, double *scrtch, int *ix )
{
#define X(I_,J_)	(x+(I_)*(nvar1)+(J_))
#define Y(I_,J_)	(y+(I_)*(nvar2)+(J_))
#define OFF(I_,J_)	(off+(I_)*(nvar2)+(J_))
#define BETA(I_,J_,K_)	(beta+(I_)*(nvar2)*(nvar1)+(J_)*(nvar1)+(K_))
#define VARNC(I_,J_)	(varnc+(I_)*(nvar2)+(J_))
#define W1(I_,J_)	(w1+(I_)*(nvar1)+(J_))
#define W2(I_,J_)	(w2+(I_)*(nvar2)+(J_))
#define B(I_,J_)	(b+(I_)*(nvar1)+(J_))
#define U(I_,J_)	(u+(I_)*(npat)+(J_))
#define V(I_,J_)	(v+(I_)*(npat)+(J_))
#define SS(I_,J_)	(ss+(I_)*(nvar2)+(J_))
#define R2(I_,J_)	(r2+(I_)*(nvar2)+(J_))
#define PRESS(I_,J_,K_)	(press+(I_)*(ncomp)*(nvar2+1)+(J_)*(nvar2+1)+(K_))
#define CR2(I_,J_,K_)	(cr2+(I_)*(ncomp)*(nvar2)+(J_)*(nvar2)+(K_))
#define YPRED(I_,J_)	(ypred+(I_)*(nvar2)+(J_))
#define RES(I_,J_)	(res+(I_)*(nvar2+1)+(J_))
 /*
int npat,nvar1,nvar2,ncomp; 	// changed from long for prototyping 7/13/95
int i, ib, iboot1, ic, icomp, icros1, iex, iout, j, jj, kl, l,
	nn, novariance ;	// changed from long for prototyping 7/13/95
double atl, pmin,      s, sum, varmax, xnpat, xnpatt;
*/
/*     MAIN PLS DRIVER */

/* NPAT,NVAR1,NVAR2,NCOMP, */

int npat = *pnpat;
    int nvar1= *pnvar1;
    int nvar2 = *pnvar2;
    int ncomp = *pncomp;

(*ierr) = 0;
int icros1 = 0;
int iboot1 = 0;

int nn = npat*(nvar1 +nvar2 +1) + 1;

/*  calculate quantiles for the q-q plot */
qq(&npat,weyt,scrtch,ierr);
if( (*ierr) != 0 ) return 0;

for(int i=0; i < npat; i+=1 ) *RES(i,nvar2) = scrtch[i];

/*  set flags for no cross-validation and no bootstrap */
if( (*icros) == 0 ){ (*icros) = 1; icros1 = 1; }
if( (*iboot) == 0 ){ (*iboot) = 1; iboot1 = 1; }

/*  bootstrap monster loop */
for(int ib=0; ib < *iboot; ib+=1 ) {
    for(int icomp=0; icomp < ncomp; icomp+=1 )
	for(int j=0; press && j < nvar2+1; j+=1 )
	    press[ib*ncomp*(nvar2+1)+(icomp)*(nvar2+1)+j] = 0.;

	/*  no bootstrap, copy the weights */
	if( iboot1 == 1 ) for (int i = 0; i < npat; i++) weytb[i] = weyt[i];
	/*  draw a bootstrap sample , new weights in weytb */
	    else randomp(weytb,weyt,ix,&npat);

	/*  set up an auxiliary vector for cross-validation pointing to non zero weight */
    int kl = 0;
	if( icros1 != 1 ){
	    for(int i=0; i < npat; i+=1 ) if( weytb[i] > 0. ) {
			kl = kl + 1;
			ix[kl-1] = i;
	    }
	    ranums(scrtch,&kl);
	    sort(scrtch,ix,1,&kl);
	    for(int i=0; i < npat; i+=1 ) scrtch[i] = 0.;

		/*  calculate how many samples to delete */
	}
    int iout = kl/(*icros);
    int iex = lmod(kl,(*icros));

	/*  calculate sum of squares */
    double atl = 0.0;
    double xnpatt = 0.;
	for(int j=0; j < nvar2; j+=1 ){
	    ssy[j] = 0.;
	    for(int i=0; i < npat; i+=1 ){
		    xnpatt = xnpatt + weytb[i];
		    atl = atl +  *Y( i, j ) *weytb[i];
        }
	    atl = atl/xnpatt;
	    for(int i=0; i < npat; i+=1 )
		    ssy[j] += powi((double)( *Y(i,j) - atl ),2) *weytb[i];
	}

	/*  cross-validation loop */
	for(int ic=0; ic < *icros; ic+=1 ){
			/*  no cross-validation, copy the weights */
	    if( icros1 == 1 ) for (int i = 0; i < npat; i++) scrtch[i] = weyt[i];
		/*  set weights for cross-validation in scrtch */
		else crossp(scrtch,weytb,ix,&kl,&iout,&iex,&ic);

		/*  SUM OF WEIGHTS */
	    double xnpat = 0.;
	    for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. ) xnpat += scrtch[i];
	    if( xnpat <= 0. ){ (*ierr) = 2; return 0; }

		/*  MEAN VALUES AND SCALE */
		/*  1 - X;  2 - X-XBAR;   3 - (X-XBAR)/XSCAL;    4 - X/XSCAL; */
		/*  CALCULATE XBAR THE MEAN */
	    if( (*icent) > 1 ){
		    for(int j=0; j < nvar1; j+=1 ){
		        double s = 0.;
		        for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			        s += *X(i,j) * scrtch[i];
		        xbar[j] = s/xnpat;
		    }
		    for(int j=0; j < nvar2; ++j) {
		        double s = 0.;
		        for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			        s += *Y(i,j) * scrtch[i];
		        int jj = j + nvar1;
		        xbar[jj] = s/xnpat;
		    }
	    }
		/*   CALCULATE XSCAL THE VARIANCE */
	    if( (*icent) >= 3 ) {
		    int novariance = 0;
		    for(int j=0; j < nvar1; j+=1 ){
		        double s = 0.;
		        for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			        s += powi((double)(*X(i,j)- xbar[j]),2)*scrtch[i];
                xscal[j] = sqrt(s/(xnpat -1));
                if( s <= 0. ) novariance = novariance + 1;
            }
            if( novariance != 0 ) {
                double varmax = 1.e20;
                for(int j=0; j < nvar1; j+=1 )
                    if( xscal[j] > 0.0 && xscal[j] < varmax )
				            varmax = xscal[j];
                if( varmax == 1.e20 ){ (*ierr) = 2; return 0; }
                varmax = varmax/1000.;
                for(int j=0; j < nvar1; j+=1 ) if( xscal[j] <= 0.0 )
				        xscal[j] = varmax;
            }
            for(int j=0; j < nvar2; j+=1 ){
                int jj = j + nvar1;
                double s = 0.;
                for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			            s += powi((double)( *Y(i,j) - xbar[jj]), 2) * scrtch[i];
                xscal[jj] = sqrt(s/(xnpat -1));
                if( xscal[j] == 0 ){ (*ierr) = 2; return 0; }
            }
        }
		/*  ICENT = 1 NO CENTERING OR SCALING */
	    if( (*icent) == 1 ){
		    for(int i=0; i < npat; i+=1 ){
		        for(int j=0; j < nvar1; j+=1 ) {
                    int l = i*nvar1 + j +npat;
                    scrtch[l] = *X(i,j);
                }
		        for(int j=0; j < nvar2; j+=1 ){
			        int l = npat + npat*nvar1 + (i )*nvar2 + j;
			        scrtch[l] =  y[(i)*nvar2+j];
		        }
		    }

		    for(int j=0; j < nvar1+nvar2; j+=1 ){
		        xbar[j] = 0.;
		        xscal[j] = 1.;
		    }
	    }
		/*  ICENT = 2 */
	    if((*icent) == 2 ) for(int i=0; i < npat; i+=1 ){
		    for(int j=0; j < nvar1; j+=1 ){
			    int l = (i )*nvar1 + j + npat;
			    scrtch[l] =  *X(i,j) - xbar[j];
		    }
		    for(int xj=0; xj < nvar2; xj+=1 ){
			    int jj = xj + nvar1;
			    int xl = npat + npat*nvar1 + (i )*nvar2 + xj;
			    scrtch[xl] =  *Y(i,xj) - xbar[jj];
		    }
		    for(int j=0; j < nvar1+nvar2; j+=1 ) xscal[j] = 1.;
	    }
		/*  ICENT = 3 */
	    if( (*icent) == 3 ) for(int i=0; i < npat; i+=1 ){
		    for(int j=0; j < nvar1; j+=1 ){
			    int l = (i )*nvar1 + j + npat;
			    scrtch[l] = ( *X(i,j) - xbar[j])/xscal[j];
		    }
		    for(int j=0; j < nvar2; j+=1 ){
			    int jj = j + nvar1;
			    int l = npat + npat*nvar1 + (i )*nvar2 + j;
			    scrtch[l] = ( *Y(i,j) - xbar[jj])/xscal[jj];
		    }
	    }
		/*  ICENT = 4 */
	    if( (*icent) == 4 ){
		    for(int i=0; i < npat; i+=1 ){
		        for(int j=0; j < nvar1; j+=1 ){
			        int l = (i )*nvar1 + j + npat;
			        scrtch[l] =  *X(i,j)/ xscal[j];
		        }
		        for(int j=0; j < nvar2; j+=1 ){
		    	    int jj = j + nvar1;
		    	    int l = npat + npat*nvar1 + (i )*nvar2 + j;
		    	    scrtch[l] =  *Y(i,j) / xscal[jj];
		        }
		    }
		    for(int j=0; j < nvar1+nvar2; j+=1 ) xbar[j] = 0.;
	    }

		/*  CALCULATE MODEL OF NCOMP COMPONENTS */
	    for(int icomp=0; icomp < ncomp; icomp+=1 ){
		    if (!pcpls(&nvar1,&nvar2,&npat,&icomp,nitmax,eps,&scrtch[npat],
			  &scrtch[npat*(nvar1+1)],scrtch,w1,w2,b,
			  ro,u,v,&scrtch[nn-1],ierr)) return 0;
		    if( (*ierr) != 0 ){
		        int ncomp = icomp - 1;
		        if( ncomp < 0 ) return 0;
		        (*ierr) = 0;
		        std::cout << "No more than " << ncomp << "components can be calculated\n";
		        goto L_222;
		    }
		    for(int j=0; j < nvar2; j+=1 ){
		        double sum = 0.;
		        int jj = j + nvar1;
		        for(int i=0; i < npat; i+=1 ) {
			        int l = npat + npat*nvar1 + (i )*nvar2 + j;
			        sum += powi((double)scrtch[l],2)*scrtch[i];
		        }
		        if( (*icent) >= 3 ) {
			        sum *= powi((double)xscal[jj],2);
			        *VARNC(icomp,j) = 1. - sum/ssy[j];
		        }
		    }
	    }
		L_222:
			;
		/*  load the data again for residual calculation */
	    for(int i=0; i < npat; i+=1 ){
		    for(int j=0; j < nvar1; j+=1 ){
		        int l = (i )*nvar1 + j + npat;
		        scrtch[l] = *X(i,j );
		    }
		    for(int j=0; j < nvar2; j+=1 ){
		        int l = (i )*nvar2 + j + npat*(nvar1 +1);
		        scrtch[l] =  *Y( i, j );
		    }
	    }
	/*  samples with negative weights are predicted, so set all of them to - */
	    if( icros1 == 1 ) for(int i=0; i < npat; i+=1 ) scrtch[i] =  - scrtch[i];

		/*  REGRESSION MATRIX AND PREDICTION */
	    for(int icomp=0; icomp < ncomp; icomp+=1 ){
		    if( iboot1 != 1 ){
			    matrix(&nvar1,&nvar2,&icomp,&beta[(ib)*nvar1*nvar2 ],
				  w1,w2,b,ro,&scrtch[nn-1]);
			    resid(&npat,&nvar1,&nvar2,&scrtch[npat],&scrtch[npat*(nvar1+1)+1-1],
                      scrtch,xbar,xscal,&off[nvar2*ib],
                      &beta[(ib)*nvar1*nvar2],res,ypred,sss,ierr);
		    }
		    else{
			    matrix(&nvar1,&nvar2,&icomp,&beta[(ic)*nvar2*nvar1],
				  w1,w2,b,ro,&scrtch[nn-1]);
			    resid(&npat,&nvar1,&nvar2,&scrtch[npat],
                 &scrtch[npat*(nvar1+1)],scrtch,xbar,xscal,&off[nvar2*ic],
                 &beta[(ic)*nvar1*nvar2],res,ypred,sss,ierr);
		    }
		    if( (*ierr) != 0 ) return 0;
		    if( icros1 == 0 )
				/*  sum of cross-validated squared residuals */
		        for(int j=0; press && j < nvar2; j+=1 )
			        *PRESS(ib,icomp,j) +=  sss[j];
	    }
	    if( icros1 == 1 ){
			/* sum of squared residuals and r squared */
		    for(int j=0;  j < nvar2; j+=1 ){
		        *SS(ib,j) = sss[j];
		        *R2(ib,j) = 1. - *SS(ib,j) /ssy[j];
		    }
		    for(int j=0; j < nvar2; j+=1 ) *SS( ib,j ) /= xnpatt;
	    }
	}

	/*  find the optimal number of components */
	if( icros1 != 1 ){
	    double pmin = 1.e10;
	    for(int icomp=0; icomp < ncomp; icomp+=1 ){
		double s = 0.;
		for(int j=0; j < nvar2; j+=1 )
		    s += *PRESS(ib,icomp,j);
            *PRESS(ib,icomp,nvar2+1) = s/( (double) nvar2);
            if (*PRESS(ib,icomp,nvar2+1) < pmin ){
				pmin = *PRESS(ib,icomp,nvar2+1);
				nopt[ib] = icomp;
		    }
	    }
	}
	/*  CROSS-VALIDATED R SQUARED */
	if( icros1 != 1 ){
	    for(int icomp=0; icomp < ncomp; icomp+=1 )
		    for(int j=0; j < nvar2; j+=1 )
		        *CR2(ib,icomp,j) = 1. - *PRESS(ib,icomp,j) /( ssy[ j]);
	        for(int j=0; j < nvar2+1; j+=1 ){
		        double undo;
		        if (npat-ncomp > 0) undo = (double) npat-ncomp;
			        else {undo = 1.0; ncomp = npat - 2; }
		        for(int icomp=0; icomp < ncomp; icomp+=1 )
                    *PRESS(ib,icomp,j) = *PRESS(ib,icomp,j) /xnpatt * ( undo / (double)(npat-icomp));
	        }
    }
}
if( iboot1 == 1 ) (*iboot) = 0;
if( icros1 == 1 ) (*icros) = 0;
/*  calculate the mean and the variance of the bootstrapped quantities */
if( iboot1 != 1 )
	bootpls(ss,r2,press,cr2,off,beta,nopt,&ncomp,iboot,&nvar2,&nvar1);

int tmpa_55 = ((*iboot) +2)*nvar2;
    int tmpa_56 = (nvar2 +1)*ncomp*((*iboot) +2);
see(ss,press,&tmpa_55,&tmpa_56,&npat,&ncomp);	/*SEE (sqrt) *///

return 1;
#undef RES
#undef YPRED
#undef CR2
#undef PRESS
#undef R2
#undef SS
#undef V
#undef U
#undef B
#undef W2
#undef W1
#undef VARNC
#undef BETA
#undef OFF
#undef Y
#undef X
} /* end of func */

/*
int main() {
// for testing
   double xs[9] = {3.0,2.0,3.0, 4.0,5.0,5.0, 6.0,7.0,6.0};
   double ys[3] = {3.0,4.0,5.0};
   int nOKMol = 3;
   int nOKfVal = 3;

   double ywts[nOKMol];
   int ierr;
   int nys = 1;
   int mxncomp = 2;
   int maxiter = 100;
   double eps =1.0e-4;
   int iboot = 0;
   int icros = 0;
   int icent = 3;
   double *weytb = new double[nOKMol];
   double *xbar = new double[nOKMol+1];
   double *xscal = new double[nOKMol+1];
   double *intcpt = new double[3];
   double *coeff = new double[2*nOKfVal];
   double *varnc = new double[mxncomp*nOKfVal];
   double *wtx = new double[mxncomp*nOKfVal];
   double *wty = new double[mxncomp];
   double *loadings = new double[mxncomp*nOKfVal];
   double *latentx = new double[mxncomp*nOKMol];
   double *latenty = new double[mxncomp*nOKMol];
   double *inner = new double[mxncomp];
   double *ssqRsdl = new double[mxncomp*6];
   double *r2 = new double[2];
   double *sdep = new double[nOKfVal*mxncomp];
   double *q2 = new double[2*mxncomp];
   int *optncomp = new int[2];
   double *ypred = new double[nOKMol];
   double *resdl = new double[2*nOKMol];
   double *sss = new double[2];
   double *ssy = new double[ nOKMol ];
   int *iscratch = new int[nOKMol];
   double *scratch = new double[nOKMol * (nOKfVal + 2) + (nOKfVal+1)];
   for (int iy = 0; iy < 3; ++iy) intcpt[iy] = 0.0;
   for (int iy = 0; iy < 2; ++iy) 
	{optncomp[iy] = 0; r2[iy]=sss[iy]=ssqRsdl[iy]=0.0;}
   for (int iy = 0; iy < nOKMol; ++iy) {iscratch[iy] = 0; 
	ssy[iy] = weytb[iy] = ypred[iy] = 0.0; ywts[iy] = 1.0;}
   for (int iy = 0; iy < 2*nOKMol; ++iy) resdl[iy] = 0.0;
   for (int iy = 0; iy < 2*nOKfVal; ++iy) coeff[iy] = 0.0; 
   for (int iy = 0; iy < mxncomp*nOKfVal; ++iy) 
	varnc[iy] = wtx[iy] = loadings[iy] = sdep[iy] = 0.0;
   for (int iy = 0; iy < mxncomp; ++iy) wty[iy] = inner[iy] = 0.0;
   for (int iy = 0; iy < mxncomp*nOKMol; ++iy) latentx[iy] = latenty[iy] = 0.0;
std::cout << nOKMol << ' ' << nOKfVal << ' ' << nOKMol * (nOKfVal + 2) + (nOKfVal+1) << std::endl;

   plsjer( &nOKMol, &nOKfVal, &nys, &mxncomp, &maxiter, &eps, &iboot,
        &icros, &icent, xs, ys, ywts, weytb, xbar, xscal, intcpt, coeff,
        varnc, wtx, wty, loadings, latentx, latenty, inner, ssqRsdl, r2,
           NULL, NULL, optncomp, ypred, resdl, sss, ssy, &ierr, scratch, iscratch);
std::cout << *r2 << ' ' << *ssqRsdl << "\n";
}

1201447-30072015
00001ZY40aiTiyS7T31eq1dj1m5ERE
Ml2EhOhukYoBZ47EDq5YKURQrGHC0F
NYhaxpIWRlGLX2UnxDWG8gVraYn7r0
*/

float rbAttn = 0.85f;
//ofstream TClog ;

void wcsv( std::string fn, double *vals, int nrow, int ncol ) {
    ofstream f;
    f.open(fn);
    f.precision(5);
    for (int r = 0; r < nrow; ++r) {
	for (int c = 0; c < ncol; ++c) {
	    if (c > 0) f << ',';
	    f << vals[r*ncol+c];
	}
	f << endl;
    }
    f.close();
}

void wcsvF( std::string fn, float *vals, int nrow, int ncol ) {
    ofstream f;
    f.open(fn);
    f.precision(5);
    float *v = vals;
    for (int r = 0; r < nrow; ++r) {
	for (int c = 0; c < ncol; ++c) {
	    if (c > 0) f << ',';
	    f << *v;
	    v++;
	}
	f << endl;
    }
    f.close();
}

void grid2txt( string fname, float *fld, int *npts, float *lo )
// ******************************************
{
   
  ofstream gcout;
  gcout.open(fname);

   int xdim = npts[0];
   int ydim = npts[1];
   int zdim = npts[2];
  char buffer[80];

  int ip = 0;
  for (unsigned int iz=0; iz<zdim; ++iz)
    for (unsigned int iy=0; iy<ydim; ++iy)
      for (unsigned int ix=0; ix<xdim; ++ix) {
          float x = lo[0] + 2.0f * ((float)ix ) - 1.0f;
          float y = lo[0] + 2.0f * ((float)iy ) - 1.0f;
          float z = lo[0] + 2.0f * ((float)iz ) - 1.0f;
        sprintf(buffer,"%6.2f %6.2f %6.2f %-12.6f\n", x, y, z, fld[ip] );
	    gcout << buffer ;
        ip++;
      }
  gcout.close();
}

void savemodel(string modelfn, int nFpt, double intcpt, float *beta,
               int *npts, float *lo,  int ncomp, double sdep, double q2,
            double r2, double s) {
// ********************************************
    ofstream f;
    f.open(modelfn);
    f << "Stats: #comp=" << ncomp << " (q2=" << q2 << ";SDEP=" << sdep << ") r2=" << r2 << "; s=" << s << endl;
    for (int i = 0; i < 3; i++) f << npts[i] << endl;
    for (int i = 0; i < 3; i++) f << lo[i] << endl;

    f << intcpt << endl;
    for (int i = 0; i < nFpt; i++) f << beta[i] << endl;
    f.close();
}

void addRBAtm( OEGraphMol mol, int nowAtm, bool *atAdded, vector<int> &NextLyr, 
	vector<int> &NextLyrFrom, bool *visited, double *wts, double wt ) {
// ******************************************

//cout << nowAtm << " add\n";
  for (OEIter<OEAtomBase>
        nbor = mol.GetAtom(OEHasAtomIdx((unsigned int) nowAtm))->GetAtoms(); nbor; ++nbor) {
// cout << ' ' << nbor->GetIdx() ;
     if (visited[nbor->GetIdx()]) continue;
     *atAdded = true;
     wts[nbor->GetIdx()] = wt;
     if (nbor->GetDegree() == 1) continue;
     OEBondBase *cBond = nbor->GetBond(mol.GetAtom(OEHasAtomIdx((unsigned int) nowAtm)));
     if (!cBond->IsInRing() && cBond->GetOrder() == 1 && 
		    !visited[(int) nbor->GetIdx()] ) {
         bool newnxt = true;
	    if ((int) NextLyr.size() > 0)
	    for (int i = 0; i < (int) NextLyr.size(); ++i)
		    if (NextLyr[i] == (int) nbor->GetIdx()) { newnxt = false; break;}
	    if (newnxt) {
	        NextLyr.push_back(nbor->GetIdx());
	        NextLyrFrom.push_back(nowAtm);
	    }
     }
     else {
	    visited[nbor->GetIdx()] = true;
	    wts[nbor->GetIdx()] = wt;
        addRBAtm(mol,nbor->GetIdx(),atAdded,NextLyr,NextLyrFrom,visited,wts,wt);
     }
  } 
}

void getRbAtWts( double *wts, OEGraphMol mol ) {
// ********************************************
  double currWt = 1.0;
  bool visited[(int) mol.NumAtoms() ];
  for (int i = 0; i < (int) mol.NumAtoms(); ++i) visited[i] = false;
  vector<int> NowLayer;
  vector<int> NextLayer;
  vector<int> NextLayerFrom;

  int rootAt = 0;
  OEStringToNumber(OEGetSDData(mol,"AnchorAts"), rootAt);

  wts[rootAt] = currWt; 
  visited[rootAt] = true;
  NowLayer.push_back(rootAt);
  while (true) {
	bool atAdded = false;
	NextLayer.clear();
	NextLayerFrom.clear();
	for (int i = 0; i < (int) NowLayer.size(); ++i) 
	    addRBAtm(mol,NowLayer[i],&atAdded,NextLayer,NextLayerFrom,visited,wts,currWt); 
	if (!atAdded) break;
	NowLayer.clear();
	for (int i = 0; i < (int) NextLayer.size(); ++i) {
		NowLayer.push_back( NextLayer[i]);
		visited[ NextLayer[i]] = true;
		wts[ NextLayer[i]] = currWt;
	}
	currWt *= rbAttn;
    NextLayer.clear();
  }
}
static double scutoff[16] = {9999.,   0.,   2.,   4.,   6.,   8.,  10.,  12., 14.,  16.,  18.,  20.,  22.,  24.,  26.,  30.  };
static double scutValues[16] = {9999.,   -0.1,   1.,   3.,   5.,   7.,  9.,  11., 13.,  15.,  17., 19.,  21.,  23.,  25.,  29.  };

static double ecutoff[16] = {9999.,   -12.9999, -10.9999, -8.9999, -6.9999, -4.9999, -2.9999, -0.9999, 1.0001,   3.0001,   5.0001,   7.0001,   9.0001,  11.0001,  13.0001, 15.0001 };
static double ecutValues[16] = {9999.,   -13.9999, -11.9999, -9.9999, -7.9999 , -5.9999, -3.9999, -1.9999, 0.0001,   2.0001,   4.0001,   6.0001,   8.0001,  10.0001,  12.0001, 14.0001 };

static double LookupBinnedValue(double value, int isElectroStatic) {
// **********************************************************
double *cutoff ;
double *cutoffValues ;
int i ;

   if ( isElectroStatic ) {
       cutoff = ecutoff ;
       cutoffValues = ecutValues ;
    }
    else {
       cutoff = scutoff ;
       cutoffValues = scutValues ;
    }

    for ( i = 1 ; i < 16 ; i++ )
        if ( value <= cutoff[i] )
            return cutoffValues[i] ;

    return 0.00 ;
}

bool calcField( OEGraphMol mol, float *ster, float *elec, int *nstep, float *basept,
	int molID ) {
// *****************************************************
#define STERIC_MAX 30.0f
#define Q2KC 332.0f
#define MIN_SQ_DISTANCE 1.0e-4
#define ELECTRO_MAX 14.0f
#define ELECTRO_MIN -14.0f

  double rbAtWts[ (int) mol.NumAtoms() ];
//    for (unsigned int i = 0; i < mol.NumAtoms(); i++) rbAtWts[i] = 0.0;
  getRbAtWts( rbAtWts, mol );

  double AtChg[ (int) mol.NumAtoms() ];
  if (!OEGasteigerPartialCharges(mol)) {
	cout << "Partial charge calculations failed for structure " << molID << endl;
	return false;
  }
  for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom)
        AtChg[atom->GetIdx()] = atom->GetPartialCharge() ;

  double vdwA[ (int) mol.NumAtoms() ];
  double vdwB[ (int) mol.NumAtoms() ];

  for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
	double radius = 0.0, eps = 0.0;
	int anum[12] = {1,6,7,8,9,11,14,15,16,17,35,53};
	double epsval[12] = 
{0.042,0.107,0.095,0.116,0.400,0.0,0.042,0.314,0.314,0.4,0.434,0.6};
	double radval[12] = 
{1.5,1.7,1.55,1.52,1.47,1.2,2.1,1.8,1.8,1.75,1.85,1.98};
	bool haveAW = false;
	for (int naw = 0; naw < 12; ++naw) if (((int)atom->GetAtomicNum()) == anum[naw]) {
	    haveAW = true;
	    radius = radval[naw];
	    eps = epsval[naw];
	}
	if (!haveAW) {radius = 1.7; eps = .107;}
	radius += 1.7;
	eps = sqrt( eps * .107 );
	vdwA[ atom->GetIdx() ] = 2.0 * eps * rbAtWts[ atom->GetIdx() ] * 
		pow( radius, 6.0 );
	vdwB[ atom->GetIdx() ] = eps * rbAtWts[ atom->GetIdx() ] * 
		pow( radius, 12.0 );
  }
  double coo[ 3 * (int) mol.NumAtoms() ];
  OEGetPackedCoords(mol,coo);

  int nstepx = nstep[0];
  int nstepy = nstep[1];
  int nstepz = nstep[2];
  double stepx = 2.0;
  double stepy = 2.0;
  double stepz = 2.0;
  double lowx = basept[0];
  double lowy = basept[1];
  double lowz = basept[2];
  double *coord = coo;
  int max_steps = (int) (4.0 / stepx);
  if ( max_steps <= 0 || ((double) max_steps * stepx ) < 4.0 ) max_steps += 1;
  for (int nat = 0; nat < (int) mol.NumAtoms(); ++nat) {
        double curr_x = *coord;
        double curr_y = *(coord+1);
        double curr_z = *(coord+2);
        coord += 3;

        int iz = (int) ( fabs(curr_z - lowz + 0.5) / stepz);
        int iy = (int) ( fabs(curr_y - lowy + 0.5) / stepy);
        int ix = (int) ( fabs(curr_x - lowx + 0.5) / stepx);

        int curr_iz = iz - max_steps;
        int curr_iy = iy - max_steps;
        int curr_ix = ix - max_steps;
        int curr_nstepsz = iz + max_steps + 1;
        int curr_nstepsy = iy + max_steps + 1;
        int curr_nstepsx = ix + max_steps + 1;
        if ( curr_iz < 0 ) curr_iz = 0;
        if ( curr_iy < 0 ) curr_iy = 0;
        if ( curr_ix < 0 ) curr_ix = 0;
        if ( curr_iz >= nstepz ) curr_iz = nstepz - 1;
        if ( curr_iy >= nstepy ) curr_iy = nstepy - 1;
        if ( curr_ix >= nstepx ) curr_ix = nstepx - 1;
        if ( curr_nstepsz > nstepz ) curr_nstepsz = nstepz;
        if ( curr_nstepsy > nstepy ) curr_nstepsy = nstepy;
        if ( curr_nstepsx > nstepx ) curr_nstepsx = nstepx;
        double maxw = STERIC_MAX * rbAtWts[ nat ];
	double z;
	double y;
	double x;
	int st;
        for ( iz=0, z=lowz; iz < nstepz; iz++, z += stepz ) {
            double zd = z - curr_z;
            zd = zd*zd;
            for (iy=0, y=lowy; iy < nstepy; iy++, y += stepy ) {
                double yd = y - curr_y;
                yd = yd*yd;
                for (ix=0, x=lowx; ix < nstepx; ix++, x += stepx ) {
                    st = ( iz * nstepy * nstepx ) + (iy * nstepx ) + ix;
                    double sum_steric = (double) ster[st];
                    double xd = x - curr_x;
                    double dis2 = xd*xd + yd + zd;
                    if (iz >= curr_iz && iz < curr_nstepsz &&
			iy>=curr_iy && iy < curr_nstepsy &&
                            ix >= curr_ix && ix < curr_nstepsx) {
		      double atm_steric;
                      if ( dis2 >= MIN_SQ_DISTANCE ) {
                        double dis6 = dis2 * dis2 * dis2;
                        double dis12= dis6 * dis6 ;
                        atm_steric = vdwB[nat]/dis12 - vdwA[nat]/dis6;
                        if (  atm_steric > maxw ) atm_steric = maxw;
                      }
                      else atm_steric = maxw;
                      sum_steric += atm_steric ;
                      ster[st] = sum_steric > STERIC_MAX ?
				STERIC_MAX : (float) sum_steric;
                    }
                    if ( dis2 == 0.0 ) dis2 = 0.000001;
                    elec[st] += (float) (AtChg[nat] / dis2) * rbAtWts[nat] ;
                    if (nat == (int) mol.NumAtoms() - 1) {
                        elec[st] *= Q2KC;
                        if (elec[st] > ELECTRO_MAX) elec[st] = ELECTRO_MAX;
                           else if (elec[st] < ELECTRO_MIN) elec[st] = ELECTRO_MIN;
                    }
		}  // X
	    } // Y
	} // Z
    } //atom

    for (int np = 0; np < nstepx * nstepy * nstepz; np ++ ) {
        ster[np] = (float) LookupBinnedValue( (double) ster[np], false );
        elec[np] = (float) LookupBinnedValue( (double) elec[np], true);
    }

    grid2txt(to_string(molID)+"ste.fd",ster, nstep, basept);
    grid2txt(to_string(molID)+"ele.fd",elec,nstep, basept);
    /*
     * grid2txt("ste.oe.fd",ster,true);
grid2txt("ele.oe.fd",elec,true);
*/
  return true;
}

void mpredict( int nFpt, string predfn, double intcpt, float *beta,
              float *ster,  float *elec, int *nstep, float *basept ) {
// **********************************************************************
    int nMol = 0;
    int nOKMol = 0;
    OEGraphMol molC;
    double valBio;
    oemolistream ifs;

    double sumRes = 0.0;
    double sumSsq = 0.0;
    ofstream fResid;
    fResid.open("predx.txt");

    if (!ifs.open(predfn))
        OEThrow.Fatal("Unable to read from '" + predfn + "'");
    while (OEReadMolecule(ifs, molC)) {
        ++nMol;

        memset(ster, 0, nFpt * sizeof(float));
        memset(elec, 0, nFpt * sizeof(float));
        if (!calcField(molC, ster, elec, nstep, basept, nMol)) {
            cout << "Field calculation failed for structure " << nMol << endl;
            fResid << nMol << endl;
            continue;
        }

        float predB = (float) intcpt;
        int npt = nFpt/2;
        for (int i = 0; i < nFpt; i++)
            predB += beta[i] * (i < npt ? ster[i] : elec[i-npt]);

        if (OEHasSDData(molC, "LOGBIO")) {
            OEStringToNumber(OEGetSDData(molC, "LOGBIO"), valBio);
            float resid = predB - (float) valBio;
            sumRes += resid;
            sumSsq += resid * resid;
            nOKMol++;
            fResid << nMol << ',' << predB << ',' << valBio << ',' << resid << endl;
        }
        else fResid << nMol << ',' << predB << endl;
    }
    fResid.close();
}

void fpredict( string modelfn, string predfn ){
// ***************************************************
    ifstream mdl;
    mdl.open(modelfn);
    if (!mdl) {
        cout << "Cannot open model file '" << modelfn << "' for reading\n";
        return;
    }
    string stats;
    getline(mdl,stats);
    float lo[3];
    int npts[3];

    for (int i = 0; i < 3; i++) mdl >> npts[i];
    for (int i = 0; i < 3; i++) mdl >> lo[i];

    int nFpts = npts[0] * npts[1] * npts[2];
    float *ster = new float[ nFpts ];
    float *elec = new float[ nFpts ];

    double intcpt;
    mdl >> intcpt;
    float *beta = new float[nFpts];
    for (int i = 0; i < 2*nFpts; i++) mdl >> beta[i];

    mpredict(nFpts, predfn, intcpt, beta, ster, elec, npts, lo);

}

int main(int argc, char** argv) {
// ******************************************************

    OEInterface itf(InterfaceData, argc, argv);
    rbAttn = (float) itf.Get<double>("-rbatt");
    bool mkmodel = itf.Get<bool>("-mkmodel");
    bool predict = itf.Get<bool>("-predict");
    string modelfn = itf.Get<std::string>("-model");
    string predfn = itf.Get<std::string>("-pred");
    rbAttn = (float) itf.Get<double>("-rbatt");
    int ncomp = itf.Get<int>("-ncomp");

    if (!mkmodel && !predict) {
    cout << "Nothing to do: Neither making a model nor predicting requested\n";
    return 0;
}

if (mkmodel) {

    // so we are to make a new model
    oemolistream ifs;
//    string candfn = itf.Get<std::string>("-in");
    string candfn = "in.sdf";
    if (!ifs.open(candfn))
        OEThrow.Fatal("Unable to read (tcfa-produced) structures from '" + candfn + "'");

    bool nuModel = true;
    unsigned int numMol = 0;
    OEGraphMol mol;
    float lo[3];
    int npts[3];
    if (nuModel) {
        float ctr[3], ext[3];
        float box[6] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MIN, FLT_MIN, FLT_MIN};
        while (OEReadMolecule(ifs, mol)) {
            ++numMol;
            OEGetCenterAndExtents(mol, ctr, ext);
            box[0] = min(box[0], ctr[0] - ext[0]);
            box[1] = min(box[1], ctr[1] - ext[1]);
            box[2] = min(box[2], ctr[2] - ext[2]);
            box[3] = max(box[3], ctr[0] + ext[0]);
            box[4] = max(box[4], ctr[1] + ext[1]);
            box[5] = max(box[5], ctr[2] + ext[2]);
        }
        ifs.close();

        for (unsigned int i = 0; i < 3; ++i) {
            box[i] -= 2.0f;
            box[i + 3] += 2.0f;
            npts[i] = (int) ((box[i + 3] - box[i]) / 2.0f);
            lo[i] = box[i];
        }
    }

    vector<float *> fields;
    int npt = npts[0] * npts[1] * npts[2];
    int nFpt = 2*npt;
    float *ster = new float[npt];
    float *elec = new float[npt];
    float *sumFv = new float[nFpt];
    float *ssqFv = new float[nFpt];
    float *hiFv = new float[nFpt];
    float *loFv = new float[nFpt];
    for (int i = 0; i < nFpt; ++i) {
        sumFv[i] = ssqFv[i] = 0.0f;
        hiFv[i] = -1.0e05f;
        loFv[i] = 1.0e10f;
    }
    vector<float> bioVs;
    float valBio;

    bool *OKmol = new bool[numMol];
    for (unsigned int i = 0; i < numMol; ++i) OKmol[i] = true;

    if (!ifs.open(candfn))
        OEThrow.Fatal("Unable to read from '" + candfn + "'");
    int nMol = 0;
    int nOKMol = 0;
    OEGraphMol molC;
    while (OEReadMolecule(ifs, molC)) {
        ++nMol;

        if (!OEHasSDData(molC, "LOGBIO")) {
            cout << "No LOGBIO value for structure " << nMol << endl;
            OKmol[nMol] = false;
            continue;
        }

        memset(ster, 0, npt * sizeof(float));
        memset(elec, 0, npt * sizeof(float));
        if (!calcField(molC, ster, elec, npts, lo, nMol)) {
            cout << "Field calculation failed for structure " << nMol << endl;
            OKmol[nMol] = false;
            continue;
        }
        ++nOKMol;

        OEStringToNumber(OEGetSDData(molC, "LOGBIO"), valBio);
        bioVs.push_back(valBio);

        float *nowFld = new float[nFpt];
        memcpy(nowFld, ster, sizeof(float) * npt);
        memcpy(nowFld + npt, elec, sizeof(float) * npt);
        fields.push_back(nowFld);
        for (int i = 0; i < nFpt; i++) {
            sumFv[i] += nowFld[i];
            ssqFv[i] += nowFld[i] * nowFld[i];;
            if (nowFld[i] < loFv[i]) loFv[i] = nowFld[i];
            if (nowFld[i] > hiFv[i]) hiFv[i] = nowFld[i];
        }
    }
    wcsvF("sum.fd", sumFv, nFpt, 1);
    wcsvF("ssq.fd", ssqFv, nFpt, 1);
    cout << "Fields done\n";
// identifying columns where no field variance was experienced
    bool *OKfVal = new bool[nFpt];
    for (int i = 0; i < nFpt; ++i) OKfVal[i] = false;

    int nOKfVal = 0;
    for (int i = 0; i < nFpt; ++i)
        if (ssqFv[i] > 0.0000001 && hiFv[i] != loFv[i]) {
            ++nOKfVal;
            OKfVal[i] = true;
        }
// dropping missing rows (structures) ..
    double *ys = new double[nOKMol];
    int iy = 0;
    for (int i = 0; i < nMol; i++)
        if (OKmol[i]) {
            ys[iy] = (double) bioVs[iy];
            ++iy;
        }
    double *xs = new double[nOKMol * nOKfVal];
    double *xsig = new double[nOKfVal + 2];
    double *xmean = new double[nOKfVal + 2];

    int nster = 0;
    int ifd = 0;
    int ix = 0;
    for (int i = 0; i < nMol; i++)
        if (OKmol[i]) {
            float *Vals = fields[ix];
// .. dropping signal-less lattice points
            for (int j = 0; j < nFpt; j++)
                if (OKfVal[j]) {
                    xs[ifd] = (double) Vals[j];
                    // and totaling by field type
                    xmean[nOKfVal + (j < nFpt / 2 ? 0 : 1)] += xs[ifd];
                    xsig[nOKfVal + (j < nFpt / 2 ? 0 : 1)] += xs[ifd] * xs[ifd];
                    if (i == 0 && j < nFpt / 2) nster++;
                    ++ifd;
                }
            ++ix;
        }
    // CoMFA field weights are by the field, not the lattice point
    xsig[nOKfVal] = ((xsig[nOKfVal] - ((xmean[nOKfVal] * xmean[nOKfVal]) / (double) (nster * nOKMol))) /
                     ((double) ((nster * nOKMol) - 1)));
    xsig[nOKfVal] = sqrt(xsig[nOKfVal] * ((double) nster));
    xsig[nOKfVal + 1] = (
            (xsig[nOKfVal + 1] - ((xmean[nOKfVal + 1] * xmean[nOKfVal + 1]) / (double) ((nOKfVal - nster) * nOKMol))) /
            ((double) ((nOKfVal - nster) * nOKMol - 1)));
    xsig[nOKfVal + 1] = sqrt(xsig[nOKfVal + 1] * ((double) (nOKfVal - nster)));
    xmean[nOKfVal] /= (double) (nOKMol * nster);
    xmean[nOKfVal + 1] /= (double) (nOKMol * (nOKfVal - nster));
    for (int j = 0; j < nOKfVal; j++) {
        xsig[j] = j < nster ? xsig[nOKfVal] : xsig[nOKfVal + 1];
        xmean[j] = j < nster ? xmean[nOKfVal] : xmean[nOKfVal + 1];
    }
    for (int i = 0; i < nOKMol; i++)
        for (int j = 0; j < nOKfVal; j++)
            xs[i * nOKfVal + j] = (xs[i * nOKfVal + j] - xmean[j]) / xsig[j];
    wcsv("xavg.csv", xmean, nOKfVal, 1);
    wcsv("xsdv.csv", xsig, nOKfVal, 1);
    wcsv("x.csv", xs, nOKMol, nOKfVal);
    wcsv("y.csv", ys, nOKMol, 1);
    cout << "Field values filtered and scaled\n";
    double *ywts = new double[nOKMol];
    int ierr;
    int nys = 1;
    double *resdl = new double[2 * nOKMol];
    double q2;
    double qsdep;
    int mxncomp = ncomp != -1 ? ncomp :
                  sampls(20, nOKMol, ys, nOKfVal, xs, &qsdep, &q2, resdl);
//    int mxncomp = 3;
    int maxiter = 100;
    double eps = 1.0e-4;
    int iboot = 0;
    int icros = 0;
    int icent = 2;
    double *weytb = new double[nOKMol];
    double *xbar = new double[nOKfVal + 1];
    double *xscal = new double[nOKfVal + 1];
    double *intcpt = new double[3];
    double *coeff = new double[2 * nOKfVal];
    double *varnc = new double[mxncomp * nOKfVal];
    double *wtx = new double[mxncomp * nOKfVal];
    double *wty = new double[mxncomp];
    double *loadings = new double[mxncomp * nOKfVal];
    double *latentx = new double[mxncomp * nOKMol];
    double *latenty = new double[mxncomp * nOKMol];
    double *inner = new double[mxncomp];
    double *ssqRsdl = new double[2];
    double *r2 = new double[2];
    double *sdep = new double[nOKfVal * mxncomp];
//   double *q2 = new double[2*mxncomp];
    int *optncomp = new int[2];
    double *ypred = new double[nOKMol];
    double *sss = new double[2];
    double *ssy = new double[nOKMol];
    int *iscratch = new int[nOKMol];
    double *scratch = new double[nOKMol * (nOKfVal + 2) + (nOKfVal + 1)];
    for (iy = 0; iy < nOKMol * (nOKfVal + 2) + (nOKfVal + 1); iy++) scratch[iy] = 0.0;
    for (iy = 0; iy < 3; ++iy) intcpt[iy] = 0.0;
    for (iy = 0; iy < 2; ++iy) {
        optncomp[iy] = 0;
        r2[iy] = sss[iy] = ssqRsdl[iy] = 0.0;
    }
    for (iy = 0; iy < nOKMol; ++iy) {
        iscratch[iy] = 0;
        ssy[iy] = weytb[iy] = ypred[iy] = 0.0;
        ywts[iy] = 1.0;
    }
    for (iy = 0; iy < 2 * nOKMol; ++iy) resdl[iy] = 0.0;
    for (iy = 0; iy < 2 * nOKfVal; ++iy) coeff[iy] = 0.0;
    for (iy = 0; iy < mxncomp * nOKfVal; ++iy)
        varnc[iy] = wtx[iy] = loadings[iy] = sdep[iy] = 0.0;
    for (iy = 0; iy < mxncomp; ++iy) wty[iy] = inner[iy] = 0.0;
    for (iy = 0; iy < mxncomp * nOKMol; ++iy) latentx[iy] = latenty[iy] = 0.0;
    int plsOK = plsjer(&nOKMol, &nOKfVal, &nys, &mxncomp, &maxiter, &eps, &iboot,
                       &icros, &icent, xs, ys, ywts, weytb, xbar, xscal, intcpt, coeff,
                       varnc, wtx, wty, loadings, latentx, latenty, inner, ssqRsdl, r2,
                       NULL, NULL, optncomp, ypred, resdl, sss, ssy, &ierr, scratch, iscratch);
    if (plsOK != 1) {
        cout << "PLS error: ID =" << ierr << endl;
        return 1;
    }
    for (iy = 0; iy < nOKfVal; iy++) coeff[iy] /= xsig[iy];
    for (iy = 0; iy < nOKfVal; iy++) intcpt[0] -= xmean[iy] * coeff[iy];

// expand model; write out the contours
    float *beta = new float[nFpt];
    int nf = 0;
    for (int j = 0; j < nFpt; j++)
        if (OKfVal[j]) {
            beta[j] = (float) coeff[nf];
            nf++;
        }
        else beta[j] = 0.0;

    float *cntrSter = new float[npt];
    float *cntrElec = new float[npt];
    nf = 0;
    for (int j = 0; j < nFpt; j++)
        if (OKfVal[j]) {
            if (j < npt)
                cntrSter[j] = (float) (coeff[nf] * xsig[nf]);
            else
                cntrElec[j - npt] = (float) (coeff[nf] * xsig[nf]);
            nf++;
        } else {
            cntrSter[j] = 0.0;
            cntrElec[j] = 0.0;
        }

    savemodel(modelfn, nFpt, intcpt[0], beta, npts, lo, mxncomp, qsdep, q2, r2[0], resdl[0]);
    grid2txt("stdev*coeff.ste", cntrSter, npts, lo);
    grid2txt("stdev*coeff.ele", cntrElec, npts, lo);

    if (predict)
        mpredict(nFpt, predfn, intcpt[0], beta, ster, elec, npts, lo);
}

    if (predict && !mkmodel)
        fpredict(modelfn, predfn);

    return 0;
}
