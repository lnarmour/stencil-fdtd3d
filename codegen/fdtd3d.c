// This file is generated from test alphabets program by code generator in alphaz
// To compile this code, use -lm option for math library.

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "external_functions.h"

// Common Macros
#define max(x, y)   ((x)>(y) ? (x) : (y))
#define MAX(x, y)	((x)>(y) ? (x) : (y))
#define min(x, y)   ((x)>(y) ? (y) : (x))
#define MIN(x, y)	((x)>(y) ? (y) : (x))
#define CEILD(n,d)  (int)ceil(((double)(n))/((double)(d)))
#define ceild(n,d)  (int)ceil(((double)(n))/((double)(d)))
#define FLOORD(n,d) (int)floor(((double)(n))/((double)(d)))
#define floord(n,d) (int)floor(((double)(n))/((double)(d)))
#define CDIV(x,y)    CEILD((x),(y))
#define div(x,y)    CDIV((x),(y))
#define FDIV(x,y)    FLOORD((x),(y))
#define LB_SHIFT(b,s)  ((int)ceild(b,s) * s)
#define MOD(i,j)   ((i)%(j))
#define mallocCheck(v,s,d) if ((v) == NULL) { printf("Failed to allocate memory for %s : size=%lu\n", "sizeof(d)*(s)", sizeof(d)*(s)); exit(-1); }
// Reduction Operators
#define RADD(x,y)    ((x)+=(y))
#define RMUL(x,y)    ((x)*=(y))
#define RMAX(x,y)    ((x)=MAX((x),(y)))
#define RMIN(x,y)    ((x)=MIN((x),(y)))

// Common functions for min and max
//functions for integer max
inline int __max_int(int x, int y){
	return ((x)>(y) ? (x) : (y));
}

inline short __max_short(short x, short y){
	return ((x)>(y) ? (x) : (y));
}

inline long __max_long(long x, long y){
	return ((x)>(y) ? (x) : (y));
}

inline unsigned int __max_unsigned_int(unsigned int x, unsigned int y){
	return ((x)>(y) ? (x) : (y));
}

inline unsigned short __max_unsigned_short(unsigned short x, unsigned short y){
	return ((x)>(y) ? (x) : (y));
}

//function for float max
inline float __max_float(float x, float y){
	return ((x)>(y) ? (x) : (y));
}

//function for double max
inline double __max_double(double x, double y){
	return ((x)>(y) ? (x) : (y));
}

//function for integer min
inline int __min_int(int x, int y){
	return ((x)>(y) ? (y) : (x));
}

inline short __min_short(short x, short y){
	return ((x)>(y) ? (y) : (x));
}

inline long __min_long(long x, long y){
	return ((x)>(y) ? (y) : (x));
}

inline unsigned int __min_unsigned_int(unsigned int x, unsigned int y){
	return ((x)>(y) ? (y) : (x));
}

inline unsigned short __min_unsigned_short(unsigned short x, unsigned short y){
	return ((x)>(y) ? (y) : (x));
}

inline unsigned long __min_unsigned_long(unsigned long x, unsigned long y){
	return ((x)>(y) ? (y) : (x));
}

inline float __min_float(float x, float y){
	return ((x)>(y) ? (y) : (x));
}

inline double __min_double(double x, double y){
	return ((x)>(y) ? (y) : (x));
}






//Memory Macros
#define ex(i,j,k) ex[i][j][k]
#define ey(i,j,k) ey[i][j][k]
#define ez(i,j,k) ez[i][j][k]
#define Ex(t,i,j,k) Ex[MOD(t,2)][i][j][k]
#define Ey(t,i,j,k) Ey[MOD(t,2)][i][j][k]
#define Hx(t,i,j,k) Hx[MOD(t,2)][i][j][k]
#define Hy(t,i,j,k) Hy[MOD(t,2)][i][j][k]
#define Hz(t,i,j,k) Hz[MOD(t,2)][i][j][k]
#define Cexhy(i) Cexhy[i]
#define Cexhz(i) Cexhz[i]
#define Ceyhx(i) Ceyhx[i]
#define Ceyhz(i) Ceyhz[i]
#define Cezhx(i) Cezhx[i]
#define Cezhy(i) Cezhy[i]
#define S(t,i,j) S[t][i][j]
#define Ez(t,i,j,k) Ez[MOD(t,2)][i][j][k]

void fdtd3d(long nt, long Nx, long Ny, long Nz, long Sx, long Sy, long Sz, double* f0, double* Lf, double* CF, double* er_x, double* er_y, double* er_z, double* Lx, double* Ly, double* Lz, double* e0, double* u0, double* c0, double* dx, double* dy, double* dz, double* dt, double**** Ez){
	///Parameter checking
	if (!((Nx >= 10 && Ny >= 10 && Nz >= 10 && nt >= 1))) {
		printf("The value of parameters are not valid.\n");
		exit(-1);
	}
	//Memory Allocation
	int mz1, mz2, mz3, mz4;
	
	double* _lin_ex = (double*)malloc(sizeof(double)*((Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_ex, ((Nx+1) * (Ny+1) * (Nz+1)), double);
	double*** ex = (double***)malloc(sizeof(double**)*(Nx+1));
	mallocCheck(ex, (Nx+1), double**);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		ex[mz1] = (double**)malloc(sizeof(double*)*(Ny+1));
		mallocCheck(ex[mz1], (Ny+1), double*);
		for (mz2=0;mz2 < Ny+1; mz2++) {
			ex[mz1][mz2] = &_lin_ex[(mz1*((Ny+1) * (Nz+1))) + (mz2*(Nz+1))];
		}
	}
	
	double* _lin_ey = (double*)malloc(sizeof(double)*((Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_ey, ((Nx+1) * (Ny+1) * (Nz+1)), double);
	double*** ey = (double***)malloc(sizeof(double**)*(Nx+1));
	mallocCheck(ey, (Nx+1), double**);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		ey[mz1] = (double**)malloc(sizeof(double*)*(Ny+1));
		mallocCheck(ey[mz1], (Ny+1), double*);
		for (mz2=0;mz2 < Ny+1; mz2++) {
			ey[mz1][mz2] = &_lin_ey[(mz1*((Ny+1) * (Nz+1))) + (mz2*(Nz+1))];
		}
	}
	
	double* _lin_ez = (double*)malloc(sizeof(double)*((Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_ez, ((Nx+1) * (Ny+1) * (Nz+1)), double);
	double*** ez = (double***)malloc(sizeof(double**)*(Nx+1));
	mallocCheck(ez, (Nx+1), double**);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		ez[mz1] = (double**)malloc(sizeof(double*)*(Ny+1));
		mallocCheck(ez[mz1], (Ny+1), double*);
		for (mz2=0;mz2 < Ny+1; mz2++) {
			ez[mz1][mz2] = &_lin_ez[(mz1*((Ny+1) * (Nz+1))) + (mz2*(Nz+1))];
		}
	}
	
	double* _lin_Ex = (double*)malloc(sizeof(double)*(2 * (Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_Ex, (2 * (Nx+1) * (Ny+1) * (Nz+1)), double);
	double**** Ex = (double****)malloc(sizeof(double***)*(2));
	mallocCheck(Ex, (2), double***);
	for (mz1=0;mz1 < 2; mz1++) {
		Ex[mz1] = (double***)malloc(sizeof(double**)*(Nx+1));
		mallocCheck(Ex[mz1], (Nx+1), double**);
		for (mz2=0;mz2 < Nx+1; mz2++) {
			Ex[mz1][mz2] = (double**)malloc(sizeof(double*)*(Ny+1));
			mallocCheck(Ex[mz1][mz2], (Ny+1), double*);
			for (mz3=0;mz3 < Ny+1; mz3++) {
				Ex[mz1][mz2][mz3] = &_lin_Ex[(mz1*((Nx+1) * (Ny+1) * (Nz+1))) + (mz2*((Ny+1) * (Nz+1))) + (mz3*(Nz+1))];
			}
		}
	}
	
	double* _lin_Ey = (double*)malloc(sizeof(double)*(2 * (Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_Ey, (2 * (Nx+1) * (Ny+1) * (Nz+1)), double);
	double**** Ey = (double****)malloc(sizeof(double***)*(2));
	mallocCheck(Ey, (2), double***);
	for (mz1=0;mz1 < 2; mz1++) {
		Ey[mz1] = (double***)malloc(sizeof(double**)*(Nx+1));
		mallocCheck(Ey[mz1], (Nx+1), double**);
		for (mz2=0;mz2 < Nx+1; mz2++) {
			Ey[mz1][mz2] = (double**)malloc(sizeof(double*)*(Ny+1));
			mallocCheck(Ey[mz1][mz2], (Ny+1), double*);
			for (mz3=0;mz3 < Ny+1; mz3++) {
				Ey[mz1][mz2][mz3] = &_lin_Ey[(mz1*((Nx+1) * (Ny+1) * (Nz+1))) + (mz2*((Ny+1) * (Nz+1))) + (mz3*(Nz+1))];
			}
		}
	}
	
	double* _lin_Hx = (double*)malloc(sizeof(double)*(2 * (Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_Hx, (2 * (Nx+1) * (Ny+1) * (Nz+1)), double);
	double**** Hx = (double****)malloc(sizeof(double***)*(2));
	mallocCheck(Hx, (2), double***);
	for (mz1=0;mz1 < 2; mz1++) {
		Hx[mz1] = (double***)malloc(sizeof(double**)*(Nx+1));
		mallocCheck(Hx[mz1], (Nx+1), double**);
		for (mz2=0;mz2 < Nx+1; mz2++) {
			Hx[mz1][mz2] = (double**)malloc(sizeof(double*)*(Ny+1));
			mallocCheck(Hx[mz1][mz2], (Ny+1), double*);
			for (mz3=0;mz3 < Ny+1; mz3++) {
				Hx[mz1][mz2][mz3] = &_lin_Hx[(mz1*((Nx+1) * (Ny+1) * (Nz+1))) + (mz2*((Ny+1) * (Nz+1))) + (mz3*(Nz+1))];
			}
		}
	}
	
	double* _lin_Hy = (double*)malloc(sizeof(double)*(2 * (Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_Hy, (2 * (Nx+1) * (Ny+1) * (Nz+1)), double);
	double**** Hy = (double****)malloc(sizeof(double***)*(2));
	mallocCheck(Hy, (2), double***);
	for (mz1=0;mz1 < 2; mz1++) {
		Hy[mz1] = (double***)malloc(sizeof(double**)*(Nx+1));
		mallocCheck(Hy[mz1], (Nx+1), double**);
		for (mz2=0;mz2 < Nx+1; mz2++) {
			Hy[mz1][mz2] = (double**)malloc(sizeof(double*)*(Ny+1));
			mallocCheck(Hy[mz1][mz2], (Ny+1), double*);
			for (mz3=0;mz3 < Ny+1; mz3++) {
				Hy[mz1][mz2][mz3] = &_lin_Hy[(mz1*((Nx+1) * (Ny+1) * (Nz+1))) + (mz2*((Ny+1) * (Nz+1))) + (mz3*(Nz+1))];
			}
		}
	}
	
	double* _lin_Hz = (double*)malloc(sizeof(double)*(2 * (Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_Hz, (2 * (Nx+1) * (Ny+1) * (Nz+1)), double);
	double**** Hz = (double****)malloc(sizeof(double***)*(2));
	mallocCheck(Hz, (2), double***);
	for (mz1=0;mz1 < 2; mz1++) {
		Hz[mz1] = (double***)malloc(sizeof(double**)*(Nx+1));
		mallocCheck(Hz[mz1], (Nx+1), double**);
		for (mz2=0;mz2 < Nx+1; mz2++) {
			Hz[mz1][mz2] = (double**)malloc(sizeof(double*)*(Ny+1));
			mallocCheck(Hz[mz1][mz2], (Ny+1), double*);
			for (mz3=0;mz3 < Ny+1; mz3++) {
				Hz[mz1][mz2][mz3] = &_lin_Hz[(mz1*((Nx+1) * (Ny+1) * (Nz+1))) + (mz2*((Ny+1) * (Nz+1))) + (mz3*(Nz+1))];
			}
		}
	}
	
	double Chxey;
	double Chxez;
	double* Cexhy = (double*)malloc(sizeof(double)*(Ny+1));
	mallocCheck(Cexhy, (Ny+1), double);
	
	double* Cexhz = (double*)malloc(sizeof(double)*(Ny+1));
	mallocCheck(Cexhz, (Ny+1), double);
	
	double Chyex;
	double Chyez;
	double* Ceyhx = (double*)malloc(sizeof(double)*(Ny+1));
	mallocCheck(Ceyhx, (Ny+1), double);
	
	double* Ceyhz = (double*)malloc(sizeof(double)*(Ny+1));
	mallocCheck(Ceyhz, (Ny+1), double);
	
	double Chzex;
	double Chzey;
	double* Cezhx = (double*)malloc(sizeof(double)*(Ny+1));
	mallocCheck(Cezhx, (Ny+1), double);
	
	double* Cezhy = (double*)malloc(sizeof(double)*(Ny+1));
	mallocCheck(Cezhy, (Ny+1), double);
	
	double* _lin_S = (double*)malloc(sizeof(double)*((Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_S, ((Nx+1) * (Ny+1) * (Nz+1)), double);
	double*** S = (double***)malloc(sizeof(double**)*(Nx+1));
	mallocCheck(S, (Nx+1), double**);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		S[mz1] = (double**)malloc(sizeof(double*)*(Ny+1));
		mallocCheck(S[mz1], (Ny+1), double*);
		for (mz2=0;mz2 < Ny+1; mz2++) {
			S[mz1][mz2] = &_lin_S[(mz1*((Ny+1) * (Nz+1))) + (mz2*(Nz+1))];
		}
	}
	
	double V;
	#define S0(i,j,k,i3,i4,i5) ex(i3,i4,i5) = *er_x
	#define S1(i,j,k,i3,i4,i5) ex(i3,i4,i5) = 1
	#define S2(i,j,k,i3,i4,i5) ey(i3,i4,i5) = *er_y
	#define S3(i,j,k,i3,i4,i5) ey(i3,i4,i5) = 1
	#define S4(i,j,k,i3,i4,i5) ez(i3,i4,i5) = *er_z
	#define S5(i,j,k,i3,i4,i5) ez(i3,i4,i5) = 1
	#define S6(i0,i1,i2,i3,i4,i5) Chxey = (*dt)/((*u0)*(*dz))
	#define S7(i0,i1,i2,i3,i4,i5) Chxez = (*dt)/((*u0)*(*dy))
	#define S8(i,j,k,i3,i4,i5) Cexhy(i4) = (*dt)/(((*e0)*(ex(i3,i4,i5)))*(*dz))
	#define S9(i,j,k,i3,i4,i5) Cexhz(i4) = (*dt)/(((*e0)*(ex(i3,i4,i5)))*(*dy))
	#define S10(i0,i1,i2,i3,i4,i5) Chyex = (*dt)/((*u0)*(*dz))
	#define S11(i0,i1,i2,i3,i4,i5) Chyez = (*dt)/((*u0)*(*dx))
	#define S12(i,j,k,i3,i4,i5) Ceyhx(i4) = (*dt)/(((*e0)*(ey(i3,i4,i5)))*(*dz))
	#define S13(i,j,k,i3,i4,i5) Ceyhz(i4) = (*dt)/(((*e0)*(ey(i3,i4,i5)))*(*dx))
	#define S14(i0,i1,i2,i3,i4,i5) Chzex = (*dt)/((*u0)*(*dy))
	#define S15(i0,i1,i2,i3,i4,i5) Chzey = (*dt)/((*u0)*(*dx))
	#define S16(i,j,k,i3,i4,i5) Cezhx(i4) = (*dt)/(((*e0)*(ez(i3,i4,i5)))*(*dy))
	#define S17(i,j,k,i3,i4,i5) Cezhy(i4) = (*dt)/(((*e0)*(ez(i3,i4,i5)))*(*dx))
	#define S18(t,i,j,k,i4,i5) S(k,i4,i5) = RampSig(i,*f0,*dt)
	#define S19(t,i,j,k,i4,i5) S(k,i4,i5) = 0
	#define S20(t,i,j,k,i4,i5) Hx(i,k,i4,i5) = 0
	#define S21(t,i,j,k,i4,i5) Hx(i,k,i4,i5) = Hx(i-1,k,i4,i5)
	#define S22(t,i,j,k,i4,i5) Hx(i,k,i4,i5) = ((Hx(i-1,k,i4,i5))+((Chxey)*((Ey(i-1,k,i4,i5+1))-(Ey(i-1,k,i4,i5)))))-((Chxez)*((Ez(i-1,k,i4+1,i5))-(Ez(i-1,k,i4,i5))))
	#define S23(t,i,j,k,i4,i5) Hy(i,k,i4,i5) = 0
	#define S24(t,i,j,k,i4,i5) Hy(i,k,i4,i5) = Hy(i-1,k,i4,i5)
	#define S25(t,i,j,k,i4,i5) Hy(i,k,i4,i5) = ((Hy(i-1,k,i4,i5))+((Chyez)*((Ez(i-1,k+1,i4,i5))-(Ez(i-1,k,i4,i5)))))-((Chyex)*((Ex(i-1,k,i4,i5+1))-(Ex(i-1,k,i4,i5))))
	#define S26(t,i,j,k,i4,i5) Hz(i,k,i4,i5) = 0
	#define S27(t,i,j,k,i4,i5) Hz(i,k,i4,i5) = Hz(i-1,k,i4,i5)
	#define S28(t,i,j,k,i4,i5) Hz(i,k,i4,i5) = ((Hz(i-1,k,i4,i5))+((Chzex)*((Ex(i-1,k,i4+1,i5))-(Ex(i-1,k,i4,i5)))))-((Chzey)*((Ey(i-1,k+1,i4,i5))-(Ey(i-1,k,i4,i5))))
	#define S29(t,i,j,k,i4,i5) Ex(i,k,i4,i5) = 0
	#define S30(t,i,j,k,i4,i5) Ex(i,k,i4,i5) = Ex(i-1,k,i4,i5)
	#define S31(t,i,j,k,i4,i5) Ex(i,k,i4,i5) = ((Ex(i-1,k,i4,i5))+((Cexhz(i4))*((Hz(i,k,i4,i5))-(Hz(i,k,i4-1,i5)))))-((Cexhy(i4))*((Hy(i,k,i4,i5))-(Hy(i,k,i4,i5-1))))
	#define S32(t,i,j,k,i4,i5) Ey(i,k,i4,i5) = 0
	#define S33(t,i,j,k,i4,i5) Ey(i,k,i4,i5) = Ey(i-1,k,i4,i5)
	#define S34(t,i,j,k,i4,i5) Ey(i,k,i4,i5) = ((Ey(i-1,k,i4,i5))+((Ceyhx(i4))*((Hx(i,k,i4,i5))-(Hx(i,k,i4,i5-1)))))-((Ceyhz(i4))*((Hz(i,k,i4,i5))-(Hz(i,k-1,i4,i5))))
	#define S35(t,i,j,k,i4,i5) Ez(i,k,i4,i5) = 0
	#define S36(t,i,j,k,i4,i5) Ez(i,k,i4,i5) = Ez(i-1,k,i4,i5)
	#define S37(t,i,j,k,i4,i5) Ez(i,k,i4,i5) = (((S(k,i4,i5))+(Ez(i-1,k,i4,i5)))+((Cezhy(i4))*((Hy(i,k,i4,i5))-(Hy(i,k-1,i4,i5)))))-((Cezhx(i4))*((Hx(i,k,i4,i5))-(Hx(i,k,i4-1,i5))))
	#define S38(t,i,j,k,i4,i5) V = PrintMe(Ez(i,k,i4,i5))
	{
		//Domain
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==1 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && 5i4>=4Ny+1 && i4>=1 && Nz>=i5 && i5>=1 && Ny>=i4 && i3>=1 && Nx>=i3}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==1 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && 4Ny>=5i4 && i4>=1 && Nz>=i5 && i5>=1 && Ny>=i4 && i3>=1 && Nx>=i3}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==1 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && 5i4>=4Ny+1 && i4>=1 && Nz>=i5 && i5>=1 && Ny>=i4 && i3>=1 && Nx>=i3}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==1 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && 4Ny>=5i4 && i4>=1 && Nz>=i5 && i5>=1 && Ny>=i4 && i3>=1 && Nx>=i3}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==1 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && 5i4>=4Ny+1 && i4>=1 && Nz>=i5 && i5>=1 && Ny>=i4 && i3>=1 && Nx>=i3}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==1 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && 4Ny>=5i4 && i4>=1 && Nz>=i5 && i5>=1 && Ny>=i4 && i3>=1 && Nx>=i3}
		//{i0,i1,i2,i3,i4,i5|i5==0 && i4==0 && i3==0 && i2==0 && i1==0 && i0==2 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1}
		//{i0,i1,i2,i3,i4,i5|i5==0 && i4==0 && i3==0 && i2==0 && i1==0 && i0==2 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==3 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1 && Nx>=i3 && i3>=1 && i4>=1 && Ny>=i4 && Nz>=i5 && i5>=1}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==3 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1 && Nx>=i3 && i3>=1 && i4>=1 && Ny>=i4 && Nz>=i5 && i5>=1}
		//{i0,i1,i2,i3,i4,i5|i5==0 && i4==0 && i3==0 && i2==0 && i1==0 && i0==2 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1}
		//{i0,i1,i2,i3,i4,i5|i5==0 && i4==0 && i3==0 && i2==0 && i1==0 && i0==2 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==4 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1 && Nx>=i3 && i3>=1 && i4>=1 && Ny>=i4 && Nz>=i5 && i5>=1}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==4 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1 && Nx>=i3 && i3>=1 && i4>=1 && Ny>=i4 && Nz>=i5 && i5>=1}
		//{i0,i1,i2,i3,i4,i5|i5==0 && i4==0 && i3==0 && i2==0 && i1==0 && i0==2 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1}
		//{i0,i1,i2,i3,i4,i5|i5==0 && i4==0 && i3==0 && i2==0 && i1==0 && i0==2 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==5 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1 && Nx>=i3 && i3>=1 && i4>=1 && Ny>=i4 && Nz>=i5 && i5>=1}
		//{i,j,k,i3,i4,i5|k==0 && j==0 && i==5 && Nx>=10 && Ny>=10 && Nz>=10 && nt>=1 && Nx>=i3 && i3>=1 && i4>=1 && Ny>=i4 && Nz>=i5 && i5>=1}
		//{t,i,j,k,i4,i5|i5==Sz && i4==Sy && k==Sx && j==1 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && Ny>=Sy && Nx>=Sx && Sy>=1 && Sz>=1 && Nz>=Sz && i>=0 && nt>=i+1 && Sx>=1}
		//{t,i,j,k,i4,i5|j==1 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=0 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4 && i5>=Sz+1 && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|j==1 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=0 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4 && i5>=1 && Sz>=i5+1 && Nz>=i5} || {t,i,j,k,i4,i5|i5==Sz && j==1 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && Sz>=1 && Nz>=Sz && i>=0 && nt>=i+1 && k>=1 && Nx>=k && i4>=Sy+1 && i4>=1 && Ny>=i4} || {t,i,j,k,i4,i5|i5==Sz && j==1 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && Sz>=1 && Nz>=Sz && i>=0 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Sy>=i4+1 && Ny>=i4} || {t,i,j,k,i4,i5|i5==Sz && i4==Sy && j==1 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && Sy>=1 && Ny>=Sy && Sz>=1 && Nz>=Sz && i>=0 && nt>=i+1 && k>=Sx+1 && k>=1 && Nx>=k} || {t,i,j,k,i4,i5|i5==Sz && i4==Sy && j==1 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && Sy>=1 && Ny>=Sy && Sz>=1 && Nz>=Sz && i>=0 && nt>=i+1 && k>=1 && Sx>=k+1 && Nx>=k}
		//{t,i,j,k,i4,i5|j==0 && i==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && Nz>=i5 && i5>=1 && Ny>=i4 && i4>=1 && k>=1 && Nx>=k}
		//{t,i,j,k,i4,i5|i5==1 && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4} || {t,i,j,k,i4,i5|i4==1 && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|i4==Ny && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|i5==Nz && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4}
		//{t,i,j,k,i4,i5|j==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && i4>=2 && Ny>=i4+1 && i5>=2 && Nz>=i5+1 && Nx>=k && k>=1 && nt>=i+1}
		//{t,i,j,k,i4,i5|j==0 && i==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && Nz>=i5 && i5>=1 && Ny>=i4 && i4>=1 && k>=1 && Nx>=k}
		//{t,i,j,k,i4,i5|i5==1 && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4} || {t,i,j,k,i4,i5|k==1 && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && i4>=1 && Ny>=i4 && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|k==Nx && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && i4>=1 && Ny>=i4 && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|i5==Nz && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4}
		//{t,i,j,k,i4,i5|j==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && k>=2 && Nx>=k+1 && i5>=2 && Nz>=i5+1 && Ny>=i4 && i4>=1 && nt>=i+1}
		//{t,i,j,k,i4,i5|j==0 && i==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && Nz>=i5 && i5>=1 && Ny>=i4 && i4>=1 && k>=1 && Nx>=k}
		//{t,i,j,k,i4,i5|i4==1 && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|k==1 && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && i4>=1 && Ny>=i4 && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|k==Nx && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && i4>=1 && Ny>=i4 && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|i4==Ny && j==0 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i5>=1 && Nz>=i5}
		//{t,i,j,k,i4,i5|j==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && k>=2 && Nx>=k+1 && i4>=2 && Ny>=i4+1 && Nz>=i5 && i5>=1 && nt>=i+1}
		//{t,i,j,k,i4,i5|j==2 && i==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && Nz>=i5 && i5>=1 && Ny>=i4 && i4>=1 && k>=1 && Nx>=k}
		//{t,i,j,k,i4,i5|i5==1 && j==2 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4} || {t,i,j,k,i4,i5|i4==1 && j==2 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i5>=1 && Nz>=i5}
		//{t,i,j,k,i4,i5|j==2 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && i4>=2 && i5>=2 && Nz>=i5 && Ny>=i4 && Nx>=k && k>=1 && nt>=i+1}
		//{t,i,j,k,i4,i5|j==2 && i==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && Nz>=i5 && i5>=1 && Ny>=i4 && i4>=1 && k>=1 && Nx>=k}
		//{t,i,j,k,i4,i5|i5==1 && j==2 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4} || {t,i,j,k,i4,i5|k==1 && j==2 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && i4>=1 && Ny>=i4 && i5>=1 && Nz>=i5}
		//{t,i,j,k,i4,i5|j==2 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && k>=2 && i5>=2 && Nz>=i5 && Ny>=i4 && i4>=1 && Nx>=k && nt>=i+1}
		//{t,i,j,k,i4,i5|j==2 && i==0 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && Nz>=i5 && i5>=1 && Ny>=i4 && i4>=1 && k>=1 && Nx>=k}
		//{t,i,j,k,i4,i5|i4==1 && j==2 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && k>=1 && Nx>=k && i5>=1 && Nz>=i5} || {t,i,j,k,i4,i5|k==1 && j==2 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && nt>=i+1 && i4>=1 && Ny>=i4 && i5>=1 && Nz>=i5}
		//{t,i,j,k,i4,i5|j==2 && t==6 && nt>=1 && Nx>=10 && Ny>=10 && Nz>=10 && i>=1 && k>=2 && i4>=2 && Nz>=i5 && i5>=1 && Ny>=i4 && Nx>=k && nt>=i+1}
		//{t,i,j,k,i4,i5|j==3 && t==6 && Nx>=10 && Ny>=10 && Nz>=10 && i>=0 && nt>=i+1 && k>=1 && Nx>=k && i4>=1 && Ny>=i4 && i5>=1 && Nz>=i5 && nt>=1}
		int c2,c4,c5,c6;
		for(c4=1;c4 <= Nx;c4+=1)
		 {
		 	for(c5=1;c5 <= (4*Ny) / 5;c5+=1)
		 	 {
		 	 	for(c6=1;c6 <= Nz;c6+=1)
		 	 	 {
		 	 	 	S1((1),(0),(0),(c4),(c5),(c6));
		 	 	 	S3((1),(0),(0),(c4),(c5),(c6));
		 	 	 	S5((1),(0),(0),(c4),(c5),(c6));
		 	 	 }
		 	 }
		 	for(c5=ceild(4*Ny+1, 5);c5 <= Ny;c5+=1)
		 	 {
		 	 	for(c6=1;c6 <= Nz;c6+=1)
		 	 	 {
		 	 	 	S0((1),(0),(0),(c4),(c5),(c6));
		 	 	 	S2((1),(0),(0),(c4),(c5),(c6));
		 	 	 	S4((1),(0),(0),(c4),(c5),(c6));
		 	 	 }
		 	 }
		 }
		S6((2),(0),(0),(0),(0),(0));
		S7((2),(0),(0),(0),(0),(0));
		S10((2),(0),(0),(0),(0),(0));
		S11((2),(0),(0),(0),(0),(0));
		S14((2),(0),(0),(0),(0),(0));
		S15((2),(0),(0),(0),(0),(0));
		for(c4=1;c4 <= Nx;c4+=1)
		 {
		 	for(c5=1;c5 <= Ny;c5+=1)
		 	 {
		 	 	for(c6=1;c6 <= Nz;c6+=1)
		 	 	 {
		 	 	 	S8((3),(0),(0),(c4),(c5),(c6));
		 	 	 	S9((3),(0),(0),(c4),(c5),(c6));
		 	 	 }
		 	 }
		 }
		for(c4=1;c4 <= Nx;c4+=1)
		 {
		 	for(c5=1;c5 <= Ny;c5+=1)
		 	 {
		 	 	for(c6=1;c6 <= Nz;c6+=1)
		 	 	 {
		 	 	 	S12((4),(0),(0),(c4),(c5),(c6));
		 	 	 	S13((4),(0),(0),(c4),(c5),(c6));
		 	 	 }
		 	 }
		 }
		for(c4=1;c4 <= Nx;c4+=1)
		 {
		 	for(c5=1;c5 <= Ny;c5+=1)
		 	 {
		 	 	for(c6=1;c6 <= Nz;c6+=1)
		 	 	 {
		 	 	 	S16((5),(0),(0),(c4),(c5),(c6));
		 	 	 	S17((5),(0),(0),(c4),(c5),(c6));
		 	 	 }
		 	 }
		 }
		if ((Nx >= Sx && Ny >= Sy && Nz >= Sz && Sx >= 1 && Sy >= 1 && Sz >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S20((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S23((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S26((6),(0),(0),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny >= Sy+1 && Nz >= Sz+1 && Sx >= 2 && Sy >= 2 && Sz >= 2)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Sy),(c6));
				 }
				S18((6),(0),(1),(Sx),(Sy),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Sy),(c6));
				 }
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny >= Sy+1 && Nz >= Sz+1 && Sy >= 2 && Sz >= 2)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Sy),(c6));
				 }
				S18((6),(0),(1),(Nx),(Sy),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Sy),(c6));
				 }
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny == Sy && Nz >= Sz+1 && Sx >= 2 && Sz >= 2)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Ny),(c6));
				 }
				S18((6),(0),(1),(Sx),(Ny),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Ny),(c6));
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny == Sy && Nz >= Sz+1 && Sz >= 2)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Ny),(c6));
				 }
				S18((6),(0),(1),(Nx),(Ny),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Ny),(c6));
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny == Sy && Sx >= 2 && Sz == 1)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	S19((6),(0),(1),(Sx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				S18((6),(0),(1),(Sx),(Ny),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Ny),(c6));
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny == Sy && Sz == 1)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	S19((6),(0),(1),(Nx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
				S18((6),(0),(1),(Nx),(Ny),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Ny),(c6));
				 }
			}
		}
		if ((Nx >= Sx+1 && Nz >= Sz+1 && Sx >= 2 && Sy == 1 && Sz >= 2)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(1),(c6));
				 }
				S18((6),(0),(1),(Sx),(1),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(1),(c6));
				 }
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Nz >= Sz+1 && Sy == 1 && Sz >= 2)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(1),(c6));
				 }
				S18((6),(0),(1),(Nx),(1),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(1),(c6));
				 }
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Sx >= 2 && Sy == 1 && Sz == 1)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				S18((6),(0),(1),(Sx),(1),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(1),(c6));
				 }
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	S19((6),(0),(1),(Sx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Sy == 1 && Sz == 1)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				S18((6),(0),(1),(Nx),(1),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(1),(c6));
				 }
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	S19((6),(0),(1),(Nx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny >= Sy+1 && Sx >= 2 && Sy >= 2 && Sz == 1)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	S19((6),(0),(1),(Sx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				S18((6),(0),(1),(Sx),(Sy),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Sy),(c6));
				 }
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	S19((6),(0),(1),(Sx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny >= Sy+1 && Sy >= 2 && Sz == 1)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	S19((6),(0),(1),(Nx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
				S18((6),(0),(1),(Nx),(Sy),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Sy),(c6));
				 }
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	S19((6),(0),(1),(Nx),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny >= Sy+1 && Nz == Sz && Sx >= 2 && Sy >= 2)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Sy),(c6));
				 }
				S18((6),(0),(1),(Sx),(Sy),(Nz));
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny >= Sy+1 && Nz == Sz && Sy >= 2)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Sy),(c6));
				 }
				S18((6),(0),(1),(Nx),(Sy),(Nz));
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny == Sy && Nz == Sz && Sx >= 2)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(Ny),(c6));
				 }
				S18((6),(0),(1),(Sx),(Ny),(Nz));
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny == Sy && Nz == Sz)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(Ny),(c6));
				 }
				S18((6),(0),(1),(Nx),(Ny),(Nz));
			}
		}
		if ((Nx >= Sx+1 && Nz == Sz && Sx >= 2 && Sy == 1)) {
			{
				for(c4=1;c4 <= Sx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Sx),(1),(c6));
				 }
				S18((6),(0),(1),(Sx),(1),(Nz));
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Sx),(c5),(c6));
				 	 }
				 }
				for(c4=Sx+1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Nz == Sz && Sy == 1)) {
			{
				for(c4=1;c4 <= Nx-1;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(Nx),(1),(c6));
				 }
				S18((6),(0),(1),(Nx),(1),(Nz));
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(Nx),(c5),(c6));
				 	 }
				 }
			}
		}
		if ((Ny >= Sy+1 && Nz >= Sz+1 && Sx == 1 && Sy >= 2 && Sz >= 2)) {
			{
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Sy),(c6));
				 }
				S18((6),(0),(1),(1),(Sy),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Sy),(c6));
				 }
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny == Sy && Nz >= Sz+1 && Sx == 1 && Sz >= 2)) {
			{
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Ny),(c6));
				 }
				S18((6),(0),(1),(1),(Ny),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Ny),(c6));
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny >= Sy+1 && Sx == 1 && Sy >= 2 && Sz == 1)) {
			{
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	S19((6),(0),(1),(1),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				S18((6),(0),(1),(1),(Sy),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Sy),(c6));
				 }
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	S19((6),(0),(1),(1),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny == Sy && Sx == 1 && Sz == 1)) {
			{
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	S19((6),(0),(1),(1),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				S18((6),(0),(1),(1),(Ny),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Ny),(c6));
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny >= Sy+1 && Nz == Sz && Sx == 1 && Sy >= 2)) {
			{
				for(c5=1;c5 <= Sy-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Sy),(c6));
				 }
				S18((6),(0),(1),(1),(Sy),(Nz));
				for(c5=Sy+1;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny == Sy && Nz == Sz && Sx == 1)) {
			{
				for(c5=1;c5 <= Ny-1;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(Ny),(c6));
				 }
				S18((6),(0),(1),(1),(Ny),(Nz));
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz >= Sz+1 && Sx == 1 && Sy == 1 && Sz >= 2)) {
			{
				for(c6=1;c6 <= Sz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(1),(c6));
				 }
				S18((6),(0),(1),(1),(1),(Sz));
				for(c6=Sz+1;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(1),(c6));
				 }
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Sz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz == Sz && Sx == 1 && Sy == 1)) {
			{
				for(c6=1;c6 <= Nz-1;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(1),(c6));
				 }
				S18((6),(0),(1),(1),(1),(Nz));
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Sx == 1 && Sy == 1 && Sz == 1)) {
			{
				S18((6),(0),(1),(1),(1),(1));
				for(c6=2;c6 <= Nz;c6+=1)
				 {
				 	S19((6),(0),(1),(1),(1),(c6));
				 }
				for(c5=2;c5 <= Ny;c5+=1)
				 {
				 	S19((6),(0),(1),(1),(c5),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(0),(1),(1),(c5),(c6));
				 	 }
				 }
				for(c4=2;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(0),(1),(c4),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny >= Sy && Nz >= Sz && Sx >= 1 && Sy >= 1 && Sz >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S29((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S32((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S35((6),(0),(2),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S38((6),(0),(3),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny >= Sy && Sx >= 1 && Sy >= 1 && Sz <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S20((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S23((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S26((6),(0),(0),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S29((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S32((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S35((6),(0),(2),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S38((6),(0),(3),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny >= Sy && Nz <= Sz-1 && Sx >= 1 && Sy >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S20((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S23((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S26((6),(0),(0),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S29((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S32((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S35((6),(0),(2),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S38((6),(0),(3),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Sx >= 1 && Sy <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S20((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S23((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S26((6),(0),(0),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Nz >= Sz+1 && Sx >= 1 && Sy <= 0 && Sz >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Nz <= Sz && Sx >= 1 && Sy <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Sx >= 1 && Sy <= 0 && Sz <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Sx >= 1 && Sy <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S29((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S32((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S35((6),(0),(2),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S38((6),(0),(3),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Sx >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S20((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S23((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S26((6),(0),(0),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Nz >= Sz+1 && Sx >= 1 && Sz >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Nz <= Sz && Sx >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Sx >= 1 && Sz <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Sx >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S29((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S32((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S35((6),(0),(2),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S38((6),(0),(3),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Sx <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S20((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S23((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S26((6),(0),(0),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz >= Sz+1 && Sx <= 0 && Sz >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz <= Sz && Sx <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Sx <= 0 && Sz <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Sx <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S29((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S32((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S35((6),(0),(2),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S38((6),(0),(3),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S20((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S23((6),(0),(0),(c4),(c5),(c6));
				 	 	 	S26((6),(0),(0),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1 && Nz >= Sz+1 && Sz >= 1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1 && Nz <= Sz)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1 && Sz <= 0)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(0),(1),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1)) {
			{
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S29((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S32((6),(0),(2),(c4),(c5),(c6));
				 	 	 	S35((6),(0),(2),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
				for(c4=1;c4 <= Nx;c4+=1)
				 {
				 	for(c5=1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S38((6),(0),(3),(c4),(c5),(c6));
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny >= Sy+1 && Nz >= Sz+1 && Sx >= 2 && Sy >= 2 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Sy),(c6));
				 	 }
				 	S18((6),(c2),(1),(Sx),(Sy),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Sy),(c6));
				 	 }
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny >= Sy+1 && Nz >= Sz+1 && Sy >= 2 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Sy),(c6));
				 	 }
				 	S18((6),(c2),(1),(Nx),(Sy),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Sy),(c6));
				 	 }
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny == Sy && Nz >= Sz+1 && Sx >= 2 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Ny),(c6));
				 	 }
				 	S18((6),(c2),(1),(Sx),(Ny),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Ny),(c6));
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny == Sy && Nz >= Sz+1 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Ny),(c6));
				 	 }
				 	S18((6),(c2),(1),(Nx),(Ny),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Ny),(c6));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny == Sy && Sx >= 2 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(Sx),(Ny),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Ny),(c6));
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny == Sy && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(Nx),(Ny),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Ny),(c6));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Nz >= Sz+1 && Sx >= 2 && Sy == 1 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(1),(c6));
				 	 }
				 	S18((6),(c2),(1),(Sx),(1),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Nz >= Sz+1 && Sy == 1 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(1),(c6));
				 	 }
				 	S18((6),(c2),(1),(Nx),(1),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Sx >= 2 && Sy == 1 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(Sx),(1),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Sy == 1 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(Nx),(1),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny >= Sy+1 && Sx >= 2 && Sy >= 2 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(Sx),(Sy),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Sy),(c6));
				 	 }
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny >= Sy+1 && Sy >= 2 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(Nx),(Sy),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Sy),(c6));
				 	 }
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny >= Sy+1 && Nz == Sz && Sx >= 2 && Sy >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Sy),(c6));
				 	 }
				 	S18((6),(c2),(1),(Sx),(Sy),(Nz));
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny >= Sy+1 && Nz == Sz && Sy >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Sy),(c6));
				 	 }
				 	S18((6),(c2),(1),(Nx),(Sy),(Nz));
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Ny == Sy && Nz == Sz && Sx >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(Ny),(c6));
				 	 }
				 	S18((6),(c2),(1),(Sx),(Ny),(Nz));
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Ny == Sy && Nz == Sz)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(Ny),(c6));
				 	 }
				 	S18((6),(c2),(1),(Nx),(Ny),(Nz));
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx+1 && Nz == Sz && Sx >= 2 && Sy == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Sx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Sx),(1),(c6));
				 	 }
				 	S18((6),(c2),(1),(Sx),(1),(Nz));
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Sx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=Sx+1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx == Sx && Nz == Sz && Sy == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(Nx),(1),(c6));
				 	 }
				 	S18((6),(c2),(1),(Nx),(1),(Nz));
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(Nx),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny >= Sy+1 && Nz >= Sz+1 && Sx == 1 && Sy >= 2 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Sy),(c6));
				 	 }
				 	S18((6),(c2),(1),(1),(Sy),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Sy),(c6));
				 	 }
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny >= Sy+1 && Nz == Sz && Sx == 1 && Sy >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Sy),(c6));
				 	 }
				 	S18((6),(c2),(1),(1),(Sy),(Nz));
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny >= Sy+1 && Sx == 1 && Sy >= 2 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c5=1;c5 <= Sy-1;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(1),(Sy),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Sy),(c6));
				 	 }
				 	for(c5=Sy+1;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny == Sy && Nz >= Sz+1 && Sx == 1 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Ny),(c6));
				 	 }
				 	S18((6),(c2),(1),(1),(Ny),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny == Sy && Nz == Sz && Sx == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Ny),(c6));
				 	 }
				 	S18((6),(c2),(1),(1),(Ny),(Nz));
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Ny == Sy && Sx == 1 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c5=1;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	S18((6),(c2),(1),(1),(Ny),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz >= Sz+1 && Sx == 1 && Sy == 1 && Sz >= 2)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c6=1;c6 <= Sz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(1),(c6));
				 	 }
				 	S18((6),(c2),(1),(1),(1),(Sz));
				 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz == Sz && Sx == 1 && Sy == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c6=1;c6 <= Nz-1;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(1),(c6));
				 	 }
				 	S18((6),(c2),(1),(1),(1),(Nz));
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Sx == 1 && Sy == 1 && Sz == 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	S18((6),(c2),(1),(1),(1),(1));
				 	for(c6=2;c6 <= Nz;c6+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S19((6),(c2),(1),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S19((6),(c2),(1),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny >= Sy && Sx >= 1 && Sy >= 1 && Sz <= 0)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny >= Sy && Nz <= Sz-1 && Sx >= 1 && Sy >= 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Nz >= Sz+1 && Sx >= 1 && Sy <= 0 && Sz >= 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Nz <= Sz && Sx >= 1 && Sy <= 0)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Sx >= 1 && Sy <= 0 && Sz <= 0)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Nz >= Sz+1 && Sx >= 1 && Sz >= 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Nz <= Sz && Sx >= 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx >= Sx && Ny <= Sy-1 && Sx >= 1 && Sz <= 0)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz >= Sz+1 && Sx <= 0 && Sz >= 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nz <= Sz && Sx <= 0)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Sx <= 0 && Sz <= 0)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1 && Nz >= Sz+1 && Sz >= 1)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Sz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	for(c6=Sz+1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1 && Nz <= Sz)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
		if ((Nx <= Sx-1 && Sz <= 0)) {
			{
				for(c2=1;c2 <= nt-1;c2+=1)
				 {
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(1),(c6));
				 	 	S24((6),(c2),(0),(1),(1),(c6));
				 	 	S27((6),(c2),(0),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(c5),(1));
				 	 	S24((6),(c2),(0),(1),(c5),(1));
				 	 	S27((6),(c2),(0),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(1),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(1),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(1),(c5),(Nz));
				 	 	S24((6),(c2),(0),(1),(c5),(Nz));
				 	 	S27((6),(c2),(0),(1),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(1),(Ny),(c6));
				 	 	S24((6),(c2),(0),(1),(Ny),(c6));
				 	 	S27((6),(c2),(0),(1),(Ny),(c6));
				 	 }
				 	for(c4=2;c4 <= Nx-1;c4+=1)
				 	 {
				 	 	S21((6),(c2),(0),(c4),(1),(1));
				 	 	S24((6),(c2),(0),(c4),(1),(1));
				 	 	S27((6),(c2),(0),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(1),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(1),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(1),(Nz));
				 	 	S24((6),(c2),(0),(c4),(1),(Nz));
				 	 	S27((6),(c2),(0),(c4),(1),(Nz));
				 	 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(1));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 	 {
				 	 	 	 	S22((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S25((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 	S28((6),(c2),(0),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 	S21((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S24((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 	S28((6),(c2),(0),(c4),(c5),(Nz));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(1));
				 	 	S24((6),(c2),(0),(c4),(Ny),(1));
				 	 	S27((6),(c2),(0),(c4),(Ny),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S21((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S25((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 	S27((6),(c2),(0),(c4),(Ny),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S24((6),(c2),(0),(c4),(Ny),(Nz));
				 	 	S27((6),(c2),(0),(c4),(Ny),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(1),(c6));
				 	 	S24((6),(c2),(0),(Nx),(1),(c6));
				 	 	S27((6),(c2),(0),(Nx),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny-1;c5+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(c5),(1));
				 	 	S24((6),(c2),(0),(Nx),(c5),(1));
				 	 	S27((6),(c2),(0),(Nx),(c5),(1));
				 	 	for(c6=2;c6 <= Nz-1;c6+=1)
				 	 	 {
				 	 	 	S22((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S24((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 	S27((6),(c2),(0),(Nx),(c5),(c6));
				 	 	 }
				 	 	S21((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S24((6),(c2),(0),(Nx),(c5),(Nz));
				 	 	S27((6),(c2),(0),(Nx),(c5),(Nz));
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S21((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S24((6),(c2),(0),(Nx),(Ny),(c6));
				 	 	S27((6),(c2),(0),(Nx),(Ny),(c6));
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S19((6),(c2),(1),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c6=1;c6 <= Nz;c6+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(1),(c6));
				 	 	S33((6),(c2),(2),(1),(1),(c6));
				 	 	S36((6),(c2),(2),(1),(1),(c6));
				 	 }
				 	for(c5=2;c5 <= Ny;c5+=1)
				 	 {
				 	 	S30((6),(c2),(2),(1),(c5),(1));
				 	 	S33((6),(c2),(2),(1),(c5),(1));
				 	 	S36((6),(c2),(2),(1),(c5),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S31((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S33((6),(c2),(2),(1),(c5),(c6));
				 	 	 	S36((6),(c2),(2),(1),(c5),(c6));
				 	 	 }
				 	 }
				 	for(c4=2;c4 <= Nx;c4+=1)
				 	 {
				 	 	S30((6),(c2),(2),(c4),(1),(1));
				 	 	S33((6),(c2),(2),(c4),(1),(1));
				 	 	S36((6),(c2),(2),(c4),(1),(1));
				 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S34((6),(c2),(2),(c4),(1),(c6));
				 	 	 	S36((6),(c2),(2),(c4),(1),(c6));
				 	 	 }
				 	 	for(c5=2;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	S30((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S33((6),(c2),(2),(c4),(c5),(1));
				 	 	 	S37((6),(c2),(2),(c4),(c5),(1));
				 	 	 	for(c6=2;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S31((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S34((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 	S37((6),(c2),(2),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 	for(c4=1;c4 <= Nx;c4+=1)
				 	 {
				 	 	for(c5=1;c5 <= Ny;c5+=1)
				 	 	 {
				 	 	 	for(c6=1;c6 <= Nz;c6+=1)
				 	 	 	 {
				 	 	 	 	S38((6),(c2),(3),(c4),(c5),(c6));
				 	 	 	 }
				 	 	 }
				 	 }
				 }
			}
		}
	}
	#undef S0
	#undef S1
	#undef S2
	#undef S3
	#undef S4
	#undef S5
	#undef S6
	#undef S7
	#undef S8
	#undef S9
	#undef S10
	#undef S11
	#undef S12
	#undef S13
	#undef S14
	#undef S15
	#undef S16
	#undef S17
	#undef S18
	#undef S19
	#undef S20
	#undef S21
	#undef S22
	#undef S23
	#undef S24
	#undef S25
	#undef S26
	#undef S27
	#undef S28
	#undef S29
	#undef S30
	#undef S31
	#undef S32
	#undef S33
	#undef S34
	#undef S35
	#undef S36
	#undef S37
	#undef S38
	
	//Memory Free
	free(_lin_ex);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		free(ex[mz1]);
	}
	free(ex);
	
	free(_lin_ey);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		free(ey[mz1]);
	}
	free(ey);
	
	free(_lin_ez);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		free(ez[mz1]);
	}
	free(ez);
	
	free(_lin_Ex);
	for (mz1=0;mz1 < 2; mz1++) {
		for (mz2=0;mz2 < Nx+1; mz2++) {
			free(Ex[mz1][mz2]);
		}
		free(Ex[mz1]);
	}
	free(Ex);
	
	free(_lin_Ey);
	for (mz1=0;mz1 < 2; mz1++) {
		for (mz2=0;mz2 < Nx+1; mz2++) {
			free(Ey[mz1][mz2]);
		}
		free(Ey[mz1]);
	}
	free(Ey);
	
	free(_lin_Hx);
	for (mz1=0;mz1 < 2; mz1++) {
		for (mz2=0;mz2 < Nx+1; mz2++) {
			free(Hx[mz1][mz2]);
		}
		free(Hx[mz1]);
	}
	free(Hx);
	
	free(_lin_Hy);
	for (mz1=0;mz1 < 2; mz1++) {
		for (mz2=0;mz2 < Nx+1; mz2++) {
			free(Hy[mz1][mz2]);
		}
		free(Hy[mz1]);
	}
	free(Hy);
	
	free(_lin_Hz);
	for (mz1=0;mz1 < 2; mz1++) {
		for (mz2=0;mz2 < Nx+1; mz2++) {
			free(Hz[mz1][mz2]);
		}
		free(Hz[mz1]);
	}
	free(Hz);
	
	
	
	free(Cexhy);
	free(Cexhz);
	
	
	free(Ceyhx);
	free(Ceyhz);
	
	
	free(Cezhx);
	free(Cezhy);
	free(_lin_S);
	for (mz1=0;mz1 < Nx+1; mz1++) {
		free(S[mz1]);
	}
	free(S);
	
}

//Memory Macros
#undef ex
#undef ey
#undef ez
#undef Ex
#undef Ey
#undef Hx
#undef Hy
#undef Hz
#undef Cexhy
#undef Cexhz
#undef Ceyhx
#undef Ceyhz
#undef Cezhx
#undef Cezhy
#undef S
#undef Ez


//Common Macro undefs
#undef max
#undef MAX
#undef min
#undef MIN
#undef CEILD
#undef ceild
#undef FLOORD
#undef floord
#undef CDIV
#undef FDIV
#undef LB_SHIFT
#undef MOD
#undef RADD
#undef RMUL
#undef RMAX
#undef RMIN
