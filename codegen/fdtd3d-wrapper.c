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
#include <time.h>
#include <sys/time.h>
#include <sys/errno.h>


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
#define EPSILON 1.0E-9
#define M_PI 3.14159265358979323846






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

#define Ez_verify(t,i,j,k) Ez_verify[MOD(t,2)][i][j][k]
#define var_Ez(t,i,j,k) Ez(t,i,j,k)
#define var_Ez_verify(t,i,j,k) Ez_verify(t,i,j,k)

//function prototypes
void fdtd3d(long, long, long, long, long, long, long, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double****);
void fdtd3d_verify(long, long, long, long, long, long, long, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double****);

//main
int main(int argc, char** argv) {
	//Check number of args
	if (argc <= 7) {
		printf("Number of argument is smaller than expected.\n");
		printf("Expecting nt,Nx,Ny,Nz,Sx,Sy,Sz\n");
		exit(0);
	}
	
	char *end = 0;
	char *val = 0;
	//Read Parameters
	//Initialisation of nt
	errno = 0;
	end = 0;
	val = argv[1];
	long nt = strtol(val,&end,10);
	if ((errno == ERANGE && (nt == LONG_MAX || nt == LONG_MIN)) || (errno != 0 && nt == 0)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}
	if (end == val) {
		fprintf(stderr, "No digits were found for nt\n");
		exit(EXIT_FAILURE);
	}
	if (*end != '\0'){
		printf("For parameter nt: Converted part: %ld, non-convertible part: %s\n", nt, end);
		exit(EXIT_FAILURE);
	}
	
	//Initialisation of Nx
	errno = 0;
	end = 0;
	val = argv[2];
	long Nx = strtol(val,&end,10);
	if ((errno == ERANGE && (Nx == LONG_MAX || Nx == LONG_MIN)) || (errno != 0 && Nx == 0)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}
	if (end == val) {
		fprintf(stderr, "No digits were found for Nx\n");
		exit(EXIT_FAILURE);
	}
	if (*end != '\0'){
		printf("For parameter Nx: Converted part: %ld, non-convertible part: %s\n", Nx, end);
		exit(EXIT_FAILURE);
	}
	
	//Initialisation of Ny
	errno = 0;
	end = 0;
	val = argv[3];
	long Ny = strtol(val,&end,10);
	if ((errno == ERANGE && (Ny == LONG_MAX || Ny == LONG_MIN)) || (errno != 0 && Ny == 0)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}
	if (end == val) {
		fprintf(stderr, "No digits were found for Ny\n");
		exit(EXIT_FAILURE);
	}
	if (*end != '\0'){
		printf("For parameter Ny: Converted part: %ld, non-convertible part: %s\n", Ny, end);
		exit(EXIT_FAILURE);
	}
	
	//Initialisation of Nz
	errno = 0;
	end = 0;
	val = argv[4];
	long Nz = strtol(val,&end,10);
	if ((errno == ERANGE && (Nz == LONG_MAX || Nz == LONG_MIN)) || (errno != 0 && Nz == 0)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}
	if (end == val) {
		fprintf(stderr, "No digits were found for Nz\n");
		exit(EXIT_FAILURE);
	}
	if (*end != '\0'){
		printf("For parameter Nz: Converted part: %ld, non-convertible part: %s\n", Nz, end);
		exit(EXIT_FAILURE);
	}
	
	//Initialisation of Sx
	errno = 0;
	end = 0;
	val = argv[5];
	long Sx = strtol(val,&end,10);
	if ((errno == ERANGE && (Sx == LONG_MAX || Sx == LONG_MIN)) || (errno != 0 && Sx == 0)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}
	if (end == val) {
		fprintf(stderr, "No digits were found for Sx\n");
		exit(EXIT_FAILURE);
	}
	if (*end != '\0'){
		printf("For parameter Sx: Converted part: %ld, non-convertible part: %s\n", Sx, end);
		exit(EXIT_FAILURE);
	}
	
	//Initialisation of Sy
	errno = 0;
	end = 0;
	val = argv[6];
	long Sy = strtol(val,&end,10);
	if ((errno == ERANGE && (Sy == LONG_MAX || Sy == LONG_MIN)) || (errno != 0 && Sy == 0)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}
	if (end == val) {
		fprintf(stderr, "No digits were found for Sy\n");
		exit(EXIT_FAILURE);
	}
	if (*end != '\0'){
		printf("For parameter Sy: Converted part: %ld, non-convertible part: %s\n", Sy, end);
		exit(EXIT_FAILURE);
	}
	
	//Initialisation of Sz
	errno = 0;
	end = 0;
	val = argv[7];
	long Sz = strtol(val,&end,10);
	if ((errno == ERANGE && (Sz == LONG_MAX || Sz == LONG_MIN)) || (errno != 0 && Sz == 0)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}
	if (end == val) {
		fprintf(stderr, "No digits were found for Sz\n");
		exit(EXIT_FAILURE);
	}
	if (*end != '\0'){
		printf("For parameter Sz: Converted part: %ld, non-convertible part: %s\n", Sz, end);
		exit(EXIT_FAILURE);
	}
	
	
	///Parameter checking
	if (!((Nx >= 10 && Ny >= 10 && Nz >= 10 && nt >= 1))) {
		printf("The value of parameters are not valid.\n");
		exit(-1);
	}
	
	// These are the parameters that control the physical properties of the simulated env    
	double f0 = 1e6;
	double Lf = 10;
	double CF = 0.99;
	double er_x = 6;
	double er_y = 6;
	double er_z = 6;
	double Lx = Nx/Lf;
	double Ly = Ny/Lf;
	double Lz = Nz/Lf;
	double e0 = 8.854e-12;
	double u0 = 4 * M_PI * (1e-07);
	double c0 = sqrt(1/(e0*u0));
	double dx = c0/f0/sqrt(er_x)/Lf;
	double dy = c0/f0/sqrt(er_y)/Lf;
	double dz = c0/f0/sqrt(er_z)/Lf ;
	double dt = (sqrt(1/((1/(dx*dx))+(1/(dy*dy))+(1/(dz*dz)))) / c0);
	dt = dt*CF;
	double T = 1/(f0*dt);
	dt = dt*T/ceil(T);
	
	printf("nt: %ld\n", nt);
	printf("Nx: %ld\n", Nx);
	printf("Ny: %ld\n", Ny);
	printf("Nz: %ld\n", Nz);
	printf("Sx: %ld\n", Sx);
	printf("Sy: %ld\n", Sy);
	printf("Sz: %ld\n\n", Sz);
	
	printf("er_x: %E\n", er_x);
	printf("er_y: %E\n", er_y);
	printf("er_z: %E\n", er_z);
	printf("e0: %E\n", e0);
	printf("u0: %E\n", u0);
	printf("dx: %E\n", dx);
	printf("dy: %E\n", dy);
	printf("dz: %E\n", dz);
	printf("dt: %E\n", dt);

	//Memory Allocation
	int mz1, mz2, mz3, mz4;	

	double* _lin_Ez = (double*)malloc(sizeof(double)*(2 * (Nx+1) * (Ny+1) * (Nz+1)));
	mallocCheck(_lin_Ez, (2 * (Nx+1) * (Ny+1) * (Nz+1)), double);
	double**** Ez = (double****)malloc(sizeof(double***)*(2));
	mallocCheck(Ez, (2), double***);
	for (mz1=0;mz1 < 2; mz1++) {
		Ez[mz1] = (double***)malloc(sizeof(double**)*(Nx+1));
		mallocCheck(Ez[mz1], (Nx+1), double**);
		for (mz2=0;mz2 < Nx+1; mz2++) {
			Ez[mz1][mz2] = (double**)malloc(sizeof(double*)*(Ny+1));
			mallocCheck(Ez[mz1][mz2], (Ny+1), double*);
			for (mz3=0;mz3 < Ny+1; mz3++) {
				Ez[mz1][mz2][mz3] = &_lin_Ez[(mz1*((Nx+1) * (Ny+1) * (Nz+1))) + (mz2*((Ny+1) * (Nz+1))) + (mz3*(Nz+1))];
			}
		}
	}
	#ifdef VERIFY
		double* _lin_Ez_verify = (double*)malloc(sizeof(double)*(2 * (Nx+1) * (Ny+1) * (Nz+1)));
		mallocCheck(_lin_Ez_verify, (2 * (Nx+1) * (Ny+1) * (Nz+1)), double);
		double**** Ez_verify = (double****)malloc(sizeof(double***)*(2));
		mallocCheck(Ez_verify, (2), double***);
		for (mz1=0;mz1 < 2; mz1++) {
			Ez_verify[mz1] = (double***)malloc(sizeof(double**)*(Nx+1));
			mallocCheck(Ez_verify[mz1], (Nx+1), double**);
			for (mz2=0;mz2 < Nx+1; mz2++) {
				Ez_verify[mz1][mz2] = (double**)malloc(sizeof(double*)*(Ny+1));
				mallocCheck(Ez_verify[mz1][mz2], (Ny+1), double*);
				for (mz3=0;mz3 < Ny+1; mz3++) {
					Ez_verify[mz1][mz2][mz3] = &_lin_Ez_verify[(mz1*((Nx+1) * (Ny+1) * (Nz+1))) + (mz2*((Ny+1) * (Nz+1))) + (mz3*(Nz+1))];
				}
			}
		}
	#endif

	//Initialization of rand
	srand((unsigned)time(NULL));
	 
	//Timing
	struct timeval time;
	double elapsed_time;
	
	//Call the main computation
	gettimeofday(&time, NULL);
	elapsed_time = (((double) time.tv_sec) + ((double) time.tv_usec)/1000000);
	
	fdtd3d(nt, Nx, Ny, Nz, Sx, Sy, Sz, &f0, &Lf, &CF, &er_x, &er_y, &er_z, &Lx, &Ly, &Lz, &e0, &u0, &c0, &dx, &dy, &dz, &dt, Ez);

	gettimeofday(&time, NULL);
	elapsed_time = (((double) time.tv_sec) + ((double) time.tv_usec)/1000000) - elapsed_time;

	// timing information
	printf("Execution time : %lf sec.\n", elapsed_time);
	
	#ifdef TIMING
		FILE * fp = fopen( "trace.dat","a+");
		if (fp == NULL) {
				printf("I couldn't open trace.dat for writing.\n");
				exit(EXIT_FAILURE);
		}
		fprintf(fp, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%lf\n",nt,Nx,Ny,Nz,Sx,Sy,Sz,elapsed_time);
		fclose(fp);
	#endif
	
	//Verification Run
	#ifdef VERIFY
		#ifdef TIMING
			gettimeofday(&time, NULL);
			elapsed_time = (((double) time.tv_sec) + ((double) time.tv_usec)/1000000);
		#endif
    	fdtd3d_verify(nt, Nx, Ny, Nz, Sx, Sy, Sz, &f0, &Lf, &CF, &er_x, &er_y, &er_z, &Lx, &Ly, &Lz, &e0, &u0, &c0, &dx, &dy, &dz, &dt, Ez_verify);
    	#ifdef TIMING
    		gettimeofday(&time, NULL);
			elapsed_time = (((double) time.tv_sec) + ((double) time.tv_usec)/1000000) - elapsed_time;
			
			FILE * fp_verify = fopen( "trace_verify.dat","a+");
			if (fp == NULL) {
					printf("I couldn't open trace_verify.dat for writing.\n");
					exit(EXIT_FAILURE);
			}
			fprintf(fp, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%lf\n",nt,Nx,Ny,Nz,Sx,Sy,Sz,elapsed_time);
			fclose(fp_verify);
		#endif
	#endif
    	
	#ifdef CHECKING
    	//Print Outputs
		
		{
			#ifdef NO_PROMPT
				#define S0(t,i,j,k) printf("%0.2lf\n",var_Ez(t,i,j,k))
			#else
				#define S0(t,i,j,k) printf("Ez(%ld,%ld,%ld,%ld)=",(long) t,(long) i,(long) j,(long) k);printf("%0.2lf\n",var_Ez(t,i,j,k))
			#endif
			int c1,c2,c3,c4;
			for(c1=0;c1 <= nt-1;c1+=1)
			 {
			 	for(c2=1;c2 <= Nx;c2+=1)
			 	 {
			 	 	for(c3=1;c3 <= Ny;c3+=1)
			 	 	 {
			 	 	 	for(c4=1;c4 <= Nz;c4+=1)
			 	 	 	 {
			 	 	 	 	S0((c1),(c2),(c3),(c4));
			 	 	 	 }
			 	 	 }
			 	 }
			 }
			#undef S0
		}
	#elif VERIFY
		//Compare outputs for verification
		{
			//Error Counter
			int _errors_ = 0;
			#define S0(t,i,j,k) if (fabs(1.0 - var_Ez_verify(t,i,j,k)/var_Ez(t,i,j,k)) > EPSILON) _errors_++;
			int c1,c2,c3,c4;
			for(c1=0;c1 <= nt-1;c1+=1)
			 {
			 	for(c2=1;c2 <= Nx;c2+=1)
			 	 {
			 	 	for(c3=1;c3 <= Ny;c3+=1)
			 	 	 {
			 	 	 	for(c4=1;c4 <= Nz;c4+=1)
			 	 	 	 {
			 	 	 	 	S0((c1),(c2),(c3),(c4));
			 	 	 	 }
			 	 	 }
			 	 }
			 }
			#undef S0
			if(_errors_ == 0){
				printf("TEST for Ez PASSED\n");
			}else{
				printf("TEST for Ez FAILED. #Errors: %d\n", _errors_);
			}
		}
    #endif
    
	//Memory Free
	free(_lin_Ez);
	for (mz1=0;mz1 < 2; mz1++) {
		for (mz2=0;mz2 < Nx+1; mz2++) {
			free(Ez[mz1][mz2]);
		}
		free(Ez[mz1]);
	}
	free(Ez);
	#ifdef VERIFY
		free(_lin_Ez_verify);
		for (mz1=0;mz1 < 2; mz1++) {
			for (mz2=0;mz2 < Nx+1; mz2++) {
				free(Ez_verify[mz1][mz2]);
			}
			free(Ez_verify[mz1]);
		}
		free(Ez_verify);
	#endif
	
	return EXIT_SUCCESS;
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
#undef EPSILON
