#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*********** debug ***********/
//#define debug

#if defined (_WIN32) || defined(_WIN64)
    #define path "output/"
#elif __linux__
    #define path "output\\"
#elif defined __APPLE__
    #define path "output/"
#else
    #error "Unknown OS..."
#endif

// Defining the built-in functions to get the max of two values and min of two values.
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//#define DEBUGPRINT(...)       printf(__VA_ARGS__);   fflush(stdout)
#define FLUID
/// Constants defined to calculate the profiles
#define NSTAT 2500 // Number of divisions in the profile in function of z
#define dz 0.1L // For now it is a constant, but can be funtion of z it self
const long LISTAT = (int)(3.0L/dz); // Inferior limit to negative values of z
#define MAXGRAINS 25000 // Number of maximum grains

#ifdef FLUID
const double long PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230;
#endif

/*********** declaration functions ***********/
double long ran0();
double long ran2(long *idum);
double gaussian(long *idum);
//double gaussian2(long *idum);
double long distperio(double long x);
void findneighbours();
void findneighbours2();
void closeall();
void init();
void loadparameters(const char *);
void loadgrains(const char *);
void savegrains(const char *);
void savegrainsprofile(const char *);
void extractwall();
void prepareconfig();
void predictor();
void detectcontacts();
void forcecalculation();
void corrector();
void cumulatestats();
void timemeasurement();
void periodsave();
void saveconfiguration();
#ifdef FLUID
void loadfluid(const char *);
void savefluid(const char *);
void calculateprofiles();
void solvefluidvelocity(); // Function that solves the fluid velocity, a partial differential equation (a non-linear equation in diffusion family).
void solvefluidvelocity2(); // Function that solves the fluid velocity, a partial differential equation (a non-linear equation in diffusion family).
void fluidforce(long); // Function that calculate the force that the fluid exert over a grain. The parameter passed is the identification number of the grain.
void sheargrainswithfluid(); // Function that imposes the fluid velocity, increasing the slope of a linear fluid profile from the very bottom, when reach the middle of the simulation, it deslocates the linear fluid profile up til the top of granular layer
long double U(long double, long double, long double, long double); // Function to be integrated by Runge-Kutta to give u from du/dz
long double L(long double, long double, long double, long double); // Function to be integrated by Runge-Kutta to give l from dl/dz
long double RK4(long double, long double, long double, long double, long double (*)(long double, long double, long double, long double));
long double volumefraction(long, long double, long double);
long double areafraction(long, long double, long double);
void fraction(long, long double, long double, long double &, long double &);
long double Phibase();
#endif

/*********** global variables ***********/

char periodsaveout[30];
long it,isem = 987654321;

long extract;
long ncontfree,ncontot,nl,nlt,ifreqv,icont0,ncontglis;
long isave; //period over which configuration (state variables) AND cumulated statistics are saved
long itimemeasurements; //period over which time dependent quantities are measured and saved
long iperiodsave; //period over which are saved fields needed for Ep, \Gamma calculation
long icustats; //period over which quantities are evaluated to get cumulated statistics
double long ivmass[MAXGRAINS],ivmomine[MAXGRAINS],sumfn[MAXGRAINS*10],reacn[MAXGRAINS*10], reactsave[MAXGRAINS*10],react[MAXGRAINS*10], ax[MAXGRAINS],ay[MAXGRAINS],ome[MAXGRAINS],domedt[MAXGRAINS],phi[MAXGRAINS];
double long z[NSTAT], Zglis[NSTAT], Zcont[NSTAT], Phi[NSTAT], Vx[NSTAT], Vy[NSTAT], Fx[NSTAT], Fy[NSTAT], P[NSTAT], interp[NSTAT];
double long kn; //normal spring constant
double long kt; //tangential spring constant
double long dt; //time step
double long gn; //normal damping factor
double long frott; //friction coefficient
double long halfdt2,halfdt,tolff;
long nwall; //number of grains in each wall;
double long aa; //fit parameter fo P/phi for menu==4
double long bb; //fit parameter fo P/phi for menu==4
long nfg; //number of free grains
long ngtot; //total number of grains
long menu; //select the configuration
double long totalweight;
double long width; //cell width used for periodic boundary conditions
double long halfwidth; //half the width
double long press; // pressure imposed to the wall, when applicable
double long shear; // shear stress imposed to the wall, when applicable
double long theta; // angle along which gravity is applied, when applicable
double long costheta,sintheta;
double long gravity;
double long muwall; // pseudo friction rescaled by 2Pi
double long wallv;  //wall velocity, when applicable
double long timewall,wallvzero; //timescale over which the wall velocity is brought to 0
double long height;  //cell height, when applicable
double long bottom;  //cell bottom position, when applicable
double long wallx;  //wall displacement
double long walla,wallap,wallxp,wallvp;
double long dhdt,ddhdtt,hp,dhdtp,ddhdttp;  //cell height time derivative
double long vx[MAXGRAINS], vy[MAXGRAINS], radius[MAXGRAINS], rx[MAXGRAINS], ry[MAXGRAINS], rpx[MAXGRAINS], rpy[MAXGRAINS], vpx[MAXGRAINS], vpy[MAXGRAINS], xwall[500], ywall[500];
long ior[MAXGRAINS*10], iex[MAXGRAINS*10];
long listi[MAXGRAINS*10], listj[MAXGRAINS*10];
double long apx[MAXGRAINS], apy[MAXGRAINS], phip[MAXGRAINS], omep[MAXGRAINS], domedtp[MAXGRAINS], nZcont[MAXGRAINS*10], nZglis[MAXGRAINS*10], eij[MAXGRAINS*10], xnij[MAXGRAINS*10], ynij[MAXGRAINS*10];
double long fpx[MAXGRAINS], fpy[MAXGRAINS], gam[MAXGRAINS];
long io[MAXGRAINS*10];
long endreached;
long nq,igr;
double long deltar; // width of the uniform distribution used for the radius
long islu; //seed for the random number generator
long nperiodsave=0; //index of the current file saving configurations
double long fydown,fxdown,fyup,fxup; //Normal and tangiential stresses on the walls
double long springxx,springyy,springxy;
// The 2 variables, density and viscosity, are used to describle the fliud. All considerations were based on Phys. Fluids 24, 103306
#ifdef FLUID
double long viscosity;
// The next 3 variables are control parameters to the system:
double long reynolds, density_ratio, shields; // Reynolds number in grain dimensions, funciton of density ratio (fluid and grain), gravity, grain diameter and viscosity. Shields number, funciton of density (fluid and grain), the velocity field of the fluid, gravity and grain diameter
long H = NSTAT; // Counter to find the shear applied position
// The next 2 constants are the only constants that defines the fluid it self
const double long dragd = 0.5L; // Drag coefficient of the grain in the turbulent regime (Granular Reynolds -> inf)
const double long reynoldsc = 24.0L; // Reynolds number where drag coefficient is almost constant
const double long RvD = 26.0L; // van Driest's Reynolds number
// Some other variables to calculate the force and other parameters of the interaction between fluid and grain
double long ustar; // Characteristic velocity of the fluid - something like a shear velocity
double long z0; // Characteristic length of turbulence effect
const double long kapa = 0.4L; // Characteristic dimensionless constant of the fluid velocity field
double long zb; // The height that bed forms
double long reynoldsg; // The Reynolds number on the grain interaction
double long grain_density; // The density of the grain (Associated with the mass of the grain)
double long fluid_density; // The density of the fluid
double long u[NSTAT]; // Fluid velocity profile.
double long fxfluid_arch[MAXGRAINS], fyfluid_arch[MAXGRAINS], fxfluid_drag[MAXGRAINS], fyfluid_drag[MAXGRAINS]; // The force the grains exert on the fluid, in x and z directions
double long fxFd[NSTAT], fyFd[NSTAT]; // Drag stress that fluid exert on grains, in x and z directions
double long fxFa[NSTAT], fyFa[NSTAT]; // Archimedes stress that fluid exert on grains, in x and z directions
double long tau[NSTAT]; // Shear stress divided by density along z direction
double long dpressuredz[NSTAT]; // Derivative of the pressure along z direction
#endif
FILE *out;
FILE *out1;
FILE *outest;
long nbiteration;
double long max_x, max_y, min_x, min_y; // Greater and lower positions of particles in x and y directions

void loadparameters(const char * name){
// This function only loads what is inside the configuration parameters file ("configin.txt")

    FILE *in;
    char comment[255];

    if ((in = fopen(name, "rt"))== NULL){
        fprintf (out1,"#no input file (%s)...........\n", name);
        exit (1);
    }
    else
        fprintf(out1,"#input file detected (%s)\n", name);
    fscanf (in,"%ld\t%s\n", &nfg,comment);
    fprintf(out1,"#free grain number=%ld \n",nfg);
    fscanf (in,"%ld\t%s\n", &nwall,comment);
    fprintf(out1,"#number of grains in each wall=%ld \n",nwall);
    fscanf (in,"%Le\t%s\n", &width,comment);
    fprintf(out1,"#cell width (unit d) width=%Le \n",width);
    halfwidth=0.5L*width;
    fscanf(in,"%Le\t%s\n", &deltar,comment);
    fprintf(out1,"#Radius picked at random with a uniform law of width %Le\n",deltar);
    fscanf (in,"%Le\t%s\n", &kn,comment);
    fprintf(out1,"#normal spring constant(?) kn=%Le \n",kn);
    fscanf (in,"%Le\t%s\n", &kt,comment);
    fprintf(out1,"#kt=%Le \n",kt);
    fscanf (in,"%Le\t%s\n", &gn,comment);
    fprintf(out1,"#damping coefficient=%Le normally smaller than %Le \n",gn,2.0L*sqrt(kn));
    fscanf (in,"%Le\t%s\n", &frott,comment);
    fprintf(out1,"#friction coefficient frott=%Le \n",frott);
    fscanf (in,"%ld\t%s\n", &itimemeasurements,comment);
    fprintf(out1,"#time dependent measurement period=%ld \n",itimemeasurements);
    fscanf (in,"%ld\t%s\n", &iperiodsave,comment);
    fprintf(out1,"#fraction of dot_gamma for periodic meseasurements=%ld \n",iperiodsave);
    fscanf (in,"%ld\t%s\n", &ifreqv,comment);
    fprintf(out1,"#neighbour list refreshing period ifreqv=%ld \n",ifreqv);
    fscanf (in,"%ld\t%s\n", &icustats,comment);
    fprintf(out1,"#statistical quantities evaluation period=%ld \n",icustats);
    fscanf (in,"%ld\t%s\n", &isave,comment);
    fprintf(out1,"#configuration and statistical quantities saving period=%ld \n",isave);
    fscanf (in,"%ld\t%s\n", &menu,comment);
    fscanf (in,"%ld\t%s\n", &extract,comment);
    fprintf(out1,"#Extraction=%ld \n",extract);
    fscanf (in,"%ld\t%s\n", &nbiteration,comment);
    fprintf(out1,"#nbiteration=%ld \n",nbiteration);
#ifdef FLUID
    if ((menu >= 30)&&(menu < 45)){
        gravity = 1.0L;
        density_ratio = 2.0L;
        reynolds = 10.0L;
        shields = 0.0L;
        fscanf(in,"%Le\t%s\n", &density_ratio, comment);
        fprintf(out1,"#Density ratio of the grain over fluid=%Le\n", density_ratio);
        fscanf(in,"%Le\t%s\n", &reynolds, comment);
        fprintf(out1,"#Fluid's Reynolds number=%Le\n", reynolds);
        fscanf(in,"%Le\t%s\n", &shields, comment);
        fprintf(out1, "#Shields number=%Le\n", shields);
        fscanf(in,"%Le\t%s\n", &height, comment);
        fprintf(out1, "#Height to impose shear after bed=%Le\n", height);
        if (menu == 30)
            shields = 0.0L;
    }
#endif
    if (menu<15){
        fscanf (in,"%Le\t%s\n", &theta,comment);
        if ((menu==2) || (menu==3) || (menu==4))
            fscanf (in,"%Le\t%s\n", &muwall,comment);
        if (menu==4){
            fscanf (in,"%Le\t%s\n", &aa,comment);
            fscanf (in,"%Le\t%s\n", &bb,comment);
        }
    } else if ((menu>14)&&(menu<28)) {
        fscanf (in,"%Le\t%s\n", &wallv,comment);
        if (menu==25)
            fscanf (in,"%Le\t%s\n", &gravity,comment);
        if (menu==26)
            fscanf (in,"%Le\t%s\n", &gravity,comment);
        if (menu==27){
            wallvzero=wallv;
            fscanf (in,"%Le\t%s\n", &timewall,comment);
        }
    } else if (menu==28)
        fscanf (in,"%Le\t%s\n", &wallv,comment);
    else if (menu==29)
        fscanf (in,"%Le\t%s\n", &shear,comment);
    fclose(in);
}

void loadgrains(const char * name){
// This function loads all grains and its contacts to the system, from binary file or from text file
    FILE *in2;

    if ((in2 = fopen(name, "rt")) == NULL){
        long i,ntemp;
        fprintf (out1,"#no binary file (%s) -> old text file...........\n", name);
        if ((in2 = fopen("datain.txt","r")) == NULL){
            fprintf (out1,"#no input file (%s)...........\n", "datain.txt");
            exit (1);
        }
        fscanf(in2,"%ld\n",&ntemp);
        if (ntemp!=ngtot){
            fprintf(out1,"\n\n\n\n\n\nProblème grave dans le nombre de points: ntemp=%li ngtot=%li\n\n\n\n\n\n\n\n\n\n",ntemp, ngtot);
            fprintf(stderr, "Problème grave dans le nombre de points - Configuration file: %li Grain file %li\n", ngtot, ntemp);
        }
        for (i = 1; i<=ngtot; i++)
            fscanf(in2,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", &radius[i],&rx[i], &ry[i],&vx[i], &vy[i],&ome[i],&ax[i], &ay[i],&domedt[i]);
        fscanf(in2,"%ld\n",&ncontfree);
        fscanf(in2,"%ld\n",&ncontot);
        for (i = 1;i<= ncontot;i++)
            fscanf(in2,"%ld\t%ld\t%Le\n", &ior[i],&iex[i],&react[i]);
        fprintf(out1,"\n\nEn entrée, ngtot=%ld\n",ntemp);
        fprintf(out1," ncontfree=%ld  et    %ld \n",ncontfree,ncontot);
    } else {
        long ntemp;
        fprintf (out1,"#Binary file detected...........\n");
        fread(&ntemp,sizeof(ntemp),1,in2);
        if (ntemp!=ngtot){
            fprintf(out1,"\n\n\n\n\n\n#Problème grave dans le nombre de points: ntemp=%li ngtot=%li\n\n\n\n\n\n\n\n\n\n",ntemp, ngtot);
            fprintf(stderr, "Problème grave dans le nombre de points - Configuration file: %li Grain file %li\n", ngtot, ntemp);
        }
        fprintf(out1,"\n\n#En entrée, ngtot=%ld\n",ntemp);
        fread(radius,sizeof(double long),ngtot+1,in2);
        fread(rx,sizeof(double long),ngtot+1,in2);
        fread(ry,sizeof(double long),ngtot+1,in2);
        fread(vx,sizeof(double long),ngtot+1,in2);
        fread(vy,sizeof(double long),ngtot+1,in2);
        fread(ome,sizeof(double long),ngtot+1,in2);
        fread(ax,sizeof(double long),ngtot+1,in2);
        fread(ay,sizeof(double long),ngtot+1,in2);
        fread(domedt,sizeof(double long),ngtot+1,in2);
        fread(&ncontfree,sizeof(long),1,in2);
        fread(&ncontot,sizeof(long),1,in2);
        fprintf(out1,"# ncontfree=%ld  et    %ld \n",ncontfree,ncontot);
        fread(ior,sizeof(long),ncontot+1,in2);
        fread(iex,sizeof(long),ncontot+1,in2);
        fread(react,sizeof(double long),ncontot+1,in2);
    }
    fclose(in2);
}

#ifdef FLUID
void loadfluid(const char * name){
    FILE *config_fluid;
    char tmp[200];
    long j = 0, k = 0;

    if ((config_fluid = fopen(name,"rt")) == NULL){
        fprintf (out1,"#no binary file (%s) -> old text file...........\n", name);
        if ((config_fluid = fopen("fluidin.txt","r")) == NULL){
            fprintf(out1, "#no input file (%s)...........\n", "fluidin.txt");
            exit(1);
        }
        fgets(tmp, 200, config_fluid);
        while(!feof(config_fluid)){
            fscanf(config_fluid, "%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", &z[j], &u[j], &tau[j], &dpressuredz[j], &fxFd[j], &fyFd[j]);
            dpressuredz[j] *= fluid_density;
            if (Phi[j] != 0.0L){
                fxFd[j] *= fluid_density*(1.0L-Phi[j])/Phi[j];
                fxFd[j] *= fluid_density*(1.0L-Phi[j])/Phi[j];
            } else {
                fxFd[j] = 0.0L;
                fxFd[j] = 0.0L;
            }
            j+= (j < NSTAT) ? 1 : 0;
        }
    } else {
        long ntemp;
        fprintf (out1,"#Binary fluid file detected...........\n");
        fread(&ntemp,sizeof(ntemp),1,config_fluid);
        fprintf(out1,"\n\n#En entrée, H=%ld\n",ntemp);
        fread(u,sizeof(double long),ntemp,config_fluid);
        fread(tau,sizeof(double long),ntemp,config_fluid);
        fread(dpressuredz,sizeof(double long),ntemp,config_fluid);
        fread(fxFd,sizeof(double long),ntemp,config_fluid);
        fread(fyFd,sizeof(double long),ntemp,config_fluid);
        fread(Phi,sizeof(double long),ntemp,config_fluid);
        j = ntemp;
    }
    fclose(config_fluid);
    if (j != H){
        fprintf(stderr, "WARNING - Fluid config does not match with the fluid's mesh. Number of lines read %li, size of the mesh %li, maximum mesh points %i.\n", j, H, NSTAT);
        if (j > H){
            fprintf(stderr, "WARNING - Imposing shear at %Le and ignoring last %li points.\n",(H-LISTAT)*dz, j-H);
        } else {
            fprintf(stderr, "WARNING - Imposing shear at %Le and extrapolating %li points.\n",(H-LISTAT)*dz, H-j);
            for (k = j; k < H; k++){
                u[k] = RK4((k-LISTAT)*dz, u[k-1], 0.0L, 0.0L, &U);
                tau[k] = ustar*ustar;
                dpressuredz[k] = -fluid_density*gravity;
            }
        }
    }
}
#endif

void init(){
/* In this simulation program, the menu is used to choice which kind of simulation will be done:
From 0 to 14, one bottom inclined wall with gravity is used to simulate the sistem.
From 15 to 29, two walls, one at bottom and other at top, without gravity, are used to simulate shear systems.
Form 30 to 44, one bottom horizontal wall with gravity is used to simulate fluid systems.
Funcitons 0, 15, and 30 create systems. */
    long i, k;

    if ((out1 = fopen ("output.tsv","w")) == NULL){
        fprintf(stderr,"Probleme d'ouverture de %s\n","output.tsv");
        exit (1);
    }
    srand(time(NULL));
    islu=rand();

    loadparameters("configin.txt");

    if ((menu>14) && (menu < 30))
        ngtot=nfg+2*nwall;
    else
        ngtot=nfg+nwall;
    fprintf(out1,"#Total number of grains=%ld \n",ngtot);
    if (menu<15){
        fprintf(out1,"#theta=%Le \n",theta);
        costheta=cos(theta);
        sintheta=sin(theta);
        fprintf(out1,"#Inclined plane configuration\n");
        if ((menu==2) || (menu==3) || (menu==4)){
            fprintf(out1,"#Includes pseudo-wall friction\n");
            fprintf(out1,"#muwall=%Le \n",muwall);
            muwall/=(2.0L*M_PI);
        }
        if (menu==4) {
            fprintf(out1,"#Includes fit expression of P/phi in pseudo-wall friction\n");
            fprintf(out1,"#aa=%Le \n",aa);
            fprintf(out1,"#bb=%Le \n",bb);
        }
    } else if ((menu>14)&&(menu<28)) {
        fprintf(out1,"#Pressure-velocity controlled configuration\n");
        fprintf(out1,"#wallv U=%Le \n",wallv);
        if (menu==25) {
            fprintf(out1,"#Gravity/anti-gravity added in five layers on both sides\n");
            fprintf(out1,"#gravity=%Le \n",gravity);
        } else if (menu==26){
            fprintf(out1,"#Gravity added\n");
            fprintf(out1,"#gravity=%Le \n",gravity);
        } else if (menu==27)
            fprintf(out1,"#time-scale decay=%Le \n",timewall);
    } else if (menu==28) {
        fprintf(out1,"#Volume-velocity controlled configuration\n");
        fprintf(out1,"#wallv U=%Le \n",wallv);
    } else if (menu==29) {
        fprintf(out1,"#Pressure-Stress controlled configuration\n");
        fprintf(out1,"#shear Tau=%Le \n",shear);
    }

// Again, all the configurations to be tested, this time to see if it's to prepare or to load
    if ((menu==0)||(menu==15)||(menu==30)) {
        if (menu==0)
            fprintf(out1,"#menu 0: prepare from scratch with bottom wall) \n");
        else if (menu==15)
            fprintf(out1,"#menu 15: prepare from scratch with two walls) \n");
#ifdef FLUID
        else if (menu==30)
            fprintf(out1,"#menu 30: prepare from scratch with bottom wall) \n");
#endif
        fprintf(out1,"-------------------------------------------------\n");
        //CHECK ADD HERE ROUGHNESS PARAMETERS
        prepareconfig();
    } else {
        fprintf(out1,"#Start from a saved configuration\n");
        loadgrains("datain.bin");
    }

    totalweight=0.0L;

/********* Mass normalisation*/
#ifndef FLUID
    for (i=1;i<=nfg;i++) {
        totalweight+= 4.0L*radius[i]*radius[i];
        ivmass[i] = 0.25L / (radius[i]*radius[i]); //inverse of the grain mass
        ivmomine[i] = 2.0L * ivmass[i] /(radius[i]*radius[i]); //inverse of the grain moment of inertia
    }
#else
    for (i=1;i<=nfg;i++){
        ivmass[i] = 1.0L/(8.0L*radius[i]*radius[i]*radius[i]); // Inverse of the grain mass (normalizated)
        totalweight+=1.0L/ivmass[i]; // Density of grain != 1 (6/PI)
        ivmomine[i] = 2.5L * ivmass[i] /(radius[i]*radius[i]); // Inverse of the grain moment of inertia of a sphere
    }
    grain_density = 6.0L/PI; // Considering a sphere as volume, with diameter=1 and mass=1
    fluid_density = grain_density/density_ratio;
#endif
    fprintf(out1,"#Mass =%Le\n",totalweight);
// Initialization of Granular variables of momentum
    for (i=nfg+1;i<=ngtot;i++) {
        rpx[i]=rx[i];
        rpy[i]=ry[i];
        vpx[i]=0.0L;
        vpy[i]=0.0L;
        omep[i]=0.0L;
        vx[i]=0.0L;
        vy[i]=0.0L;
        ome[i]=0.0L;
        phi[i] = 0.0L;
    }
// Initialization of Fluid variables and constants
#ifdef FLUID
    zb = 0.0L;
    long double zbmax = 1e20, zbase;
    max_y = 0.0L;
    for (i = 1; i <= ngtot; i++){
        max_y = max(max_y, ry[i]+radius[i]);
        if (sqrt(vx[i]*vx[i]+vy[i]*vy[i]) < 0.01L){
            if ((zb < ry[i])&&(ry[i] < zbmax))
                zb = ry[i];
        } else {
            if (zbmax > ry[i])
                zbmax = ry[i];
        }
    }
    if ((menu >= 30) && (menu < 45)) {
        viscosity = sqrt((density_ratio -1.0L)*gravity)/reynolds;
        ustar = sqrt(shields*(density_ratio -1.0L)*gravity);
// Finding a position to apply shear into the fluid
        for (k = 0; k < NSTAT; k++){
            Phi[k] = 0.0L;
        }
// First, let us construct packing fraction profile
        for (i=1;i<=ngtot;i++){
            long begin = floor((ry[i]-radius[i])/dz);
            long end = ceil((ry[i]+radius[i])/dz);
            for (k = begin; k < end; k++){
                long double volumeratio = volumefraction(i, k*dz, (k+1)*dz)/(width*dz);
                // To be sure that k is not invalidating any data:
                if ((k >= -LISTAT) && (k < NSTAT-LISTAT)){
                    Phi[LISTAT+k] += volumeratio;
/*if(isinf(Phi[LISTAT+k])){
DEBUGPRINT("- ERROR, it %li phi[%li] is INF, width %Le and dz %Le\n", it, LISTAT+k, width, dz);
getchar();
}
if(isnan(Phi[LISTAT+k])){
DEBUGPRINT("- ERROR, it %li phi[%li] is NAN, width %Le and dz %Le\n", it, LISTAT+k, width, dz);
getchar();
}
if((Phi[LISTAT+k] >= 1.0L)&&(k>LISTAT+1.0L/dz)){
DEBUGPRINT("- ERROR, it %li phi[%li] is %Le\n", it, LISTAT+k, Phi[LISTAT+k]);
getchar();
}*/
                } else {
                    fprintf(stderr, "ERROR - Invalid profile of particle %li, z: %li %Le, it:%li\n", i, k, ry[i]-radius[i], it);
                    exit(1);
                }
            }
        }
// Now, integrate packing fraction profile
        long double iPhidz = 0.0L, Phib = Phibase();
        for (k = LISTAT; k < NSTAT; k++)
            iPhidz += Phi[k]*dz;
// And discover Phi of base (where it is nearly constant)
        zbase = iPhidz/Phib-0.5L; // Removing 1 average grain radius from zbase will give a prediction of where phib/2 of the compact packing happen
        if (zbase + height > zb)
            H = ceil((zbase + height)/dz) +LISTAT;
        if (zbase + height < max_y){
            H = ceil((max_y + height)/dz) +LISTAT;
            fprintf(stderr, "WARNING - Height changed, since there are grains above this imposed height!\n");
            fprintf(out1, "#WARNING - Height changed, since there are grains above this imposed height!\n");
            height = max_y + height - zbase;
        }
        if (H > NSTAT){
            fprintf(stderr, "ERROR - Number of points to do meshes is lesser than stored in memory! H = %li NSTAT = %i\n", H, NSTAT);
//DEBUGPRINT("H %li iphidz %Le phib %Le height %Le\n", H, iPhidz, Phibase(), height);
            exit(1);
        }

        bottom=0.0L;
        for (i=nfg+1;i<=nfg+nwall;i++)
            bottom+=ry[i]+radius[i];
        bottom/=(double long)nwall;
        fprintf(out1, "#Grain diameter: %Le\n", 1.0L);
        fprintf(out1, "#Gravity: %Le\n", gravity);
        fprintf(out1, "#Density ratio (Grain/Fluid): %Le\n", density_ratio);
        fprintf(out1, "#Reynolds number of the fluid: %Le\n", reynolds);
        fprintf(out1, "#Shields number: %Le\n", shields);
        fprintf(out1, "#Fluid height: %Le\n", height);
        fprintf(out1, "#Viscosity: %Le\n", viscosity);
        fprintf(out1, "#Shear velocity: %Le\n", ustar);
        fprintf(out1, "#Height after bed: %Le\n", height);
        fprintf(out1, "#Base position: %Le\n#Bed position: %Le\n#Static Bed position: %Le\n#Turbulence effect: %Le\n", bottom, zb, zbase, z0 = 0.1L*viscosity/ustar);
        fflush(out1);
    }
#endif
    if (menu<15) {
        fprintf(out1,"#Pressure =%Le\n",totalweight*costheta/(double long)width);
        fprintf(out1,"#Shear stress =%Le\n",totalweight*sintheta/(double long)width);
    } else if (menu < 30) {
        bottom=0.0L;
        for (i=nfg+1;i<=nfg+nwall;i++)
            bottom+=ry[i];
        bottom/=(double long)nwall;
        height=0.0L;
        for (i=nfg+nwall+1;i<=ngtot;i++)
            height+=ry[i];
        height/=(double long)nwall;
        fprintf(out1,"#Bottom=%Le  Height=%Le\n",bottom,height);
        for (i=1;i<=nwall;i++) {
            ywall[i]=ry[nfg+nwall+i]-height;
            xwall[i]=rx[nfg+nwall+i];
            vx[nfg+nwall+i]=wallv;
            vpx[nfg+nwall+i]=wallv;
        }
        dhdt=0.0L;
        ddhdtt=0.0L;
        wallx=0.0L;
    }
#ifdef FLUID
    else if (menu < 45){
        fprintf(out1,"#Pressure cause by grains on the wall = %Le\n",totalweight/(double long)width);
        fflush(out1);
    }
#endif


// Opening the measurement file of the energy ("timeout.tsv")
    if ((out = fopen ("timeout.tsv","w")) == NULL){
        fprintf(stderr,"Probleme d'ouverture de %s\n","timeout.tsv");
        exit (1);
    }
#ifdef FLUID
    fprintf(out,"#Time\tEK\tER\tEE\tEG\tEF\t<u>\tphiB\tq\tn\tv\tncontot\tncontglis\tnsi\tnsi1\tviol_max\t<viol>\tfxdown\tfydown\n");
#else
    fprintf(out,"#Time\tEK\tER\tEE\tEG\tphiB\tq\tn\tv\tncontot\tncontglis\tnsi\tnsi1\tviol\t<viol>\tfxdown\tfydown\n");
#endif

/** Initialisation of statistical profiles*/
/*#ifdef FLUID
    z[i] = -3.0L;
    for (i = 1; i < H; i++){
        if (z[i-1] < zb-2.5L)
            z[i] = z[i-1] +0.3L;
        else if (z[i] < zb+5.0L)
            z[i] = z[i-1] +0.02L;
        else {
            z[i] = z[i-1] +1.0L;
        }
    }
#else
    for (i = 0 ;i < NSTAT;i++)
        z[i] = z[i-1] +0.1L;
#endif*/
    for (i = 0; i < H; i++){
        P[i] = 0.0L;
        Vx[i] = 0.0L;
        Vy[i] = 0.0L;
        Fx[i] = 0.0L;
        Fy[i] = 0.0L;
        Zcont[i] = 0.0L;
        Zglis[i] = 0.0L;
        interp[i] = 0.0L;
#ifdef FLUID
        fxFa[i] = 0.0L;
        fyFa[i] = 0.0L;
        fxFd[i] = 0.0L;
        fyFd[i] = 0.0L;
        u[i] = 0.0L;
        tau[i] = 0.0L;
        dpressuredz[i] = -fluid_density*gravity;
#endif
    }

#ifdef FLUID
    if ((menu == 35)||(menu == 36)){
    // First try to set the fluid's profile near to the stationary regime.
        long double h_imposedshear = 0.0L; // Position to impose the shear

        if (shields > 0.1L){
            h_imposedshear = zbase -3.0L*shields;
        } else {
            h_imposedshear = zbase;
        }

        long begin = ceil((h_imposedshear-0.5L)/dz);
        for (i = begin+LISTAT; i < H; i++){
            long double l = 0.0L;// kapa*(i-begin-LISTAT)*dz*(1.0L-exp(-((i-begin-LISTAT)*dz*ustar)/(RvD*viscosity)));
            u[i] = RK4((i-begin-LISTAT-1)*dz, u[i-1], l, 0.0L, &U); // Runge-Kutta solver of 4th order to set the fluid's equations near to the stationary regime
            tau[i] = ustar*ustar;
        }
        tau[H] = ustar*ustar;
    }

    if (menu >= 37)
        loadfluid("fluidin.bin");
#endif

/*********** Determination time step  ***********/
    dt = 0.02L/sqrt(kn); /* ! xx/50 -> arbitraire!*/
    fprintf(out1, "#dt = %Le\n", dt);

    if (menu==27)
        timewall=dt/timewall;

/********** Coefficients of the predictor-corrector algorithm  ***********/
    halfdt2 = 0.5L*dt*dt;//CHECK that this is indeed predictor-corrector
    halfdt = 0.5L*dt;

    fflush(out1);
}

void extractwall() {
    long i,j;
    double long rxlo,rylo,radiuslo,minlo,depth=10.0L;
    long malo;

    if ((menu>14) && (menu<30)) {
        height=0.0L;
        for (i=nfg+nwall+1;i<=ngtot;i++)
            height+=ry[i];
        height/=(double long)nwall;
        height-=depth;
    }

#ifdef FLUID
    if ((menu >= 30) && (menu < 45)){
        long double iPhidz = 0.0L;
        for (i = LISTAT; i <= H; i++)
            iPhidz += Phi[i]*dz;
        depth = iPhidz/Phibase()-height-0.5L;
        if (depth < 0.0L){
            fprintf(stderr, "ERROR - Impossible to extract such wall! Granular height: %Le Height to remove: %Le \n", iPhidz/Phibase(), height);
            exit(1);
        } else if (depth < height)
            fprintf(stderr, "WARNING - Remaining depth is smaller than removed! Granular height: %Le Height to remove: %Le \n", iPhidz/Phibase(), height);
    }
#endif

    i=0;
    nwall=0;
    nfg=ngtot;
    while (i<nfg) {
        i++;
        if (ry[i]<depth) {
            rxlo=rx[i];
            rylo=ry[i];
            radiuslo=radius[i];
            for (j=i;j<ngtot;j++) {
                rx[j]=rx[j+1];
                ry[j]=ry[j+1];
                radius[j]=radius[j+1];
            }
            if (rylo<depth-2.0L) {
                ngtot--;
            } else {
                nwall++;
                rx[ngtot]=rxlo;
                ry[ngtot]=rylo;
                radius[ngtot]=radiuslo;
            }
            i--;
            nfg=ngtot-nwall;
        }
    }

    if ((menu>14) && (menu<30)) {
        fprintf(out1,"#nfg tempo=%ld\n",nfg);
        for (i=1;i<=nwall;i++) {
            minlo=1.0E40L;
            for (j=1;j<=nfg-i+1;j++)
                if ((ry[j]>height)&&(ry[j]<minlo)) {
                    minlo=ry[j];
                    malo=j;
                }
            rxlo=rx[nfg-i+1];
            rylo=ry[nfg-i+1];
            radiuslo=radius[nfg-i+1];
            rx[nfg-i+1]=rx[nfg+i];
            ry[nfg-i+1]=ry[nfg+i];
            radius[nfg-i+1]=radius[nfg+i];
            rx[nfg+i]=rx[malo];
            ry[nfg+i]=ry[malo];
            radius[nfg+i]=radius[malo];
            rx[malo]=rxlo;
            ry[malo]=rylo;
            radius[malo]=radiuslo;
        }
        nfg=ngtot-2*nwall;
        i=0;
        while (i<nfg) {
            i++;
            if (ry[i]>height) {
                for (j=i;j<ngtot;j++) {
                    rx[j]=rx[j+1];
                    ry[j]=ry[j+1];
                    radius[j]=radius[j+1];
                }
                ngtot--;
                nfg--;
                i--;
            }
        }
    }
    for (i=1;i<=ngtot;i++) {
        ry[i]-= depth;
        vx[i] = 0.0L;
        vy[i] = 0.0L;
        ax[i] = 0.0L;
        ay[i] = 0.0L;
        phi[i] = 0.0L;
        ome[i] = 0.0L;
        domedt[i] = 0.0L;
    }
    return;
}


/*********** Prepare a new configuration  ***********/
void prepareconfig() {
    long ncpt;
    double long dx,xlo;
    long i,j;

    /*Pick up the radius at random according to a uniform law*/
    for (i=1;i<=ngtot;i++)
        radius[i] = 0.5L + deltar*(ran2(&islu)-0.5L);

    ncpt = 0;
    i=-(halfwidth/2)-1;
    j=1;
    for (ncpt=1;ncpt<=nfg;ncpt++) {
        i+=1;
        if (i>(halfwidth/2-1)) {
            i=-(halfwidth/2);
            j++;
        }
        rx[ncpt] =2.0L*(double long)i+0.3L*ran2(&islu)+(double long)(j%2);
        ry[ncpt] =2.0L*(double long)j+0.3L*ran2(&islu);
    }

    //The boundary is consituted by disks
    dx = width / (double long)(nwall);
    xlo=-halfwidth;
    for (i=1;i<=nwall;i++) {
        rx[nfg+i] = xlo+0.2*ran2(&islu);
        ry[nfg+i] = 0.2L*ran2(&islu);
        xlo = xlo + dx;
    }

    if (menu==15) {
        height=2.0L*(double long)(j+1.0L);
        xlo=-halfwidth;
        for (i=1;i<=nwall;i++) {
            rx[nfg+nwall+i] = xlo+0.2*ran2(&islu);
            ry[nfg+nwall+i] =height + 0.2L*ran2(&islu);
            xlo = xlo + dx;
        }
    }

    for (i=1;i<=ngtot;i++) {
        vx[i] = 0.0L;
        vy[i] = 0.0L;
        ax[i] = 0.0L;
        ay[i] = 0.0L;
        phi[i] = 0.0L;
        ome[i] = 0.0L;
        domedt[i] = 0.0L;
    }
    return;
}

/*********** Find neighbours in a crude way (in N^2)  ***********/
/*********** Require to know ior and iex ***********/
void findneighbours() {
    long icon; // index of contacts evoluting with neighbour index
    long ineigh;  // neighbour index
    long i,j;
    double long x1,x12,xi,yi, xij,yij,dij;

    x1 = 1.1L*(1.0L+deltar); // contacts must remain amongst neighbours when they open
    x12 = x1*x1;
    ior[ncontot+1] = 0;
    iex[ncontot+1] = 0;
    icon = 1;
    ineigh = 0;

    /*********** neighbours between free grains  ***********/
    for (i=1;i<=nfg;i++) {
        xi = rx[i];
        yi = ry[i];
        for (j=i+1;j<=nfg;j++) {
            xij = distperio(rx[j] - xi);
            yij = ry[j] - yi;
            dij = (xij*xij + yij*yij);
            if (dij < x12) {
                ineigh++;
                listi[ineigh] = i;
                listj[ineigh] = j;
                if ((ior[icon]==i)&&(iex[icon]==j)) { // We know in advance that the next one previously in contact must remain in the neighbourhood
                    io[ineigh] = 1; //Already in contact before
                    icon++;
                }
                else
                    io[ineigh] = 0;
            }
        }
    }
    nl = ineigh; // index reached at the end of free-free neighbourhood
    /*********** neighbours between  wall grains and free grains ***********/
    for (i=1;i<=nfg;i++) {
        xi = rx[i];
        yi = ry[i];
        for (j =nfg+1;j<=ngtot;j++) {
            xij = distperio(rx[j] - xi);
            yij = ry[j] - yi;
            dij = (xij*xij + yij*yij);
            if (dij < x12) {
                ineigh++;
                listi[ineigh] = i;
                listj[ineigh] = j;
                if ((ior[icon]==i)&&(iex[icon]==j)) {
                    io[ineigh] = 1;
                    icon ++;
                } else
                    io[ineigh] = 0;
            }
        }
    }
    nlt = ineigh;// index reached at the end of free-wall neighbourhood

    return;
}

/*********** Find neighbours acording to cells (in NlogN - details at Computer Simulations of Liquids [Allen & Tildesley, 1991], Computational Granular Dynamics [Pöschel, 2004], Simulação Paralela de Materiais Granulares [Martins, G. H. B. - 2015])  ***********/
/*********** Require to know ior and iex ***********/
void findneighbours2() {
    // In case of test, creating a system without free grains, this routine should be aborted.
    if(nfg < 1)
        return;

    long icon; // index of contacts evoluting with neighbour index
    long ineigh;  // neighbour index
    double long x1,x12, xij, yij, dij;

    x1 = 1.1L*(1.0L+deltar); // contacts must remain amongst neighbours when they open
    x12 = x1*x1;

    //Finding the boundary to build the cells:
    double long maxx = rx[1] +radius[1], maxy = ry[1] +radius[1], minx = rx[1] -radius[1], miny = ry[1] -radius[1];
    long d, e, f, g, h, i, j, k, x, y;
    for (i = 2; i <= nfg; i++) {
        maxx = maxx < rx[i] +radius[i] ? rx[i] +radius[i] : maxx; // Maximum in x
        maxy = maxy < ry[i] +radius[i] ? ry[i] +radius[i] : maxy; // Maximum in y
        minx = minx > rx[i] -radius[i] ? rx[i] -radius[i] : minx; // Minimum in x
        miny = miny > ry[i] -radius[i] ? ry[i] -radius[i] : miny; // Minimum in y
//        gradius = gradius < radius[i] ? radius[i] : gradius; // Greater radius
    }

    double long L = maxx -minx; // Size of the x component, composed only by free grains
    double long M = maxy -miny; // Size of the y component, composed only by free grains
    double long l = x12*2.0L, m = x12*2.0L; // Size of each cell, based on the greater diameter, and the neighbour search size of findneighbours
    long size = (long)(50.0L*x12);
    long Ll = (long)(L/l) +1, Mm = (long)(M/m) +1; // Number of cells in x = Ll and y = Mm
    if ((Ll <= 5) || (Mm <= 5)) { // Too small system, solving by Findneighbours1
        findneighbours();
        return;
    }
    long cells[Ll][Mm][size], ccell[Ll][Mm]; // If a huge dispersion happens, the last fixed number should be changed
    for (i = 0; i < Ll; i++) {
        for (j = 0; j < Mm; j++) {
            ccell[i][j] = 0; // Counting number of grains in each cell
            for (k = 0; k < size; k++) {
                cells[i][j][k] = 0;
            }
        }
    }

    for (i = 1; i <= nfg; i++) {
        x = (long) ((rx[i] -minx) / l);
        y = (long) ((ry[i] -miny) / m);
        // Just to be sure that it won't be outside of the boundaries, but by this definition, it shouldn't happen (USIING PERIODIC BOUNDARY CONDITIONS)
        x = (x < 1) ? 1 : x;
        x = (x >= Ll-1) ? Ll -2 : x;
        y = (y < 1) ? 1 : y;
        y = (y >= Mm-1) ? Mm -2 : y;
        cells[x][y][ccell[x][y]] = i;
        ccell[x][y]++;
/*if (ccell[x][y]>=size){
DEBUGPRINT("ERROR, it:%li FN2 x:%li y:%li ccell:%li i:%li nfg:%li\n", it, x, y, ccell[x][y], i, nfg);
DEBUGPRINT("maxx:%Le maxy:%Le minx:%Le miny:%Le\n", maxx, maxy, minx, miny);
}*/
    }
    // Done! Grains are distributed in the cells, now it's time to transform the cells in neighbours...

    ior[ncontot+1] = 0;
    iex[ncontot+1] = 0;
    icon = 1;
    ineigh = 0;

    // This will set boundary condition in x and y directions, but the rotine of detect contacts will cut it off, if itn't on the studied problem
    for (x = 0; x < Ll; x++) { // Getting the element on x position of the cell
        for (y = 0; y < Mm; y++) { // Getting the element on y positon of the cell
            for (k = -1; k < 2; k++) { // Looking to the neighbours cells in x direction
                for(h = 0; h < 2; h++) { // Looking to the neighbours cells in y direction
                    if (!((k == -1) && (h == 0))) { // Eliminating the repeated search
                        for (g = 0; g < ccell[x][y]; g++) { // Getting all the elements on the cell
                            i = cells[x][y][g]; // Element one of the neighbourhood
                            // Looking for boundary conditions and unnecessary searchs:
                            d = x+k;
                            d = d < 1 ? Ll -2 : d; // Condition to the left side, taking the last x position in cell
                            d = d >= Ll-1 ? 1 : d; // Condition to the right side, taking the first x position in cell
                            e = y+h;
                            e = e < 1 ? Mm -2 : e; // Condition to the bottom side, taking the higher y position in cell, as it is builded, it should never happen
                            e = e >= Mm-1 ? 1 : e; // Condition to the top side, taking the lower y position in cell
                            for (f = ((h == 0) && (k == 0)) ? g+1 : 0; f < ccell[d][e]; f++) { // Getting all non repeated elements of the neighbour and self cell
                                j = cells[d][e][f]; // Element two of the neighbourhood

                                // As it was when find an elegible neighbour:
                                xij = distperio(rx[j] - rx[i]);
                                yij = ry[j] - ry[i];
                                dij = (xij*xij + yij*yij);
                                if (dij < x12) {
                                    ineigh++;
                                    listi[ineigh] = i;
                                    listj[ineigh] = j;
                                    if ((ior[icon]==i)&&(iex[icon]==j)) {//We know in advance that the next one previously in contact must remain in the neighbourhood
                                        io[ineigh] = 1; //Already in contact before
                                        icon++;
                                    }
                                    else
                                        io[ineigh] = 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    nl = ineigh; // index reached at the end of free-free neighbourhood

    // ¡Stupid way to do this! It uses the whole walls as neigborhood, but for now I see no other solution now. In a generic code, this will change a little bit...
    // Searching in walls on x:
    for (x = 0; x < Ll; x++) {
        for (y = 0; y < 2; y++) {
            for (g = 0; g < ccell[x][y]; g++) {
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++) {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12) {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j)) {
                            io[ineigh] = 1;
                            icon ++;
                        } else
                            io[ineigh] = 0;
                    }
                }
            }
        }
        for (y = Mm-2; y < Mm; y++) {
            for (g = 0; g < ccell[x][y]; g++) {
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++) {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12) {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j)) {
                            io[ineigh] = 1;
                            icon ++;
                        } else
                            io[ineigh] = 0;
                    }
                }
            }
        }
    }
    // Searching in walls on y
    for (y = 2; y < Mm-2; y++) {
        for (x = 0; x < 2; x++) {
            for (g = 0; g < ccell[x][y]; g++) {
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++) {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12) {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j)) {
                            io[ineigh] = 1;
                            icon ++;
                        } else
                            io[ineigh] = 0;
                    }
                }
            }
        }
        for (x = Ll-2; x < Ll; x++) {
            for (g = 0; g < ccell[x][y]; g++) {
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++) {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12) {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j)) {
                            io[ineigh] = 1;
                            icon ++;
                        } else
                            io[ineigh] = 0;
                    }
                }
            }
        }
    }

    nlt = ineigh;// index reached at the end of free-wall neighbourhood

    return;
}

/*********** definition function distperio  ***********/
double long distperio(double long x) {
    if (x>halfwidth)
         return (x-width);
    else if (x<-halfwidth)
         return (x+width);
    else return x;
}

/*********** detect contacts ***********/
/*********** require to know reactsave ***********/
void detectcontacts() {
    long i,j,il;
    double long xij,yij,dij,radius2;

    /*********** Initialization ***********/
    ncontot = 0;
    icont0 = 0;
    for (i = 1;i<=ngtot; i++)
        nZcont[i] = 0;

    /*********** Contact between all grains ***********/
    for (il=1;il<=nlt;il++) { //loop over neighbouring grains
        i=listi[il]; // get back which i is in the neighbourhood of j>i
        j=listj[il];
        xij = distperio(rpx[j] - rpx[i]);/***! periodic boundary conditions along x***/
        yij = rpy[j] - rpy[i];
        dij = xij*xij + yij*yij;
        radius2 = (radius[i]+radius[j])*(radius[i]+radius[j]);
        if (io[il]==1) //flag containing the information that the neighbourhood il (=neighbours i and j labelled by il), was corresponding at the previous time step to grains in contact
            icont0++;/***! index of previous contacts to get back the tangential force of previous time-step**/
        if (dij<radius2) {
            nZcont[i]++;
            nZcont[j]++;
            ncontot++;
            ior[ncontot] = i;
            iex[ncontot] = j;
            dij = sqrt(dij);
            eij[ncontot] = dij - radius[i] - radius[j]; //negative interpenetration distance
            xnij[ncontot] = xij / dij; //coordinates of the vector pointing from i (called or) to j (called ex)
            ynij[ncontot] = yij / dij;

            /***get back previous tangential force, if contact was already there***/
            if (io[il]==1) //flag see before...
                react[ncontot] = reactsave[icont0];
            else
                react[ncontot] = 0.0L;
            io[il]=1;// This neighbourhood il (neighbouring grains i,j) now corresponds to a contact
        } else {
            io[il] = 0; // This neighbourhood il (neighbouring grains i,j) does NOT correspond to a contact
        }

        if (il==nl)
            ncontfree=ncontot;
    }
    return;
}

/*********** Force calculation ***********/
void forcecalculation() {
    long i,j,il;
    double long xn,yn,xt,yt;
    double long fn,ft,fx,fy,vijn,vijt,ftest;

    /*********** initialize*********/
    if (menu<15) {
        for (i =1;i<=nfg;i++) {
            fpx[i] = 4.0L*radius[i]*radius[i]*sintheta;
            fpy[i] = - 4.0L*radius[i]*radius[i]*costheta; // Gravity
            gam[i] = 0.0L;
        }
    } else if (menu==25) {
        for (i =1;i<=nfg;i++) {
            fpx[i] = 0.0L;
            gam[i] = 0.0L;
            if ((rpy[i]<bottom+6.0L)&&(rpy[i]>bottom+3.0L))
                fpy[i] = -4.0L*radius[i]*radius[i]*gravity*((bottom+6.0L)-rpy[i])*(rpy[i]-(bottom+3.0L));
            else if ((rpy[i]>height-6.0L)&&(rpy[i]<height-3.0L))
                fpy[i] = 4.0L*radius[i]*radius[i]*gravity*((height-6.0L)-rpy[i])*(rpy[i]-(height-3.0L));
            else fpy[i] = 0.0L;
        }
    } else if (menu==26) {
        for (i =1;i<=nfg;i++) {
            fpx[i] = 0.0L;
            fpy[i] = 4.0L*radius[i]*radius[i]*gravity;
            gam[i] = 0.0L;
        }
    }
    #ifdef FLUID
    else if ((menu >= 30) && (menu < 45)) {
        for (i =1;i<=nfg;i++){
            fpx[i] = 0.0L;
            fpy[i] = -gravity/ivmass[i]; // Gravity, with mass = 1;
            gam[i] = 0.0L;
        }
    }
    #endif
    else {
        for (i =1;i<=nfg;i++) {
            fpx[i] = 0.0L;
            fpy[i] = 0.0L;
            gam[i] = 0.0L;
        }
    }

    for (i =1;i<=ngtot;i++)
        sumfn[i]=0.0L;
    ncontglis = 0;
    fyup=0.0L;
    fxup=0.0L;
    fydown=0.0L;
    fxdown=0.0L;
    springxx=0.0L;
    springyy=0.0L;
    springxy=0.0L;

    /*********** contact between all grains ***********/
    for (il = 1;il<= ncontot;il++) {
        i = ior[il];
        j = iex[il];
        /********normal forces ***********/
        xn = xnij[il];
        yn = ynij[il];
        vijn = xn*(vpx[j]-vpx[i])+yn*(vpy[j]-vpy[i]);
        fn = - eij[il]*kn- vijn*gn;
        sumfn[i]+=fn;
        sumfn[j]+=fn;//CHECK
        reacn[il] = fn;
        /*********** tangential forces  ***********/
        xt = -ynij[il];
        yt = xnij[il];
        vijt = xt*(vpx[j]-vpx[i])+yt*(vpy[j]-vpy[i]) - omep[i]*radius[i] - omep[j]*radius[j];
        ft = react[il] - kt*vijt*dt;
        ftest = frott*fn;
        if (fabsl(ft)>=ftest) {
            if ((nZcont[i]!=1)&&(nZcont[j]!=1))
                ncontglis=ncontglis+1;
            if (ft>0)
                ft = ftest;
            else
                ft = -ftest;
        }

        react[il] = ft;

        /*******Force resultants**********/
        fx = fn*xn + ft*xt;
        fy = fn*yn + ft*yt;
        fpx[i] -=  fx;
        fpy[i] -=  fy;
        fpx[j] += fx;
        fpy[j] += fy;
        gam[i] -= ft*radius[i];
        gam[j] -= ft*radius[j];


        /*Contact between free grains and wall*/
        if (il>ncontfree) {
            if (j>(nfg+nwall)) {
                fxup+=fx;
                fyup+=fy;
                springxx+=kn*xn*xn;
                springyy+=kn*yn*yn;
                springxy+=kn*xn*yn;
            } else {
                fxdown+=fx;
                fydown+=fy;
            }
        }
    }

    // Adding the fluid force on each free grain
    #ifdef FLUID
        if ((menu >= 30) && (menu < 45)) {
            for (i = 1; i <= ngtot; i++) {
                fluidforce(i);
            }
        }
    #endif

    /********** Add wall friction ***********/
    /* friction opposed to v[i] */
        if (menu==2) {
            for (i =1;i<=nfg;i++) {
                ftest=-muwall*sumfn[i]/(radius[i]*sqrt(1.0E-14L+vpx[i]*vpx[i]+vpy[i]*vpy[i]));
                fpx[i] += vpx[i]*ftest;
                fpy[i] += vpy[i]*ftest;
            }
        }

    /* friction opposed to <v[i]> i.e. along -X */
        if (menu==3)
            for (i =1;i<=nfg;i++)
                fpx[i] -= muwall*sumfn[i]/radius[i];

    /* friction opposed to <v[i]> i.e. along -X and proportionnal to fit of P */
        if (menu==4)
            for (i =1;i<=nfg;i++)
                fpx[i] -= muwall*(2.0L*M_PI)*(aa+bb*rpy[i]);// using a first fit on syy/sphi from a previous run.(2.0L*M_PI)*(37.172L-1.0696L*rpy[i])

    fxdown/=(double long)width;
    fydown/=(double long)width;
    fxup/=(double long)width;
    fyup/=(double long)width;
    /********** Archivage des reactions tangentielles ***********/
    for ( il = 1;il<=ncontot;il++)
        reactsave[il] = react[il];

    return;
}

#ifdef FLUID
// Write the differential equation in function of du/dz: Called U for solving RK4th order of 1st ODE
long double U(long double z, long double u, long double l, long double ifx) {
    if (l > 0.0L)
        return ((-viscosity +sqrt(viscosity*viscosity +4.0L*l*l*(ustar*ustar-ifx)))/(2.0L*l*l));
    else
        return (ustar*ustar-ifx)/viscosity;
}

// Write the differential equation in function of dl/dz: Called L for solving RK4th order of 1st ODE
long double L(long double z, long double l, long double u, long double ifx) {
    return (kapa*(1.0L -exp(-sqrt(fabs(u*l)/(7.0L*viscosity)))));
}

// Runge-Kutta method for 4th order in 1st ODE
long double RK4(long double x, long double y, long double const1, long double const2, long double (*F)(long double, long double, long double, long double)) {
    long double f = 0.0;
    long double k1 = 0.0, k2 = 0.0, k3 = 0.0, k4 = 0.0;
    k1 = (*F)(x, y, const1, const2);
    k2 = (*F)(x +dz/2.0L, y +k1*dz/2.0L, const1, const2);
    k3 = (*F)(x +dz/2.0L, y +k2*dz/2.0L, const1, const2);
    k4 = (*F)(x +dz, y +k3*dz, const1, const2);
    f = y + (k1 +2.0*k2 +2.0*k3 +k4)*dz/6.0;
    return f;
}

// This rotine solves the partial differential equation for the fluid, based on the mesh over z and dt for viscous case with implicit form and using Thomas algoritm to solve linear system.
void solvefluidvelocity() {
    long j;

    // Upadte the fluid
    // Condition to the stability of this integration:
    // dt < (dz)^2/2, for the viscous case
    // dt < (dz)^2/(2*SystemSize), for the turbulent case
    double long up[H], alpha[H], beta[H], lambda = dt*viscosity/(dz*dz); // up is the predicted fluid velocity

    // Implementing Van Driest's equation:
    alpha[0] = 0.0L;
    beta[0] = 0.0L;
    for (j = 1; j < H-1; j++){
        alpha[j] = -lambda/(1.0L+2.0L*lambda+lambda*alpha[j-1]);
        beta[j] = (u[j]-dt*fxFd[j]*Phi[j]/(fluid_density*(1.0L-Phi[j]))+lambda*beta[j-1])/(1.0L+2.0L*lambda+lambda*alpha[j-1]);
/*if((isnan(alpha[j]))){
DEBUGPRINT("ERROR, it %li alpha[%li] is NAN\n", it, j);
getchar();
}
if((isinf(alpha[j]))){
DEBUGPRINT("ERROR, it %li alpha[%li] is INF\n", it, j);
getchar();
}
if((isnan(beta[j]))){
DEBUGPRINT("ERROR, it %li beta[%li] is NAN\n", it, j);
getchar();
}
if((isnan(beta[j]))){
DEBUGPRINT("ERROR, it %li beta[%li] is INF\n", it, j);
getchar();
}*/
    }
    up[H-1] = u[H-1]+dt*(ustar*ustar-viscosity*(u[H-2]-u[H-3])/dz)/dz; // Imposed shear
//    up[H-1] = ustar*ustar*(height-dz)/viscosity; // Imposed velocity
//    long left = 0;
    for (j = H-2; j >= 0; j--){
        up[j] = beta[j]-alpha[j]*up[j+1];

// Personally I desagree to force system to go to zero if we have very small velocity, but as it is been a such small velocity, why not?
//        if (up[j] < 0.0L){
//            left = j;
//            break;
//            up[j] = 0.0L;
//        }
/*if((isnan(u[j]))){
DEBUGPRINT("ERROR, it %li u[%li] is NAN\n", it, j);
getchar();
}
if((isinf(u[j]))){
DEBUGPRINT("ERROR, it %li u[%li] is INF\n", it, j);
getchar();
}
if((isnan(up[j]))){
DEBUGPRINT("ERROR, it %li up[%li] is NAN\n", it, j);
getchar();
}
if((isinf(up[j]))){
DEBUGPRINT("ERROR, it %li up[%li] is INF, tau+1 %Le, tau %Le, fxF %Le, zb %Le\n", it, j, tau[j], tau[j], fxFd[j], zb);
getchar();
}*/
    }
//    for (j = left; j >= 0; j++){
//        up [j] = 0.0L;
//    }

    u[0] = up[0];
    tau[0] = viscosity*(u[1]-u[0])/dz;
    dpressuredz[0] = -fluid_density*gravity;
    for (j = 1; j < H-1; j++){
        u[j] = up[j];
        tau[j] = viscosity*(up[j+1]-up[j-1])/(2.0L*dz);
        dpressuredz[j] = -fluid_density*gravity -fyFd[j]*Phi[j]/((1.0L-Phi[j]));
    }
    u[H-1] = up[H-1];
    tau[H-1] = ustar*ustar;
    dpressuredz[H-1] = -fluid_density*gravity;

    return;
}

// This rotine solve the partial differential equation for the fluid explicitly, based on the mesh over z and the dt using the packing fraction for the mixing length
void solvefluidvelocity2(){
    long j;
    double long iphi[H], maxphi;
    // Definition of iphi = integral of phi(pakingfraction) along z:
    iphi[0] = min(Phi[0], 0.64L);
    for (j = 1; j < H; j++){
        maxphi = min(max(Phi[j], maxphi), 0.64L);
    }
    iphi[0] = min((maxphi-Phi[0])*dz/2.0L, 0.64L)/maxphi;
    for (j = 1; j < H-1; j++){
        iphi[j] = min(iphi[j-1]+(maxphi-Phi[j])*dz, 0.64L)/maxphi;
    }
    iphi[H-1] = min(iphi[H-2]+(maxphi-Phi[H-1])*dz/2.0L, 0.64L)/maxphi;

    // Condition to the stability of this integration:
    // dt < (dz)^2/2, for the viscous case
    // dt < (dz)^2/(2*SystemSize), for the turbulent case
    double long up[H]; // l is the mixing lenght, tau is the shear per density of the fluid

    // Implementing Van Driest's equation:
    tau[0] = viscosity*u[0]/dz;
    for (j = 0; j < H-1; j++){
        double long l = kapa*iphi[j]*(1.0L -exp(-iphi[j]*ustar/(RvD*viscosity)));
        tau[j+1] = (viscosity+ l*l*fabs(u[j+1] -u[j])/dz)*(u[j+1] -u[j])/(dz);
        up[j] = u[j] +(tau[j+1]-tau[j])*dt/dz -dt*fxFd[j]/(fluid_density*(1.0L-Phi[j]));
        dpressuredz[j+1] = -fluid_density*gravity -fyFd[j]/((1.0L-Phi[j]));
/*if(Phi[j] >= 1.0L){
    if(j < LISTAT +1.0L/dz){
        up[j] = 0.0L;
        tau[j] = ustar*ustar;
    } else {
        DEBUGPRINT("ERROR, it %li phi[%li] is %Le\n", it, j, Phi[j]);
        getchar();
    }
}
if((isnan(u[j]))){
DEBUGPRINT("ERROR, it %li u[%li] is NAN\n", it, j);
getchar();
}
if((isinf(u[j]))){
DEBUGPRINT("ERROR, it %li u[%li] is INF\n", it, j);
getchar();
}
if((isnan(up[j]))){
DEBUGPRINT("ERROR, it %li up[%li] is NAN\n", it, j);
getchar();
}
if((isinf(up[j]))){
DEBUGPRINT("ERROR, it %li up[%li] is INF, tau+1 %Le, tau %Le, fxF %Le, zb %Le\n", it, j, tau[j], tau[j], fxFd[j], zb);
getchar();
}*/
    }
    tau[H] = ustar*ustar;
    up[H-1] = u[H-1] +(tau[H]-tau[H-1])*dt/dz;
    dpressuredz[H] = -fluid_density*gravity -fyFd[H-1]/(1.0L-Phi[H-1]);

    for (j =  0; j < H; j++)
        u[j] = up[j];

    return;
}

void solvefluidvelocity3(){
    long j;

    // Condition to the stability of this integration:
    // dt < (dz)^2/2, for the viscous case
    // dt < (dz)^2/(2*SystemSize), for the turbulent case
    double long up[H], l[H]; // l is the mixing lenght, tau is the shear per density of the fluid

    // Implementing Orencio's equation:
    l[0] = 0.0L;
    tau[0] = viscosity*u[0]/dz;
    for (j = 0; j < H-1; j++){
        l[j+1] = RK4((j-LISTAT)*dz, l[j], u[j], 0.0L, &L);
        tau[j+1] = (viscosity+ l[j]*l[j]*fabs(u[j+1] -u[j])/dz)*(u[j+1] -u[j])/(dz);
        up[j] = u[j] +(tau[j+1]-tau[j])*dt/dz -dt*fxFd[j]/(fluid_density*(1.0L-Phi[j]));
        dpressuredz[j+1] = -fluid_density*gravity -fyFd[j]/((1.0L-Phi[j]));
/*if(Phi[j] >= 1.0L){
    if(j < LISTAT +1.0L/dz){
        up[j] = 0.0L;
        tau[j] = ustar*ustar;
    } else {
        DEBUGPRINT("ERROR, it %li phi[%li] is %Le\n", it, j, Phi[j]);
        getchar();
    }
}
if((isnan(u[j]))){
DEBUGPRINT("ERROR, it %li u[%li] is NAN\n", it, j);
getchar();
}
if((isinf(u[j]))){
DEBUGPRINT("ERROR, it %li u[%li] is INF\n", it, j);
getchar();
}
if((isnan(up[j]))){
DEBUGPRINT("ERROR, it %li up[%li] is NAN\n", it, j);
getchar();
}
if((isinf(up[j]))){
DEBUGPRINT("ERROR, it %li up[%li] is INF, tau+1 %Le, tau %Le, fxF %Le, zb %Le\n", it, j, tau[j], tau[j], fxFd[j], zb);
getchar();
}*/
    }
    tau[H] = ustar*ustar;
    up[H-1] = u[H-1] +(tau[H]-tau[H-1])*dt/dz;
    dpressuredz[H] = -fluid_density*gravity -fyFd[H-1]/(1.0L-Phi[H-1]);

    for (j =  0; j < H; j++)
        u[j] = up[j];

    return;
}

// The implementation of the fluid force on the grain. This rotine should be used for each grain that have a fluid interaction. One important parameter to this calcul is the zb.
void fluidforce(long i){
    // This function have 2 diferent origins: a drag component and an Archimedes component

    // Calculating the fluid velocity
    double long fluidvelocity = 0.0L; //u (velocity of the fluid) is function of z component, but it also can be function of x (or y). Drag coefficient, function of Granular Reynolds' number.

    long index_top = LISTAT+ceil(rpy[i]/dz), index_bottom = LISTAT+floor(rpy[i]/dz);
    fluidvelocity = u[index_bottom]+(u[index_top]-u[index_bottom])*(rpy[i]-(index_bottom-LISTAT)*dz)/dz;
//    vmodulus = sqrt((fluidvelocity-vx[i])*(fluidvelocity-vx[i]) + (vy[i])*(vy[i])); // Modulus of the diference of velocity between grain and fluid
    fxfluid_drag[i] = PI*radius[i]*reynoldsc*viscosity*fluid_density*(fluidvelocity-vpx[i])/4.0L; // Cosine of the vector
    fyfluid_drag[i] = PI*radius[i]*reynoldsc*viscosity*fluid_density*(-vpy[i])/4.0L; // Sine of the vector
    /// Archimedes force
    // Archimedes force has 2 components: against gravity and in the shearing direction
    // Calculating the hydrodynamic force
    fxfluid_arch[i] = fluid_density*(tau[index_top] -tau[index_bottom])/dz; // Shear from fluid
    fyfluid_arch[i] = -(dpressuredz[index_bottom]+(dpressuredz[index_top]-dpressuredz[index_bottom])*(rpy[i]-(index_bottom-LISTAT)*dz)/dz);
    fxfluid_arch[i] *= (4.0L*PI*radius[i]*radius[i]*radius[i])/3.0L;
    fyfluid_arch[i] *= (4.0L*PI*radius[i]*radius[i]*radius[i])/3.0L;
    // Calculating the hydrostatic force

    // An important hint: The fluid it self have a time constant to ajust the drag force. If the fluid time is greater than the dt, divernces may happen.
    // Force caused by fluid on grains:
    fpx[i] += fxfluid_drag[i] +fxfluid_arch[i];
    fpy[i] += fyfluid_drag[i] +fyfluid_arch[i];

    fxfluid_drag[i] /= (4.0L*PI*radius[i]*radius[i]*radius[i])/3.0L;
    fyfluid_drag[i] /= (4.0L*PI*radius[i]*radius[i]*radius[i])/3.0L;
    fxfluid_arch[i] /= (4.0L*PI*radius[i]*radius[i]*radius[i])/3.0L;
    fyfluid_arch[i] /= (4.0L*PI*radius[i]*radius[i]*radius[i])/3.0L;

    return;
}

void sheargrainswithfluid(){
    long i;
    if (it < nbiteration/2)
    // Increase the fluid velocity
        for (i = -LISTAT; i < H-LISTAT; i++)
            if (i > 0)
                u[LISTAT+i] = (i*dz)*it/nbiteration;
            else
                u[LISTAT+i] = 0.0L;
    else
    // Deslocate the fluid profile
        for (i = -LISTAT; i < H-LISTAT; i++){
            u[LISTAT+i] = 0.5L*((i*dz) -(PI*nfg/(6.0L*width*0.50L))*((it*2.0L/nbiteration)-1.0L));
            if(u[LISTAT+i] < 0.0L)
                u[LISTAT+i] = 0.0L;
        }
    return;
}
#endif

/*********** Predictor (require ax, ay, domedt)  ***********/
void predictor(){
    long i;

    if (menu==15){
        hp=height+dt*dhdt+halfdt2*ddhdtt;
        dhdtp=dhdt+dt*ddhdtt;
        ddhdttp=ddhdtt;
        wallx+=dt*wallv;
        if (wallx>width)
            wallx-=width;
        for (i=1;i<=nwall;i++){
            rpx[nfg+nwall+i]=distperio(xwall[i]+wallx);
            rpy[nfg+nwall+i]=ywall[i]+hp;
            vpy[nfg+nwall+i]=dhdtp;
        }
    } else if ((menu>15)&&(menu<29)) {
        if (menu==27){
            wallv=wallvzero*exp(-(double long)it*timewall);
            for (i=1;i<=nwall;i++)
                vpx[nfg+nwall+i]=wallv;
        }
        wallx+=dt*wallv;
        if (wallx>width)
            wallx-=width;
        for (i=1;i<=nwall;i++){
            rpx[nfg+nwall+i]=distperio(xwall[i]+wallx);
            rpy[nfg+nwall+i]=ywall[i]+height;//useless?
        }
    } else if (menu==29) {
        hp=height+dt*dhdt+halfdt2*ddhdtt;
        dhdtp=dhdt+dt*ddhdtt;
        ddhdttp=ddhdtt;
        wallxp=wallx+dt*wallv+halfdt2*walla;
        if (wallxp>0.5L*width)
            wallxp-=width;
        wallvp=wallv+dt*walla;
        wallap=walla;
        for (i=1;i<=nwall;i++){
            rpx[nfg+nwall+i]=distperio(xwall[i]+wallxp);
            rpy[nfg+nwall+i]=ywall[i]+hp;
            vpx[nfg+nwall+i]=wallvp;
            vpy[nfg+nwall+i]=dhdtp;
        }
    }

    /*********** Predictions sur les grains libres ***********/
    for (i=1;i<=nfg;i++){
        rpx[i] = distperio(rx[i] + dt*vx[i] + halfdt2*ax[i]);
        rpy[i] = ry[i] + dt*vy[i] + halfdt2*ay[i];
        phip[i] = phi[i] + dt*ome[i] + halfdt2*domedt[i];
        vpx[i] = vx[i] + dt*ax[i] ;
        vpy[i] = vy[i] + dt*ay[i] ;
        omep[i] = ome[i] + dt*domedt[i] ;
        apx[i] = ax[i] ;
        apy[i] = ay[i] ;
        domedtp[i] = domedt[i];
    }
    return;
}

/*********** Corrector ***********/
void corrector(){
    long i;

    if (menu==15){
        ddhdtt=fyup-1.0L-gn*dhdtp;
        dhdt = dhdtp + halfdt*(ddhdtt-ddhdttp);
        height = hp;
        for (i=1;i<=nwall;i++){
            rx[nfg+nwall+i]=distperio(xwall[i]+wallx);
            ry[nfg+nwall+i]=ywall[i]+height;
        }
    } else if ((menu>15)&&(menu<29)) {
        if (menu<28)
            height+=(fyup-1.0L)/springyy;
        for (i=1;i<=nwall;i++){
            rx[nfg+nwall+i]=distperio(xwall[i]+wallx);
            ry[nfg+nwall+i]=ywall[i]+height;
        }
    } else if (menu==29) {
        ddhdtt=fyup-1.0L-gn*dhdtp;
        dhdt = dhdtp + halfdt*(ddhdtt-ddhdttp);
        height = hp;
        walla=shear+fxup;
        wallv = wallvp + halfdt*(walla-wallap);
        wallx = wallxp;
        for (i=1;i<=nwall;i++){
            rx[nfg+nwall+i]=distperio(xwall[i]+wallx);
            ry[nfg+nwall+i]=ywall[i]+height;
        }
    }

/*********** Correction on free grains ***********/
    for (i = 1;i<=nfg;i++){
        ax[i] = fpx[i]*ivmass[i];
        ay[i] = fpy[i]*ivmass[i];
        domedt[i] = gam[i]*ivmomine[i];
        vx[i] = vpx[i] + halfdt*(ax[i]-apx[i]);
        vy[i] = vpy[i] + halfdt*(ay[i]-apy[i]);
        ome[i] = omep[i] + halfdt*(domedt[i]-domedtp[i]);
        rx[i] = rpx[i];
        ry[i] = rpy[i];
        phi[i] = phip[i];
    }
    return;
}

long double volumefraction(long i, long double lower, long double upper){
#ifndef FLUID
// 2D case
    if ((radius[i] > fabsl(rpy[i]-upper)) && (radius[i] > fabsl(rpy[i]-lower)))
        return (sqrt(radius[i]*radius[i]-(rpy[i]-upper)*(rpy[i]-upper))*(upper-rpy[i])-sqrt(radius[i]*radius[i]-(rpy[i]-lower)*(rpy[i]-lower))*(lower-rpy[i])+radius[i]*radius[i]*(asin((upper-rpy[i])/radius[i])-asin((lower-rpy[i])/radius[i]));
    else if (radius[i] > fabsl(rpy[i]-upper))
        return sqrt(radius[i]*radius[i]-(rpy[i]-upper)*(rpy[i]-upper))*(upper-rpy[i])+radius[i]*radius[i]*(asin((upper-rpy[i])/radius[i])+M_PI*0.5L);
    else if (radius[i] > fabsl(rpy[i]-lower))
        return (-sqrt(radius[i]*radius[i]-(rpy[i]-lower)*(rpy[i]-lower))*(lower-rpy[i]))-radius[i]*radius[i]*(asin((lower-rpy[i])/radius[i])-M_PI*0.5L);
    else
        return PI*radius[i]*radius[i];
#else
// 3D case (formulation extracted from http://mathworld.wolfram.com/SphericalSegment.html)
/*    if (rpy[i] -bottom > radius[i])
        bottom = rpy[i]-radius[i];
    if (upper -rpy[i] > radius[i])
        upper = rpy[i]+radius[i];
    long double h = upper-bottom;
    long double d = bottom - rpy[i];
    return PI*h*(radius[i]*radius[i]-d*d-d*h-h*h/3.0L);*/

    double long vs = 4.0L*PI/3.0L*radius[i]*radius[i]*radius[i],
                ha = radius[i]-fabsl(rpy[i]-lower),
                a = sqrt(radius[i]*radius[i]-(rpy[i]-lower)*(rpy[i]-lower)),
                hb = radius[i]-fabsl(rpy[i]-upper),
                b = sqrt(radius[i]*radius[i]-(rpy[i]-upper)*(rpy[i]-upper));
    if ((radius[i] > fabs(rpy[i]-upper)) && (radius[i] > fabs(rpy[i]-lower)))
        if ((lower < rpy[i]) && (rpy[i] < upper))
            return vs - PI*((3.0L*a*a+ha*ha)*ha+(3.0L*b*b+hb*hb)*hb)/6.0L;
        else
            return PI*fabsl((3.0L*a*a+ha*ha)*ha-(3.0L*b*b+hb*hb)*hb)/6.0L;
    else if (radius[i] > fabsl(rpy[i]-upper))
        if (rpy[i] > upper)
            return PI*(3.0L*b*b+hb*hb)*hb/6.0L;
        else
            return vs - PI*(3.0L*b*b+hb*hb)*hb/6.0L;
    else if (radius[i] > fabsl(rpy[i]-lower))
        if (rpy[i] < lower)
            return PI*(3.0L*a*a+ha*ha)*ha/6.0L;
        else
            return vs - PI*(3.0L*a*a+ha*ha)*ha/6.0L;
    else
        return vs;
#endif
}

long double areafraction(long i, long double lower, long double upper){
#ifdef FLUID
// 3D case (formulation extracted from http://mathworld.wolfram.com/SphericalSegment.html)
/*    if (rpy[i] -bottom > radius[i])
        bottom = rpy[i]-radius[i];
    if (upper -rpy[i] > radius[i])
        upper = rpy[i]+radius[i];
    long double h = upper-bottom;
    return 2.0L*PI*radius[i]*h;*/

    double long as = 4.0L*PI*radius[i]*radius[i],
                ha = radius[i]-fabsl(rpy[i]-lower),
                hb = radius[i]-fabsl(rpy[i]-upper);
    if ((radius[i] > fabs(rpy[i]-upper)) && (radius[i] > fabs(rpy[i]-lower)))
        if ((lower < rpy[i]) && (rpy[i] < upper))
            return as -2.0L*PI*radius[i]*(ha+hb);
        else
            return fabsl(2.0L*PI*radius[i]*(ha-hb));
    else if (radius[i] > fabsl(rpy[i]-upper))
        if (rpy[i] > upper)
            return 2.0L*PI*radius[i]*hb;
        else
            return as - 2.0L*PI*radius[i]*hb;
    else if (radius[i] > fabsl(rpy[i]-lower))
        if (rpy[i] < lower)
            return 2.0L*PI*radius[i]*ha;
        else
            return as - 2.0L*PI*radius[i]*ha;
    else
        return as;
#endif
}

void fraction(long i, long double lower, long double upper, long double &volume, long double &area){
// 3D case (formulation extracted from http://mathworld.wolfram.com/SphericalSegment.html)
/*    if (rpy[i] -bottom > radius[i])
        bottom = rpy[i]-radius[i];
    if (upper -rpy[i] > radius[i])
        upper = rpy[i]+radius[i];
    long double h = upper-bottom;
    long double d = bottom - rpy[i];
    volume = PI*h*(radius[i]*radius[i]-d*d-d*h-h*h/3.0L);
    area = 2.0L*PI*radius[i]*h;*/

    double long vs = 4.0L*PI/3.0L*radius[i]*radius[i]*radius[i],
                as = 4.0L*PI*radius[i]*radius[i],
                ha = radius[i]-fabsl(rpy[i]-lower),
                a = sqrt(radius[i]*radius[i]-(rpy[i]-lower)*(rpy[i]-lower)),
                hb = radius[i]-fabsl(rpy[i]-upper),
                b = sqrt(radius[i]*radius[i]-(rpy[i]-upper)*(rpy[i]-upper));
    if ((radius[i] > fabs(rpy[i]-upper)) && (radius[i] > fabs(rpy[i]-lower)))
        if ((lower < rpy[i]) && (rpy[i] < upper)){
            volume = vs - PI*((3.0L*a*a+ha*ha)*ha+(3.0L*b*b+hb*hb)*hb)/6.0L;
            area = as -2.0L*PI*radius[i]*(ha+hb);
            return;
        } else {
            volume = PI*fabsl((3.0L*a*a+ha*ha)*ha-(3.0L*b*b+hb*hb)*hb)/6.0L;
            area = fabsl(2.0L*PI*radius[i]*(ha-hb));
            return;
        }
    else if (radius[i] > fabsl(rpy[i]-upper))
        if (rpy[i] > upper){
            volume = PI*(3.0L*b*b+hb*hb)*hb/6.0L;
            area = 2.0L*PI*radius[i]*hb;
            return;
        } else {
            volume = vs - PI*(3.0L*b*b+hb*hb)*hb/6.0L;
            area = as - 2.0L*PI*radius[i]*hb;
            return;
        }
    else if (radius[i] > fabsl(rpy[i]-lower))
        if (rpy[i] < lower){
            volume = PI*(3.0L*a*a+ha*ha)*ha/6.0L;
            area = 2.0L*PI*radius[i]*ha;
            return;
        } else {
            volume = vs - PI*(3.0L*a*a+ha*ha)*ha/6.0L;
            area = as - 2.0L*PI*radius[i]*ha;
            return;
        }
    else {
        volume = vs;
        area = as;
        return;
    }
    return;
}

#ifdef FLUID
void calculateprofiles(){
    long i,k;

/*********** initialize*********/
    for(k = 0; k < H; k++){
        Phi[k] = 0.0L;
        fxFd[k] = 0.0L;
        fyFd[k] = 0.0L;
    }

    for (i=1;i<=ngtot;i++){
        long begin = floor((rpy[i]-radius[i])/dz);
        long end = ceil((rpy[i]+radius[i])/dz);
        for (k = begin; k < end; k++){
            long double volumeratio = 0.0L;
            volumeratio = volumefraction(i, k*dz, (k+1)*dz)/(width*dz);

            // To be sure that k is not invalidating any data:
            if ((k >= -LISTAT) && (k < H-LISTAT)){
                Phi[LISTAT+k] += volumeratio;
                //From the proportion of drag force felt by grain in the slice to shear
                fxFd[LISTAT+k] += volumeratio*fxfluid_drag[i];
                //From the proportion of drag force felt by grain in the slice to pressure
                fyFd[LISTAT+k] += volumeratio*fyfluid_drag[i];
            } else {
                fprintf(stderr, "ERROR - Invalid profile of particle %li, z: %li %Le, it:%li\n", i, k, rpy[i]-radius[i], it);
                exit(1);
            }
        }
    }
}
#endif

/*********** Correction on free grains ***********/
void cumulatestats(){
    long i,j,k,il;
    double long fx,fy,interpgrain[ngtot+1];  // weight function
// bin size =0.1d;

/*********** initialize*********/
    for (i = 1; i<= nfg; i++)
        nZglis[i] = 0;
    for (i = 1;i <= ngtot; i++)
        interpgrain[i] = 0.0L;

    for(k = 0; k < H; k++){
        Phi[k] = 0.0L;
        Vx[k] = 0.0L;
        Vy[k] = 0.0L;
        P[k] = 0.0L;
        Zcont[k] = 0.0L;
        Zglis[k] = 0.0L;
        interp[k] = 0.0L;
        Fx[k] = 0.0L;
        Fy[k] = 0.0L;
#ifdef FLUID
        fxFa[k] = 0.0L;
        fyFa[k] = 0.0L;
        fxFd[k] = 0.0L;
        fyFd[k] = 0.0L;
#endif
    }

/*********** contact between free grains ***********/
    for (il = 1;il<= ncontfree;il++){
        interpgrain[ior[il]] += -eij[ior[il]];
        interpgrain[iex[il]] += -eij[iex[il]];
        if (fabsl(react[il])>=(0.99L*reacn[il])){
            nZglis[ior[il]]+=1;
            nZglis[iex[il]]+=1;
        }
    }
    for (il=ncontfree+1;il<=ncontot;il++){
        interpgrain[ior[il]] += -eij[ior[il]];
        if (fabsl(react[il])>=(0.99L*reacn[il])){
            nZglis[ior[il]]+=1;
        }
    }

    for( il=1;il<=ncontot;il++){
        fx = reacn[il]*xnij[il] - react[il]*ynij[il];
        fy = reacn[il]*ynij[il] + react[il]*xnij[il];
        i = ior[il];
        j = iex[il];
        if (rpy[i]>rpy[j]){
            long begin = ceil(rpy[j]/dz);
            long end = floor(rpy[i]/dz);
            for (k = begin; k <= end; k++){
                Fx[LISTAT+k] += fx;
                Fy[LISTAT+k] -= fy;
            }
        } else {
            long begin = ceil(rpy[i]/dz);
            long end = floor(rpy[j]/dz);
            for (k = begin; k <= end; k++){
                Fx[LISTAT+k] -= fx;
                Fy[LISTAT+k] += fy;
            }
        }
    }
    for (i=1;i<=ngtot;i++){
        long begin = floor((rpy[i]-radius[i])/dz);
        long end = ceil((rpy[i]+radius[i])/dz);
        for (k = begin; k < end; k++){
#ifndef FLUID
// 2D case
            ratio = volumefraction(i, k*dz, (k+1)*dz)/(width*dz);
            Phi[LISTAT+k] += ratio;
            Vx[LISTAT+k] += ratio*vx[i];
            Vy[LISTAT+k] += ratio*vy[i];
            P[LISTAT+k]+=ratio*sumfn[i]/(2.0L*M_PI*radius[i]);
            Zcont[LISTAT+k] += ratio*nZcont[i];
            Zglis[LISTAT+k] += ratio*nZglis[i];
            interp[LISTAT+k] += ratio*interpgrain[i];
#else
// 3D case
            long double volumeratio = 0.0L, arearatio = 0.0L;
            fraction(i, k*dz, (k+1)*dz, volumeratio, arearatio);
            volumeratio /= width*dz;
            arearatio /= width*dz;

            // To be sure that k is not invalidating any data:
            if ((k >= -LISTAT) && (k < H-LISTAT)){
                Phi[LISTAT+k] += volumeratio;
                Vx[LISTAT+k] += volumeratio*vpx[i];
                Vy[LISTAT+k] += volumeratio*vpy[i];
                P[LISTAT+k] += volumeratio*sumfn[i]/arearatio;
                Zcont[LISTAT+k] += volumeratio*nZcont[i];
                Zglis[LISTAT+k] += volumeratio*nZglis[i];
                interp[LISTAT+k] += volumeratio*interpgrain[i];
                //From the proportion of drag force felt by grain in the slice to shear
                fxFd[LISTAT+k] += volumeratio*fxfluid_drag[i];
                //From the proportion of drag force felt by grain in the slice to pressure
                fyFd[LISTAT+k] += volumeratio*fyfluid_drag[i];
                //From the proportion of Archimedes force felt by grain in the slice to shear
                fxFa[LISTAT+k] += volumeratio*fxfluid_arch[i];
                //From the proportion of Archimedes force felt by grain in the slice to pressure
                fyFa[LISTAT+k] += volumeratio*fyfluid_arch[i];
            } else {
                fprintf(stderr, "ERROR - Invalid profile of particle %li, z: %li %Le, it:%li\n", i, k, rpy[i]-radius[i], it);
                exit(1);
            }
#endif
        }
    }
}

long double Phibase(){
    long double Phimax = Phi[0];
    long i, phibi = 0, phibf = 0;
// Find Phimax in the profile
    for (i = 1; i < H; i++)
        Phimax = max(Phimax, Phi[i]);
// Find where Phimax/2 is from bottom to top
    for (i = 0; i < H; i++){
        if (Phi[i] > Phimax/2.0L) {
            phibi = i;
            break;
        }
    }
// Find where Phimax/2 is from top to bottom
    for (i = H-1; i > 0; i--) {
        if (Phi[i] > Phimax/2.0L){
            phibf = i;
            break;
        }
    }
// Integrate Phi between Phi initial and Phi final
    long double phib = Phi[phibi]/2.0L;
    for (i = phibi+1; i < phibf; i++)
        phib += (Phi[i]);
    phib += Phi[phibf]/2.0L;
// Average the integral to get Phib
    phib /= (phibf-phibi);
    return phib;
}

/*********** timemeasurement ***********/
void timemeasurement(){
    long i, il;
    /*********** Test des violations ( = profondeurs de penetration) ***********/
    /*********** Maximal interpenetration distance ***********/
    long double viol = 0.0L, violmean = 0.0L;
    for(il = 1; il <= ncontot; il++){
        violmean += eij[il];
        if (eij[il]<viol)
            viol=eij[il];
    }
    viol = -0.5L*viol;
    violmean /= -0.5L*ncontot;

    /*********** Energies cinetique et elastictielle totales ***********/
    long double cineti = 0.0L, rotational = 0.0L;
    for (i = 1; i <= nfg; i++){
        cineti += 0.5L*(vx[i]*vx[i] + vy[i]*vy[i])/ivmass[i];
        rotational += 0.5L*(ome[i]*ome[i])/ivmomine[i];
    }
    long double  elastic = 0.0L;
    for (il = 1;il<=ncontot;il++)
        elastic += 0.5L*kn*eij[il]*eij[il];
    long double gravitational = 0.0L;
    for (i = 1; i <= nfg; i++)
#ifndef FLUID
        gravitational += gravity*ry[i]/ivmass[i];
#else
        gravitational += gravity*ry[i]/(ivmass[i]*(density_ratio-1.0L));
    long double fluid = 0.0L, ufluid = 0.0L;
    for (i = 0; i < H; i++){
        ufluid += u[i];
        fluid += u[i]*u[i];
    }
    fluid *= (1.0L-Phi[i])*fluid_density*width*dz/2.0L;
    ufluid /= H-LISTAT;

    /********************* Transport Parameters ************************/
    long double q = 0.0L, n = 0.0L, v = 0.0L, vps = 0.0L, vp = 0.0L, phib = Phibase(), zbar = 0.0L, lambda = 0.0L;
    for (i = 1; i <= nfg; i++){
        q += radius[i]*radius[i]*radius[i]*vpx[i];
        vps += vpx[i]*vpx[i];
        vp += vpx[i];
    }
    q *= 4.0L*PI/(3.0L*phib*width);
    if (vps != 0.0L)
        n = vp*vp/(width*vps);
    if (vp != 0.0L)
        v = vps/vp;
    if (q != 0.0L){
        for (i = LISTAT; i < H; i++)
            zbar += (i-LISTAT)*Vx[i];
        zbar *= dz*dz/(q*phib);
        for (i = LISTAT; i < H; i++)
            lambda += ((i-LISTAT)*dz-zbar)*((i-LISTAT)*dz-zbar)*Vx[i];
        lambda *= dz/q;
        if (lambda > 0)
            lambda = sqrt(lambda);
        else
            lambda = sqrt(-lambda);
    } else {
        zbar = 0.0L;
        lambda = 0.0L;
    }
#endif

    long nsi,nsi1;
    nsi = 0;
    nsi1 = 0;
    for (i=1;i<=nfg;i++){
        if (nZcont[i]>0)
        nsi = nsi + 1;
        if (nZcont[i]<=1)
        nsi1 = nsi1 + 1;
    }
#ifdef FLUID
    fprintf(out,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%ld\t%ld\t%ld\t%ld\t%Le\t%Le\t%Le\t%Le\n", (double long)it*dt, cineti, rotational, elastic, gravitational, fluid, ufluid, phib, q, n, v, zbar, lambda, ncontot, ncontglis, nsi, nsi1, viol, violmean, fxdown, -fydown);
#else
    fprintf(out,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%ld\t%ld\t%ld\t%ld\t%Le\t%Le\t%Le\n", (double long)it*dt, cineti, rotational, elastic, gravitational, ncontot, ncontglis, nsi, nsi1, viol, violmean, fxdown, -fydown);
#endif
    fprintf(out1,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %ld %ld %ld\n", (double long)it*dt, cineti, rotational, elastic, gravitational, -fydown, fxdown, fyup, -fxup, height, wallx, ncontglis, ncontot, nlt);
    fflush(out);
    fflush(out1);
    return;
}

/****************Save Grains***********************/
void savegrains(const char * name){
    FILE *out2;

    if ((out2 = fopen (name,"w")) == NULL)
    {
        fprintf(stderr,"Probleme d'ouverture de %s\n",name);
        exit (1);
    }
    fwrite(&ngtot,sizeof(ngtot),1,out2);
    fwrite(radius,sizeof(double long),ngtot+1,out2);
    fwrite(rx,sizeof(double long),ngtot+1,out2);
    fwrite(ry,sizeof(double long),ngtot+1,out2);
    fwrite(vx,sizeof(double long),ngtot+1,out2);
    fwrite(vy,sizeof(double long),ngtot+1,out2);
    fwrite(ome,sizeof(double long),ngtot+1,out2);
    fwrite(ax,sizeof(double long),ngtot+1,out2);
    fwrite(ay,sizeof(double long),ngtot+1,out2);
    fwrite(domedt,sizeof(double long),ngtot+1,out2);
    fwrite(&ncontfree,sizeof(long),1,out2);
    fwrite(&ncontot,sizeof(long),1,out2);
    fwrite(ior,sizeof(long),ncontot+1,out2);
    fwrite(iex,sizeof(long),ncontot+1,out2);
    fwrite(react,sizeof(double long),ncontot+1,out2);
    fflush(out2);
    fclose(out2);
}

#ifdef FLUID
/****************Save Fluid************************/
void savefluid(const char * name){
    FILE *out_fluid1;

    if ((out_fluid1 = fopen (name,"w")) == NULL) {
        fprintf(stderr,"Probleme d'ouverture de %s\n", name);
        exit (1);
    }
    fwrite(&H,sizeof(H),1,out_fluid1);
    fwrite(u,sizeof(double long),H,out_fluid1);
    fwrite(tau,sizeof(double long),H,out_fluid1);
    fwrite(dpressuredz,sizeof(double long),H,out_fluid1);
    fwrite(fxFd,sizeof(double long),H,out_fluid1);
    fwrite(fyFd,sizeof(double long),H,out_fluid1);
    fwrite(Phi,sizeof(double long),H,out_fluid1);
    fflush(out_fluid1);
    fclose(out_fluid1);
}

void savefluidprofile(const char * name){
    FILE *out_fluid1;
    long j = 0;

    if ((out_fluid1 = fopen (name,"w")) == NULL) {
        fprintf(stderr,"Probleme d'ouverture de %s\n", name);
        exit (1);
    }
    fprintf(out_fluid1, "#z\tu(z)\ttau(z)\tdpress(z)/dz\tfxFd(z)phi(z)/(1-phi(z))\tfyFd(z)phi(z)/(1-phi(z))\n");
    for (j = 0; j < H; j++)
        fprintf(out_fluid1, "%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", (j-LISTAT)*dz, u[j], tau[j], dpressuredz[j]/fluid_density, fxFd[j]*Phi[j]/(fluid_density*(1.0L-Phi[j])), fyFd[j]*Phi[j]/(fluid_density*(1.0L-Phi[j])));
    fflush(out_fluid1);
    fclose(out_fluid1);
}
#endif

/****************Save Grain Profile****************/
void savegrainsprofile(const char * name){
    FILE *out2;
    long i;

    if ((out2 = fopen (name,"w")) == NULL)
    {
        fprintf(stderr,"Probleme d'ouverture de %s\n",name);
        exit (1);
    }
    fprintf(out2,"#z\tphi\tvx\tvy\tfx\tfy\tfxarch\tfyarch\tfxdrag\tfydrag\tp\tcontact\tsliding\tinterp\n");
    for (i = 0; i < H; i++)
        fprintf(out2,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", (double long)(i-LISTAT)*dz, Phi[i], Vx[i], Vy[i], Fx[i], Fy[i], fxFa[i]*width*dz, fyFa[i]*width*dz, fxFd[i]*width*dz, fyFd[i]*width*dz, P[i], Zcont[i], Zglis[i], interp[i]);

    fflush(out2);
    fclose(out2);
}

/****************period save***********************/
void periodsave(){
    sprintf(periodsaveout, "grain%ld.tsv",nperiodsave);
    savegrainsprofile(periodsaveout);
#ifdef FLUID
    sprintf(periodsaveout, "fluid%ld.tsv",nperiodsave);
    savefluidprofile(periodsaveout);
    sprintf(periodsaveout, "fluid%ld.bin",nperiodsave);
    savefluid(periodsaveout);
#endif
    sprintf(periodsaveout, "save%ld.bin",nperiodsave);
    savegrains(periodsaveout);
    nperiodsave++;
    return;
}

/*********** Save Configuration ***********/
void saveconfiguration(){
#ifdef FLUID
    savefluid("fluidout.bin");
#endif
    savegrains("dataout.bin");
    return;
}

/*********** Close and save final file ***********/
void closeall(){
    fclose(out);
    fclose(out1);

    return;
}

/*********** Main ***********/
int main (){
    init();
    findneighbours2();
    predictor();
    detectcontacts();
    forcecalculation();
    cumulatestats();
    timemeasurement();
    periodsave();

    fprintf(out1,"#\n#\n#En entrée test, ncontfree=%ld  et    %ld \n",ncontfree,ncontot);

    it=0;
    endreached=1;
    while(endreached!=0){
        it++;
        if ((it%ifreqv)==0)
            findneighbours2();
        predictor();
        detectcontacts();
        forcecalculation();
#ifdef FLUID
        if ((menu >= 33) && (menu < 45)){
            calculateprofiles();
        }
#endif
        if ((it%icustats)==0)
            cumulatestats();
        if ((it%itimemeasurements)==0)
            timemeasurement();
        if ((it%iperiodsave)==0)
            periodsave();
        if ((it%isave)==0)
            saveconfiguration(); //Except from angular position, all state variables are saved, including tangential forces
        if (it>nbiteration){         //nbr pas de temps (en unité de dt)
            endreached=0;
        }
        corrector();
#ifdef FLUID
        if ((menu >= 30) && (menu < 45)){
            if ((menu == 30) || (menu == 31));
            // Just let grains fall under gravity
            else if (menu == 32)
            // Impose fluid velocity proflie that shear grains
                sheargrainswithfluid();
            else {
                calculateprofiles();
                if ((menu == 33)||(menu == 35)||(menu == 37))
                // Solve viscous fluid velocity using implicit scheme
                    solvefluidvelocity();
                else if ((menu == 34)||(menu == 36)||(menu == 38))
                // Solve fluid velocity using explicit Orencio's model for mixing length
                    solvefluidvelocity3();
            }
        }
#endif
    }
    if (extract==1){
        fprintf(out1,"#extract wall \n");
        extractwall();
        findneighbours2();
        predictor();
        detectcontacts();
        forcecalculation();
        fprintf(out1,"#free grain number=%ld \n",nfg);
        fprintf(out1,"#number of grains in each wall=%ld \n",nwall);
    }
    fprintf(out1,"#'end reached !'\n");
    findneighbours2();
    predictor();
    detectcontacts();
    forcecalculation();
    cumulatestats();
    timemeasurement();
    periodsave();
    saveconfiguration();
    fprintf(out1,"#\n#\n#En sortie test, ncontfree=%ld  et    %ld \n",ncontfree,ncontot);

    closeall();
    return 0;
}

/*********** Complementary functions ***********/

/*********** random number generators CHECK ***********/
double long ran0()
{
    const long iA=843314861;
    const long iB=453816693;
    const long iM=1073741824;
    double long aux,x;

    aux = 0.5L/(double long)(iM);
    isem = isem*iA + iB;
    if (isem < 0)
    {
        isem = (isem + iM )+iM;
    }
    x = isem*aux;
    return x;
}

/*********** random number generator between 0 and 1 ***********/
double long ran2(long *idum)
{
    const long IM1=2147483563, IM2=2147483399;
    const long IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
    const long IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1;
    const long NDIV=1+IMM1/NTAB;
    const double long EPS=3.0e-16L;
    const double long RNMX=1.0L-EPS;
    const double long AM=1.0L/(double long) IM1;
    static long idum2=123456789, iy=0;
    static long iv[32];
    long j,k;
    double long temp;

    if (*idum <= 0) {
        *idum = (*idum==0 ? 1 : -*idum);
        idum2 = *idum;
        for (j=NTAB+7; j>=0;j--) {
            k=*idum/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum+=IM1;
            if (j < NTAB) iv[j] = *idum;
        };
        iy=iv[0];
    };
    k=*idum/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum < 0) *idum+=IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j]=*idum;
    if(iy < 1) iy += IMM1;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
};

double gaussian(long *idum)
{
    static double long t = 0.0L;
    double long x, w1, w2, r;

    if( t == 0 )
    {
        do
        {
            w1 = 2.0L * ran2(idum) - 1.0L;
            w2 = 2.0L * ran2(idum) - 1.0L;
            r = w1 * w1 + w2 * w2;
        }
        while( r >= 1.0 );
        r = sqrt( -2.0*log(r) / r );
        t = w2 * r;
        return(w1 * r);
    }
    else
    {
        x = t;
        t = 0.0;
        return(x);
    }
}
