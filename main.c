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

#define DEBUGPRINT(...)       printf(__VA_ARGS__);   fflush(stdout)
#define FLUID
/// Constants defined to calculate the profiles
#define NSTAT 500 // Number of divisions in the profile in function of z
#define dz 0.1L // For now it is a constant, but can be funtion of z it self
const long LISTAT = (int)(3.0L/dz); // Inferior limit to negative values of z
#define MAXGRAINS 10000

#ifdef FLUID
const double long PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609;
#endif

/*********** declaration functions ***********/
double long ran0();
double long ran2(long *idum);
double gaussian(long *idum);
double gaussian2(long *idum);
void sauveconf();
double long distperio(double long x);
void findneighbours();
void findneighbours2();
void closeall();
void init();
void extractwall();
void prepareconfig();
void predictor();
void detectcontacts();
void forcecalculation();
void corrector();
void checkendreachedibrium();
void cumulatestats();
void timemeasurement();
void periodsave();
void saveconfiguration();
#ifdef FLUID
void solvefluidvelocity(); // Function that solves the fluid velocity, a partial differential equation (a non-linear equation in diffusion family).
void savestat(); // Function to save the statistics each time of saving the data.
void particlesprofile(); // Function that calculate the profile of the particles to interact with the fluid.
void fluidforce(long); // Function that calculate the force that the fluid exert over a grain. The parameter passed is the identification number of the grain.
void savefluidconfiguration();
#endif

/*********** global variables ***********/

char nomin[30]="configin.txt";
char datanameinold[30]="datain.txt";
char datanamein[30]="datain.bin";
char datanameout[30]="dataout.bin";
char nomout[30]="timeout.igr";
char nomout2[30]="statout.igr";
char output[30]="output.txt";
char periodsaveout[30];
long it,isem = 987654321;

long extract;
long ncontfree,ncontot,nl,nlt,ifreqv,icont0,ncontglis;
long isave; //period over which configuration (state variables) AND cumulated statistics are saved
long itimemeasurements; //period over which time dependent quantities are measured and saved
long iperiodsave; //period over which are saved fields needed for Ep, \Gamma calculation
long icustats; //period over which quantities are evaluated to get cumulated statistics
double long ivmass[MAXGRAINS],ivmomine[MAXGRAINS],sumfn[MAXGRAINS*10],reacn[MAXGRAINS*10], reactsave[MAXGRAINS*10],react[MAXGRAINS*10], ax[MAXGRAINS],ay[MAXGRAINS],ome[MAXGRAINS],domedt[MAXGRAINS],phi[MAXGRAINS];
double long statZglis[NSTAT],statZcont[NSTAT],statPhi[NSTAT],statVx[NSTAT],statVy[NSTAT],statFx[NSTAT],statFy[NSTAT],statP[NSTAT];
double long kn; //normal spring constant
double long kt; //tangential spring constant
double long dt; //time step
double long gn; //normal damping factor
double long frott; //friction coefficient
double long halfdt2,halfdt,tolff,viol;
long nwall; //number of grains in each wall;
double long aa; //fit parameter fo P/phi for menu==4
double long bb; //fit parameter fo P/phi for menu==4
long nfg; //number of free grains
long ngtot; //total number of grains
long menu; //select the configuration
long NbinsZ;
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
double long vx[MAXGRAINS],vy[MAXGRAINS],radius[MAXGRAINS],rx[MAXGRAINS],ry[MAXGRAINS],rpx[MAXGRAINS],rpy[MAXGRAINS], vpx[MAXGRAINS],vpy[MAXGRAINS],xwall[500],ywall[500];
long ior[MAXGRAINS*10],iex[MAXGRAINS*10];
long listi[MAXGRAINS*10],listj[MAXGRAINS*10];
double long apx[MAXGRAINS],apy[MAXGRAINS],phip[MAXGRAINS],omep[MAXGRAINS],domedtp[MAXGRAINS],Zcont[MAXGRAINS*10],Zglis[MAXGRAINS*10],eij[MAXGRAINS*10],xnij[MAXGRAINS*10],ynij[MAXGRAINS*10];
double long fpx[MAXGRAINS],fpy[MAXGRAINS],gam[MAXGRAINS];
long io[MAXGRAINS*10];
long endreached;
long nq,igr;
double long deltar; // width of the uniform distribution used for the radius
long islu; //seed for the random number generator
long ncumul;
long nperiodsave=0; //index of the current file saving configurations
double long fydown,fxdown,fyup,fxup; //Normal and tangiential stresses on the walls
double long springxx,springyy,springxy;
// The 2 variables, density and viscosity, are used to describle the fliud. All considerations were based on Phys. Fluids 24, 103306
#ifdef FLUID
double long density_ratio, viscosity;
// The next 3 variables are control parameters to the system:
double long dragc; //Drag coefficient, function of Granular Reynolds' number
double long reynolds; // Reynolds' number in grain dimensions, funciton of density ratio (fluid and grain), gravity, grain diameter and viscosity
double long shields; // Shields' number, funciton of density (fluid and grain), the velocity field of the fluid, gravity and grain diameter
// the next 2 variables are the force interaction between grain and fluid: drag and Archimedes
double long fdrag; // Drag force, function of velocity (grain and fluid), drag coefficient, fliud's density, and grain area
double long farch; // Archimedes force, function of the fluid stress on the grain. See Phys. Fluids 24, 103306 - eq. 5
// The next 2 constants are the only constants that defines the fluid it self
const double long dragd = 0.5L; // Drag coefficient of the grain in the turbulent regime (Granular Reynolds -> inf)
const double long reynoldsc = 24.0L; // Reynolds' number where drag coefficient is almost constant
// Some other variables to calculate the force and other parameters of the interaction between fluid and grain
double long ustar; // Characteristic velocity of the fluid - something like a shear velocity
const double long z0 = 0.1L; // Characteristic length of turbulence effect
const double long kapa = 0.4L; // Characteristic dimensionless constant of the fluid velocity field
double long tau_inf; // Shear stress without granular interaction, on the direction of the flowing
double long zb; // The height that bed forms
double long reynoldsg; // The Reynolds number on the grain interaction
double long grain_density; // The density of the grain (Associated with the mass of the grain)
double long u[NSTAT]; // The profile of the fluid velocity.
double long fxfluid[MAXGRAINS], fyfluid[MAXGRAINS]; // The force the grains exert on the fluid, in x and z directions
double long fxF[NSTAT], fyF[NSTAT]; // The force the fluid exert on the grains per volume, in x and z directions
double long tau[NSTAT+1]; // Derivative of the shear stress along z direction, phased shifted dz/2
double long dpressuredz[NSTAT+1]; // Derivative of the pressure along z direction, phased shifted dz/2
double long packfraction[NSTAT]; // Instantaneous packing fraction of the mesh
double long instZglis[NSTAT],instZcont[NSTAT],instPhi[NSTAT],instVx[NSTAT],instVy[NSTAT],instFx[NSTAT],instFy[NSTAT],instP[NSTAT],instFn[NSTAT],instFt[NSTAT],instInterp[NSTAT]; // Instantaneous variables of measurement
long d = 1;
#endif
FILE *out;
FILE *out1;
#ifdef FLUID
FILE *out_fluid; // File with fluid parameters, velocities, forces, packing fraction...
#endif
FILE *outest;
long nbiteration;
double long max_x, max_y, min_x, min_y; // Greater and lower positions of particles in x and y directions

void init()
{
/*
In this simulation program, the menu is used to choice which kind of simulation will be done:
From 0 to 14, one bottom inclined wall with gravity is used to simulate the sistem.
From 15 to 29, two walls, one at bottom and other at top, without gravity, are used to simulate shear systems.
Form 30 to 44, one bottom horizontal wall with gravity is used to simulate fluid systems.
Funcitons 0, 15, and 30 create systems. (In particular, 30 has the same function of 0)
*/
long i,ntemp;
char comment[255];
FILE *in;
FILE *in2;

if ((out1 = fopen (output,"w")) == NULL)
{
    fprintf(stderr,"Probleme d'ouverture de %s\n",output);
    exit (1);
}
srand(time(NULL));
islu=rand();
in = fopen(nomin, "rt") ;
if (in == NULL)
{
fprintf (out1,"no input file...........\n");
exit (1);
}
else fprintf(out1,"input file detected \n");
fscanf (in,"%ld\t%s\n", &nfg,comment);
fprintf(out1,"free grain number=%ld \n",nfg);
fscanf (in,"%ld\t%s\n", &nwall,comment);
fprintf(out1,"number of grains in each wall=%ld \n",nwall);
fscanf (in,"%Le\t%s\n", &width,comment);
fprintf(out1,"cell width (unit d) width=%Le \n",width);
halfwidth=0.5L*width;
fscanf(in,"%Le\t%s\n", &deltar,comment);
fprintf(out1,"Radius picked at random with a uniform law of width %Le\n",deltar);
fscanf (in,"%Le\t%s\n", &kn,comment);
fprintf(out1,"normal spring constant(?) kn=%Le \n",kn);
fscanf (in,"%Le\t%s\n", &kt,comment);
fprintf(out1,"kt=%Le \n",kt);
fscanf (in,"%Le\t%s\n", &gn,comment);
fprintf(out1,"damping coefficient=%Le normally smaller than %Le \n",gn,2.0L*sqrt(kn));
fscanf (in,"%Le\t%s\n", &frott,comment);
fprintf(out1,"friction coefficient frott=%Le \n",frott);
fscanf (in,"%ld\t%s\n", &itimemeasurements,comment);
fprintf(out1,"time dependent measurement period=%ld \n",itimemeasurements);
fscanf (in,"%ld\t%s\n", &iperiodsave,comment);
fprintf(out1,"fraction of dot_gamma for periodic meseasurements=%ld \n",iperiodsave);
fscanf (in,"%ld\t%s\n", &ifreqv,comment);
fprintf(out1,"neighbour list refreshing period ifreqv=%ld \n",ifreqv);
fscanf (in,"%ld\t%s\n", &icustats,comment);
fprintf(out1,"statistical quantities evaluation period=%ld \n",icustats);
fscanf (in,"%ld\t%s\n", &isave,comment);
fprintf(out1,"configuration and statistical quantities saving period=%ld \n",isave);
fscanf (in,"%ld\t%s\n", &menu,comment);
fscanf (in,"%ld\t%s\n", &extract,comment);
fprintf(out1,"Extraction=%ld \n",extract);
fscanf (in,"%ld\t%s\n", &nbiteration,comment);
fprintf(out1,"nbiteration=%ld \n",nbiteration);
#ifdef FLUID
    if ((menu >= 30)&&(menu < 45)){
        gravity = 1.0L;
        density_ratio = 2.0L;
        reynolds = 10.0L;
        shields = 0.0L;
        fscanf(in,"%Le\t%s\n", &density_ratio, comment);
        fprintf(out1,"Density ratio of the grain over fluid=%Le\n", density_ratio);
        fscanf(in,"%Le\t%s\n", &reynolds, comment);
        fprintf(out1,"Fluid's Reynolds number=%Le\n", reynolds);
        fscanf(in,"%Le\t%s\n", &shields, comment);
        fprintf(out1, "Shields number=%Le\n", shields);
        if (menu == 30){
            shields = 0.0L;
        }
    }
#endif

if ((menu>14) && (menu < 30))
    ngtot=nfg+2*nwall;
else
    ngtot=nfg+nwall;
fprintf(out1,"Total number of grains=%ld \n",ngtot);
if (menu<15)
{
    fscanf (in,"%Le\t%s\n", &theta,comment);
    fprintf(out1,"theta theta=%Le \n",theta);
    costheta=cos(theta);
    sintheta=sin(theta);
    fprintf(out1,"Inclined plane configuration\n");
	if ((menu==2) || (menu==3) || (menu==4))
	{
		fprintf(out1,"Includes pseudo-wall friction\n");
		fscanf (in,"%Le\t%s\n", &muwall,comment);
		fprintf(out1,"muwall=%Le \n",muwall);
		muwall/=(2.0L*M_PI);
	}
    if (menu==4)
    {
        fprintf(out1,"Includes fit expression of P/phi in pseudo-wall friction\n");
        fscanf (in,"%Le\t%s\n", &aa,comment);
        fprintf(out1,"aa=%Le \n",aa);
        fscanf (in,"%Le\t%s\n", &bb,comment);
        fprintf(out1,"bb=%Le \n",bb);
    }
}
else if ((menu>14)&&(menu<28))
{
    fprintf(out1,"Pressure-velocity controlled configuration\n");
    fscanf (in,"%Le\t%s\n", &wallv,comment);
    fprintf(out1,"wallv U=%Le \n",wallv);
	if (menu==25)
	{
		fprintf(out1,"Gravity/anti-gravity added in five layers on both sides\n");
		fscanf (in,"%Le\t%s\n", &gravity,comment);
		fprintf(out1,"gravity=%Le \n",gravity);
	}
	if (menu==26)
	{
		fprintf(out1,"Gravity added\n");
		fscanf (in,"%Le\t%s\n", &gravity,comment);
		fprintf(out1,"gravity=%Le \n",gravity);
	}
	if (menu==27)
	{
		wallvzero=wallv;
		fscanf (in,"%Le\t%s\n", &timewall,comment);
		fprintf(out1,"time-scale decay=%Le \n",timewall);
	}
}
else if (menu==28)
{
    fprintf(out1,"Volume-velocity controlled configuration\n");
    fscanf (in,"%Le\t%s\n", &wallv,comment);
    fprintf(out1,"wallv U=%Le \n",wallv);
}
else if (menu==29)
{
    fprintf(out1,"Pressure-Stress controlled configuration\n");
    fscanf (in,"%Le\t%s\n", &shear,comment);
    fprintf(out1,"shear Tau=%Le \n",shear);
}

if ((menu==0)||(menu==15)||(menu==30))
{
    if (menu==0)
        fprintf(out1,"menu 0: prepare from scratch with bottom wall) \n");
    else if (menu==15)
        fprintf(out1,"menu 15: prepare from scratch with two walls) \n");
#ifdef FLUID
    else if (menu==30)
        fprintf(out1,"menu 30: prepare from scratch with bottom wall) \n");
#endif
    fprintf(out1,"-------------------------------------------------\n");

	//CHECK ADD HERE ROUGHNESS PARAMETERS
    prepareconfig();
}
else
{
    fprintf(out1,"Start from a saved configuration\n");
	in2 = fopen(datanamein, "rt") ;
	if (in2 == NULL){
		fprintf (out1,"no binary file -> old text file...........\n");

    in2=fopen(datanameinold,"r");
    fscanf(in2,"%ld\n",&ntemp);
    if (ntemp!=ngtot)
        fprintf(out1,"\n\n\n\n\n\nProblème grave dans le nombre de points\n\n\n\n\n\n\n\n\n\n");
    for (i = 1; i<=ngtot; i++)
        fscanf(in2,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", &radius[i],&rx[i], &ry[i],&vx[i], &vy[i],&ome[i],&ax[i], &ay[i],&domedt[i]);
    fscanf(in2,"%ld\n",&ncontfree);
    fscanf(in2,"%ld\n",&ncontot);
    for (i = 1;i<= ncontot;i++)
        fscanf(in2,"%ld\t%ld\t%Le\n", &ior[i],&iex[i],&react[i]);
    fprintf(out1,"\n\nEn entrée, ngtot=%ld\n",ntemp);
    fprintf(out1," ncontfree=%ld  et    %ld \n",ncontfree,ncontot);
	}
	else
	{
		fprintf (out1,"Binary file detected...........\n");
		fread(&ntemp,sizeof(ntemp),1,in2);
		if (ntemp!=ngtot)
			fprintf(out1,"\n\n\n\n\n\nProblème grave dans le nombre de points: ntemp=%li ngtot=%li\n\n\n\n\n\n\n\n\n\n",ntemp, ngtot);
		fprintf(out1,"\n\nEn entrée, ngtot=%ld\n",ntemp);
		fread(radius,sizeof(double long),ngtot+1,in2);
		fread(rx,sizeof(double long),ngtot+1,in2);
		fread(ry,sizeof(double long),ngtot+1,in2);
		fread(vx,sizeof(double long),ngtot+1,in2);
		fread(vy,sizeof(double long),ngtot+1,in2);
		fread(ome,sizeof(double long),ngtot+1,in2);
		fread(ax,sizeof(double long),ngtot+1,in2);
		fread(ay,sizeof(double long),ngtot+1,in2);
		fread(domedt,sizeof(double long),ngtot+1,in2);
		fread(&ncontfree,sizeof(ncontfree),1,in2);
		fread(&ncontot,sizeof(ncontot),1,in2);
		fprintf(out1," ncontfree=%ld  et    %ld \n",ncontfree,ncontot);
		fread(ior,sizeof(double long),ncontot+1,in2);
		fread(iex,sizeof(double long),ncontot+1,in2);
		fread(react,sizeof(double long),ncontot+1,in2);
	}
	fclose(in2);
	fclose(in);
}

totalweight=0.0L;

/********* Mass normalisation*/
#ifndef FLUID
for (i=1;i<=nfg;i++)
{
    totalweight+=4.0L*radius[i]*radius[i];
    ivmass[i] = 0.25L / (radius[i]*radius[i]); //inverse of the grain mass
    ivmomine[i] = 2.0L * ivmass[i] /(radius[i]*radius[i]); //inverse of the grain moment of inertia
}
#endif
#ifdef FLUID
for (i=1;i<=nfg;i++)
{
    ivmass[i] = 1.0L/(8.0L*radius[i]*radius[i]*radius[i]); //inverse of the grain mass (normalizated)
    totalweight+=1.0L/ivmass[i]; // Density of grain != 1 (6/PI)
    ivmomine[i] = 2.5L * ivmass[i] /(radius[i]*radius[i]); //inverse of the grain moment of inertia of a sphere
}
grain_density = 6.0L/PI; // Considering a sphere as volume, with diameter=1 and mass=1
#endif
fprintf(out1,"Mass =%Le\n",totalweight);
// Initialization of Granular variables of momentum
for (i=nfg+1;i<=ngtot;i++)
{
    rpx[i]=rx[i];
    rpy[i]=ry[i];
	vpx[i]=0.0L;
	vpy[i]=0.0L;
	omep[i]=0.0L;
	vx[i]=0.0L;
	vy[i]=0.0L;
	ome[i]=0.0L;
}
for (i=1; i<=nfg; i++){
    phi[i] = 0.0L;
}
// Initialization of Fluid variables and constants
#ifdef FLUID
    if ((menu >= 30) && (menu < 45)){
        viscosity = sqrt(density_ratio -1.0L)/reynolds;
        ustar = sqrt(shields*(density_ratio -1.0L)*gravity);
        bottom=0.0L;
        for (i=nfg+1;i<=nfg+nwall;i++)
            bottom+=ry[i]+radius[i];
        bottom/=(double long)nwall;
        zb = bottom;
        fprintf(out1, "Grain diameter: %Le\n", 1.0L);
        fprintf(out1, "Gravity: %Le\n", gravity);
        fprintf(out1, "Density ratio (Grain/Fluid): %Le\n", density_ratio);
        fprintf(out1, "Reynolds number of the fluid: %Le\n", reynolds);
        fprintf(out1, "Shields number: %Le\n", shields);
        fprintf(out1, "Viscosity: %Le\n", viscosity);
        fprintf(out1, "Shear velocity: %Le\n", ustar);
        fprintf(out1, "Bed position: %Le Turbulence effect: %Le Fluid velocity constant %Le\n", zb, z0, kapa);
        fflush(out1);
    }
#endif
if (menu<15)
{
    fprintf(out1,"Pressure =%Le\n",totalweight*costheta/(double long)width);
    fprintf(out1,"Shear stress =%Le\n",totalweight*sintheta/(double long)width);
}
else if (menu < 30)
{
	bottom=0.0L;
	for (i=nfg+1;i<=nfg+nwall;i++)
		bottom+=ry[i];
	bottom/=(double long)nwall;
    height=0.0L;
    for (i=nfg+nwall+1;i<=ngtot;i++)
        height+=ry[i];
    height/=(double long)nwall;
	fprintf(out1,"Bottom=%Le  Height=%Le\n",bottom,height);
    for (i=1;i<=nwall;i++)
    {
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
        fprintf(out1,"Pressure cause by grains on the wall =%Le\n",totalweight/(double long)width);
    }
#endif
if ((out = fopen (nomout,"w")) == NULL)
{
    fprintf(stderr,"Probleme d'ouverture de %s\n",nomout);
    exit (1);
}
fprintf(out,"IGOR\n\n");
fprintf(out,"WAVES/O/R/S it\tncontot\tnsi\tnsi1\tviol\tEk\tEp\tfydown\tfxdown\tfyup\tfxup\theight\twallx\tncontglis\n");
fprintf(out,"BEGIN\n");

#ifdef debug
	{
	/** TEST TEST**/
	if ((outest = fopen ("test.igr","w")) == NULL)
	{
		fprintf(stderr,"Probleme d'ouverture de %s\n",nomout);
		exit (1);
	}
	fprintf(outest,"IGOR\n\n");
	fprintf(outest,"WAVES/O/R/S ta\tja\til\tft\tfpre\n");
	fprintf(outest,"BEGIN\n");

	/** END TEST TEST**/
	}
#endif

/** Initialisation of statistical profiles*/
ncumul=0;
for (i = 0 ;i < NSTAT;i++)
{
    statPhi[i] = 0.0L;
	statP[i] = 0.0L;
    statVx[i] = 0.0L;
    statVy[i] = 0.0L;
    statFx[i] = 0.0L;
    statFy[i] = 0.0L;
    statZcont[i] = 0.0L;
    statZglis[i] = 0.0L;
#ifdef FLUID
    fxF[i] = 0.0L;
    fyF[i] = 0.0L;
    packfraction[i] = 0.0L;
    if (menu == 33){
        FILE *config_fluid;
        char tmp[200];
        long double pos;
        long j = 0;
        if ((config_fluid = fopen("fluidin.txt","r")) == NULL){
            fprintf(stderr, "Probleme d'ouverture de %s\n","fluidin.txt");
            exit(1);
        }
        fgets(tmp, 200 , config_fluid);
        while(!feof(config_fluid)){
            fscanf(config_fluid, "%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", &pos, &u[j], &tau[j+1], &dpressuredz[j+1], &fxF[j], &fyF[j]);
            j+= (j < NSTAT) ? 1 : 0;
        }
        fclose(config_fluid);
        if (j != NSTAT){
            fprintf(stderr, "ERROR, Fluid config does not match with the fluid's mesh. Number of lines read %li.\n", j);
        }
    }
    else if (menu == 32){
    // First try to set the fluid profile near to the stationary regime.
/*
        if ((i-LISTAT)*dz < 10*viscosity)
            u[i] = (i-LISTAT)*dz-0.8L*nfg/width-0.0L > 0.0L ? ustar*ustar*(i-LISTAT)*dz/viscosity : 0.0L;
        else if ((i-LISTAT)*dz < 50*viscosity)
            u[i] = (i-LISTAT)*dz-0.8L*nfg/width-0.0L > 0.0L ? ustar*ustar*(i-LISTAT)*dz/viscosity : 0.0L;
        else
            u[i] = (i-LISTAT)*dz-0.8L*nfg/width-0.0L > 0.0L ? ustar*log(1.0L+((i-LISTAT)*dz-0.8L*nfg/width-0.0L)/z0)/kapa : 0.0L;
*/
    // Second try to set the fluid profile near to the stationary regime.
/*
        // Should be greater than the granular bed, with 3 different functions
        long double h_position = (i-LISTAT)*dz-0.8L*nfg/width-0.0L;
        if (h_position > 0.0L){
            long double A = (2.0L -log(1.0L+50.0L*viscosity/(ustar*z0)))/(log(1.0L+5.0L*viscosity/(ustar*z0))-log(1.0L+50.0L*viscosity/(ustar*z0)));
            long double B = (2.0L*log(1.0L+50.0L*viscosity/(ustar*z0))-log(1.0L+5.0L*viscosity/(ustar*z0))*log(1.0L+50.0L*viscosity/(ustar*z0)))/(log(1.0L+5.0L*viscosity/(ustar*z0))-log(1.0L+50.0L*viscosity/(ustar*z0)));
            u[i] = (1.0L-tanh(h_position-5.0L*viscosity/ustar))/2.0L * ustar*ustar * h_position/viscosity; // Viscous Regime
            u[i] += (1.0L-tanh(h_position-50.0L*viscosity/ustar))/2.0L * (1.0L+tanh(h_position-5.0L*viscosity/ustar))/2.0L * ustar * (A*log(1.0L+h_position/z0) -B)/kapa; // Viscous-Turbulent Transition
            u[i] += (1.0L+tanh(h_position-50.0L*viscosity/ustar))/2.0L * ustar * (log(1.0L+h_position/z0))/kapa; // Turbulent Regime
        } else {
            u[i] = 0.0L;
        }
        tau[i] = 0.0L;
        dpressuredz[i] = -grain_density*gravity/density_ratio;
*/
    // Third try to set the fluid profile near to the stationary regime.
        int j = 0;
        long double maxwall = 0.0L;
        for (j = nfg+1; j <= ngtot; j++){
            maxwall = maxwall < rx[j] + radius[j] ? rx[j] + radius[j] : maxwall;
        }
        long double h_position = (i-LISTAT)*dz-0.8L*nfg/width-maxwall;
        if (h_position > 0.0L){
            u[i] = ustar*log(1.0L+h_position/z0)/kapa;
        } else {
            u[i] = 0.0L;
        }
        tau[i] = ustar*ustar*grain_density/density_ratio;
        dpressuredz[i] = -grain_density*gravity/density_ratio;
    }
    else {
        u[i] = 0.0L;
        tau[i] = 0.0L;
        dpressuredz[i] = -grain_density*gravity/density_ratio;
    }
#endif
}

#ifdef FLUID
    tau[NSTAT] = ustar*ustar;
    dpressuredz[NSTAT] = -grain_density*gravity/density_ratio;
    if ((menu >= 30) && (menu < 45)){
        if ((out_fluid = fopen ("fluid.txt","w")) == NULL)
        {
            fprintf(stderr,"Probleme d'ouverture de %s\n","fluid.txt");
            exit (1);
        }
        fprintf(out_fluid, "Grain diameter: %Le\n", 1.0L);
        fprintf(out_fluid, "Gravity: %Le\n", gravity);
        fprintf(out_fluid, "Density ratio (Grain/Fluid): %Le\n", density_ratio);
        fprintf(out_fluid, "Reynolds number of the fluid: %Le\n", reynolds);
        fprintf(out_fluid, "Shields number: %Le\n", shields);
        fprintf(out_fluid, "Viscosity: %Le\n", viscosity);
        fprintf(out_fluid, "Shear velocity: %Le\n", ustar);
        fprintf(out_fluid, "Bed position: %Le Turbulence effect: %Le Fluid velocity constant %Le\n", zb, z0, kapa);
        fflush(out_fluid);
    }
#endif

/*********** Determination time step  ***********/
dt = 0.02L/sqrt(kn); /* ! xx/50 -> arbitraire!*/
#ifdef FLUID
    if ((menu > 30) && (menu < 45)){
        // Checking the condition to do the simulation over the turbulent regime:
        if ((dt > dz*dz/(NSTAT)) || (dt > dz/(ustar*ustar)) || (dt > dz*dz/viscosity)){
        // Means that the dt of the fluid is smaller than the dt of the grains, to solve the turbulent fluid with stability
            dt = dt > dz*dz/(NSTAT) ? dz*dz/(NSTAT) : dt;
            dt = dt > dz*dz/viscosity ? dz*dz/viscosity : dt;
            dt = dt > dz/(ustar*ustar) ? dz/(ustar*ustar) : dt;
        }
    }
#endif
fprintf(out1,"dt = %Le \n",dt);
if (menu==27)
	timewall=dt/timewall;

/********** Coefficients of the predictor-corrector algorithm  ***********/
halfdt2 = 0.5L*dt*dt;//CHECK that this is indeed predictor-corrector
halfdt = 0.5L*dt;
}

void extractwall()
{
long i,j;
double long rxlo,rylo,radiuslo,minlo;
long malo;

if ((menu>14) && (menu<30))
{
	height=0.0L;
	for (i=nfg+nwall+1;i<=ngtot;i++)
		height+=ry[i];
	height/=(double long)nwall;
	height-=4.0L;
}

i=0;
nwall=0;
nfg=ngtot;
while (i<nfg)
{
	i++;
	if (ry[i]<4.0L)
	{
		rxlo=rx[i];
		rylo=ry[i];
		radiuslo=radius[i];
		for (j=i;j<ngtot;j++)
		{
			rx[j]=rx[j+1];
			ry[j]=ry[j+1];
			radius[j]=radius[j+1];
		}
		if (rylo<2.0L)
		{
			ngtot--;
		}
		else
		{
			nwall++;
			rx[ngtot]=rxlo;
			ry[ngtot]=rylo;
			radius[ngtot]=radiuslo;
		}
		i--;
		nfg=ngtot-nwall;
	}
}

if ((menu>14) && (menu<30))
{
	fprintf(out1,"nfg tempo=%ld\n",nfg);
	for (i=1;i<=nwall;i++)
	{
		minlo=1.0E40L;
		for (j=1;j<=nfg-i+1;j++)
			if ((ry[j]>height)&&(ry[j]<minlo))
			{
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
	while (i<nfg)
	{
		i++;
		if (ry[i]>height)
		{
			for (j=i;j<ngtot;j++)
			{
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
for (i=1;i<=ngtot;i++)
	{
		ry[i]-=4.0L;
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
void prepareconfig()
{
long ncpt;
double long dx,xlo;
long i,j;

/*Pick up the radius at random according to a uniform law*/
for (i=1;i<=ngtot;i++)
    radius[i] = 0.5L + deltar*(ran2(&islu)-0.5L);

ncpt = 0;
i=-(halfwidth/2)-1;
j=1;
for (ncpt=1;ncpt<=nfg;ncpt++)
{
    i+=1;
    if (i>(halfwidth/2-1))
    {
        i=-(halfwidth/2);
        j++;
    }
    rx[ncpt] =2.0L*(double long)i+0.3L*ran2(&islu)+(double long)(j%2);
    ry[ncpt] =2.0L*(double long)j+0.3L*ran2(&islu);
}

//The boundary is consituted by disks
dx = width / (double long)(nwall);
xlo=-halfwidth;
for (i=1;i<=nwall;i++)
{
    rx[nfg+i] = xlo+0.2*ran2(&islu);
    ry[nfg+i] = 0.2L*ran2(&islu);
    xlo = xlo + dx;
}

if (menu==15)
{
    height=2.0L*(double long)(j+1.0L);
    xlo=-halfwidth;
    for (i=1;i<=nwall;i++)
    {
        rx[nfg+nwall+i] = xlo+0.2*ran2(&islu);
        ry[nfg+nwall+i] =height + 0.2L*ran2(&islu);
        xlo = xlo + dx;
    }
}

for (i=1;i<=ngtot;i++)
{
    vx[i] = 0.0L;
    vy[i] = 0.0L;
    ax[i] = 0.0L;
    ay[i] = 0.0L;
    phi[i] = 0.0L;
    ome[i] = 0.0L;
    domedt[i] = 0.0L;
}
}

/*********** Find neighbours in a crude way (in N^2)  ***********/
/*********** Require to know ior and iex ***********/
void findneighbours()
{
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
for (i=1;i<=nfg;i++)
{
    xi = rx[i];
    yi = ry[i];
    for (j=i+1;j<=nfg;j++)
    {
        xij = distperio(rx[j] - xi);
        yij = ry[j] - yi;
        dij = (xij*xij + yij*yij);
        if (dij < x12)
        {
            ineigh++;
            listi[ineigh] = i;
            listj[ineigh] = j;
            if ((ior[icon]==i)&&(iex[icon]==j))//We know in advance that the next one previously in contact must remain in the neighbourhood
            {
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
for (i=1;i<=nfg;i++)
{
    xi = rx[i];
    yi = ry[i];
    for (j =nfg+1;j<=ngtot;j++)
    {
        xij = distperio(rx[j] - xi);
        yij = ry[j] - yi;
        dij = (xij*xij + yij*yij);
        if (dij < x12)
        {
            ineigh++;
            listi[ineigh] = i;
            listj[ineigh] = j;
            if ((ior[icon]==i)&&(iex[icon]==j))
            {
                io[ineigh] = 1;
                icon ++;
            }
            else
                io[ineigh] = 0;
        }
    }
}
nlt = ineigh;// index reached at the end of free-wall neighbourhood

#ifdef debug
	{
		if (icon!=(ncontot+1))
			fprintf(out1,"\n\nBUG BUG BUG %ld %ld\n\n",icon,ncontot);
	}
#endif
return;
}

/*********** Find neighbours acording to cells (in NlogN - details at Computer Simulations of Liquids [Allen & Tildesley, 1991], Computational Granular Dynamics [Pöschel, 2004], Simulação de Alta Performance de Materiais Granulares [Martins, G. H. B. - 2015])  ***********/
/*********** Require to know ior and iex ***********/
void findneighbours2(){
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
    for(i = 2; i <= nfg; i++){
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
//DEBUGPRINT("%Le %i\n",l,(long)(nfg/l));
    long Ll = (long)(L/l) +1, Mm = (long)(M/m) +1; // Number of cells in x = Ll and y = Mm
    if ((Ll <= 5) || (Mm <= 5)){ // Too small system, solving by Findneighbours1
        findneighbours();
        return;
    }
    long cells[Ll][Mm][size], ccell[Ll][Mm]; // If a huge dispersion happens, the last fix number should be changed
    for (i = 0; i < Ll; i++){
        for (j = 0; j < Mm; j++){
            ccell[i][j] = 0; // Count the number of grains in each cell
            for (k = 0; k < size; k++){
                cells[i][j][k] = 0;
            }
        }
    }

    for (i = 1; i <= nfg; i++){
        x = (long) ((rx[i] -minx) / l);
        y = (long) ((ry[i] -miny) / m);
        // Just to be sure that it won't be outside of the boundaries, but by this definition, it shouldn't happen (USIING PERIODIC BOUNDARY CONDITIONS, IT CHANGE THE THINGS)
        x = (x < 1) ? 1 : x;
		x = (x >= Ll-1) ? Ll -2 : x;
		y = (y < 1) ? 1 : y;
		y = (y >= Mm-1) ? Mm -2 : y;
        cells[x][y][ccell[x][y]] = i;
        ccell[x][y]++;
if (ccell[x][y]>=size){
DEBUGPRINT("ERROR, it:%li FN2 x:%li y:%li ccell:%li i:%li nfg:%li\n", it, x, y, ccell[x][y], i, nfg);
DEBUGPRINT("maxx:%Le maxy:%Le minx:%Le miny:%Le\n", maxx, maxy, minx, miny);
}
    }
    // Done! The grains are distributed in the cells, now it's time to transform the cells in neighbours...

    ior[ncontot+1] = 0;
    iex[ncontot+1] = 0;
    icon = 1;
    ineigh = 0;

    // This will set boundary condition in x and y directions, but the rotine of detect contacts will cut it off, if itn't on the studied problem
    for (x = 0; x < Ll; x++){ // Getting the element on x position of the cell
        for (y = 0; y < Mm; y++){ // Getting the element on y positon of the cell
            for (k = -1; k < 2; k++){ // Looking to the neighbours cells in x direction
                for(h = 0; h < 2; h++){ // Looking to the neighbours cells in y direction
                    if (!((k == -1) && (h == 0))){ // Eliminating the repeated search
//DEBUGPRINT("h=%i k=%i x=%i y=%i\n",h,k,x,y);
//DEBUGPRINT("Ll=%i Mm=%i\n",Ll,Mm);
                        for (g = 0; g < ccell[x][y]; g++){ // Getting all the elements on the cell
                            i = cells[x][y][g]; // Element one of the neighbourhood
                            // Looking for boundary conditions and unnecessary searchs:
                            d = x+k;
                            d = d < 1 ? Ll -2 : d; // Condition to the left side, taking the last x position in cell
                            d = d >= Ll-1 ? 1 : d; // Condition to the right side, taking the first x position in cell
                            e = y+h;
                            e = e < 1 ? Mm -2 : e; // Condition to the bottom side, taking the higher y position in cell, as it is builded, it should never happen
                            e = e >= Mm-1 ? 1 : e; // Condition to the top side, taking the lower y position in cell
                            for (f = ((h == 0) && (k == 0)) ? g+1 : 0; f < ccell[d][e]; f++){ // Getting all non repeated elements of the neighbour and self cell
//DEBUGPRINT("d=%i e=%i f=%i g=%i i=%i j=%i\n",d,e,f,g,i,j);
                                j = cells[d][e][f]; // Element two of the neighbourhood

                                // As it was when find an elegible neighbour:
                                xij = distperio(rx[j] - rx[i]);
                                yij = ry[j] - ry[i];
                                dij = (xij*xij + yij*yij);
                                if (dij < x12){
                                    ineigh++;
//DEBUGPRINT("ineigh=%i\n",ineigh);
                                    listi[ineigh] = i;
                                    listj[ineigh] = j;
                                    if ((ior[icon]==i)&&(iex[icon]==j))//We know in advance that the next one previously in contact must remain in the neighbourhood
                                    {
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

    // Stupid way to do this, thinking it uses the neigbour, but for now I see no other solution, in a generic code...
    // Searching in walls on x:
    for (x = 0; x < Ll; x++){
        for (y = 0; y < 2; y++){
            for (g = 0; g < ccell[x][y]; g++){
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++)
                {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12)
                    {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j))
                        {
                            io[ineigh] = 1;
                            icon ++;
                        }
                        else
                            io[ineigh] = 0;
                    }
                }
            }
        }
        for (y = Mm-2; y < Mm; y++){
            for (g = 0; g < ccell[x][y]; g++){
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++)
                {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12)
                    {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j))
                        {
                            io[ineigh] = 1;
                            icon ++;
                        }
                        else
                            io[ineigh] = 0;
                    }
                }
            }
        }
    }
    // Searching in walls on y
    for (y = 2; y < Mm-2; y++){
        for (x = 0; x < 2; x++){
            for (g = 0; g < ccell[x][y]; g++){
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++)
                {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12)
                    {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j))
                        {
                            io[ineigh] = 1;
                            icon ++;
                        }
                        else
                            io[ineigh] = 0;
                    }
                }
            }
        }
        for (x = Ll-2; x < Ll; x++){
            for (g = 0; g < ccell[x][y]; g++){
                i = cells[x][y][g];
                for (j =nfg+1;j<=ngtot;j++)
                {
                    xij = distperio(rx[j] - rx[i]);
                    yij = ry[j] - ry[i];
                    dij = (xij*xij + yij*yij);
                    if (dij < x12)
                    {
                        ineigh++;
                        listi[ineigh] = i;
                        listj[ineigh] = j;
                        if ((ior[icon]==i)&&(iex[icon]==j))
                        {
                            io[ineigh] = 1;
                            icon ++;
                        }
                        else
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
double long distperio(double long x)
{
if (x>halfwidth)
     return (x-width);
else if (x<-halfwidth)
     return (x+width);
else return x;
}

/*********** detect contacts ***********/
/*********** require to know reactsave ***********/
void detectcontacts()
{
long i,j,il;
double long xij,yij,dij,radius2;

/*********** Initialization ***********/
ncontot = 0;
icont0 = 0;
for (i = 1;i<=nfg; i++)
    Zcont[i] = 0;

/*********** Contact between all grains ***********/
for (il=1;il<=nlt;il++) //loop over neighbouring grains
{
    i=listi[il]; // get back which i is in the neighbourhood of j>i
    j=listj[il];
    xij = distperio(rpx[j] - rpx[i]);/***! periodic boundary conditions along x***/
    yij = rpy[j] - rpy[i];
    dij = xij*xij + yij*yij;
    radius2 = (radius[i]+radius[j])*(radius[i]+radius[j]);
    if (io[il]==1) //flag containing the information that the neighbourhood il (=neighbours i and j labelled by il), was corresponding at the previous time step to grains in contact
        icont0++;/***! index of previous contacts to get back the tangential force of previous time-step**/
    if (dij<radius2)
    {
        Zcont[i] = Zcont[i] + 1;
        Zcont[j] = Zcont[j] + 1;
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
    }
    else
        io[il] = 0; // This neighbourhood il (neighbouring grains i,j) does NOT correspond to a contact
	if (il==nl)
		ncontfree=ncontot;
}
return;
}

/*********** Force calculation ***********/
void forcecalculation()
{
long i,j,il;
double long xn,yn,xt,yt;
double long fn,ft,fx,fy,vijn,vijt,ftest;

/*********** initialize*********/
if (menu<15)
{
    for (i =1;i<=nfg;i++)
    {
        fpx[i] = 4.0L*radius[i]*radius[i]*sintheta;
        fpy[i] = - 4.0L*radius[i]*radius[i]*costheta; // Gravity
        gam[i] = 0.0L;
    }
}

else if (menu==25)
{
	for (i =1;i<=nfg;i++)
	{
		fpx[i] = 0.0L;
		gam[i] = 0.0L;
		if ((rpy[i]<bottom+6.0L)&&(rpy[i]>bottom+3.0L))
			fpy[i] = -4.0L*radius[i]*radius[i]*gravity*((bottom+6.0L)-rpy[i])*(rpy[i]-(bottom+3.0L));
		else if ((rpy[i]>height-6.0L)&&(rpy[i]<height-3.0L))
			fpy[i] = 4.0L*radius[i]*radius[i]*gravity*((height-6.0L)-rpy[i])*(rpy[i]-(height-3.0L));
		else fpy[i] = 0.0L;
	}
}
else if (menu==26)
{
	for (i =1;i<=nfg;i++)
	{
		fpx[i] = 0.0L;
		fpy[i] = 4.0L*radius[i]*radius[i]*gravity;
		gam[i] = 0.0L;
	}
}
#ifdef FLUID
else if ((menu >= 30) && (menu < 45))
{
    for (i =1;i<=nfg;i++){
        fpx[i] = 0.0L;
        fpy[i] = -gravity/ivmass[i]; // Gravity, with density = 1;
        gam[i] = 0.0L;
    }
}
#endif
else
{
    for (i =1;i<=nfg;i++)
    {
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
for (il = 1;il<= ncontot;il++)
{
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
    if (fabsl(ft)>=ftest)
    {
        if ((Zcont[i]!=1)&&(Zcont[j]!=1))
            ncontglis=ncontglis+1;
        if (ft>0)
            ft = ftest;
        else
            ft = -ftest;
    }

	/** TEST TEST**/
#ifdef debug
	{
		if (i==71)
			fprintf(outest,"%Le\t%ld\t%ld\t%Le\t%Le\n",(double long)it*dt,j,il,react[il],ft);
		if (j==71)
			fprintf(outest,"%Le\t%ld\t%ld\t%Le\t%Le\n",(double long)it*dt,i,il,-react[il],-ft);
	}
#endif
/** TEST TEST**/

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
	if (il>ncontfree)
	{
		if (j>(nfg+nwall))
		{
			fxup+=fx;
			fyup+=fy;
			springxx+=kn*xn*xn;
			springyy+=kn*yn*yn;
			springxy+=kn*xn*yn;
		}
		else
		{
			fxdown+=fx;
			fydown+=fy;
		}
	}
}

// Adding the fluid force on each free grain
#ifdef FLUID
    if ((menu >= 30) && (menu < 45)){
        // This function calculate the profile of the fluid in function of the particles, especialy to the packing fraction
        particlesprofile();

        for (i = 1; i <= ngtot; i++){
            fluidforce(i);
        }
    }
#endif

/********** Add wall friction ***********/
/* friction opposed to v[i] */
    if (menu==2)
	{
		for (i =1;i<=nfg;i++)
		{
			ftest=-muwall*sumfn[i]/(radius[i]*sqrt(1.0E-14L+vpx[i]*vpx[i]+vpy[i]*vpy[i]));
			fpx[i] += vpx[i]*ftest;
			fpy[i] += vpy[i]*ftest;
		}
	}

/* friction opposed to <v[i]> i.e. along -X */
    if (menu==3)
    {
        for (i =1;i<=nfg;i++)
        {
            fpx[i] -= muwall*sumfn[i]/radius[i];
        }
    }

/* friction opposed to <v[i]> i.e. along -X and proportionnal to fit of P */
    if (menu==4)
    {
        for (i =1;i<=nfg;i++)
        {

            fpx[i] -= muwall*(2.0L*M_PI)*(aa+bb*rpy[i]);// using a first fit on syy/sphi from a previous run.(2.0L*M_PI)*(37.172L-1.0696L*rpy[i])
        }
    }

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
// This rotine solve the partial differential equation for the fluid, based on the mesh over z and the dt
// Remembering that tau is on the tick of the discretizaton of z and u is between the tick
void solvefluidvelocity(){
    long i, j;

    double long phib = 0.85L, iphi = 0.0L, vps = 0.0L, vp = 0.0L;
    // Definition of iphi = integral of phi(pakingfraction) along z:
    for (i = 0; i < NSTAT; i++){
        iphi += packfraction[i];
    }

    // The ratio of vp²/vps gives the average number of moving grains
    for (i = 1; i <= nfg; i++){
        vps += vx[i]*vx[i];
        vp += vx[i];
    }

    // For the log profile, the multipling coefficient should be smaller than 0.9 to get a stable zb, otherwise, zb goes to 0
    if (vps > 0.01L){
        if (0.9*vp*vp/(vps*width) < iphi)
            zb = iphi -0.9*vp*vp/(vps*width);
        else
            zb = 0.0L;
    } else {
        zb = iphi;
    }

    // Upadte the fluid
    long zbi = (long)(zb/dz)+LISTAT-1;
//DEBUGPRINT("i:%li zbi:%li zb:%Le iphi:%Le",it, zbi, zb, iphi);
//getchar();
    // Condition to the stability of this integration:
    // dt < (dz)^2/2, for the viscous case
    // dt < (dz)^2/(2*SystemSize), for the turbulent case
    double long l, up[NSTAT], reynoldsd = 26.0L; // l is the mixing lenght, tau is the shear per density of the fluid
    double long fluid_density = grain_density/density_ratio;

    // Implementing Van Driest's equation:
    tau[0] = viscosity*u[0+1]/dz;
    for (j = 0+1; j < NSTAT-1; j++){
        if(j-zbi < 0)
            l = 0.0L;
        else
            l = kapa*(j-zbi)*dz*(1.0L -exp(-(j-zbi)*dz*ustar/(reynoldsd*viscosity)));
        tau[j+1] = (viscosity+ l*l*fabs(u[j+1] -u[j])/dz)*(u[j+1] -u[j])/(dz);
        up[j] = u[j] +(tau[j+1]-tau[j])*dt/dz -dt*fxF[j]/(fluid_density*(1.0L-packfraction[j]));
        dpressuredz[j+1] = -fluid_density*gravity -fyF[j]/((1.0L-packfraction[j]));
if(packfraction[j] >= 1.0L){
    if(j < LISTAT +1.0L/dz){
        up[j] = 0.0L;
        tau[j] = ustar*ustar;
    } else {
        DEBUGPRINT("ERROR, it %li phi[%li] is %Le\n", it, j, packfraction[j]);
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
DEBUGPRINT("ERROR, it %li up[%li] is INF, tau+1 %Le, tau %Le, fxF %Le, zbi %li, zb %Le\n", it, j, tau[j+1], tau[j], fxF[j], zbi, zb);
getchar();
}
//printf("l %Le kapa %Le tau[%i] %Le tau[%i] %Le u[%i] %Le\n", l, kapa, j+1, tau[j+1], j, tau[j], i+1, u[i+1][j]);
//getchar();
    }
    tau[NSTAT] = ustar*ustar;
    up[NSTAT-1] = u[NSTAT-1] +(tau[NSTAT]-tau[NSTAT-1])*dt/dz;
    dpressuredz[NSTAT] = -fluid_density*gravity -fyF[NSTAT-1]/(1.0L-packfraction[NSTAT-1]);

    if (it % (itimemeasurements*d) == 0){
        d *= 2;
        for (j = 0; j < NSTAT; j=1.1*j+1){
//DEBUGPRINT("viscosity[%i] %Le\n",it, viscosity);
            fprintf(out_fluid, "%Le %Le %Le %Le %Le %Le %Le\n", it*dt, (j-LISTAT)*dz, u[j], tau[j+1], dpressuredz[j+1], fxF[j]/(1.0L-packfraction[j]), fyF[j]/(1.0L-packfraction[j]));
        }
    }
long double ftmpx = 0.0L, ftmpy = 0.0L;
    for (j =  0; j < NSTAT; j++){
        u[j] = up[j];
        ftmpx += fxF[j]*width*dz / (packfraction[j] != 0.0L ? packfraction[j] : 1.0L);
        ftmpy += fyF[j]*width*dz / (packfraction[j] != 0.0L ? packfraction[j] : 1.0L);
        packfraction[j] = 0.0L;
        fxF[j] = 0.0L;
        fyF[j] = 0.0L;
    }
    for (j = 1; j <= ngtot; j++){
        ftmpx -= fxfluid[j];
        ftmpy -= fyfluid[j];
    }
//if (it%1000==0){
//DEBUGPRINT("Momentum balance %Le %Le\r", ftmpx, ftmpy);
//}

    return;
}

// This rotine will find the bed position
void particlesprofile(){
//    if(it % ifreqv == 0)
    {
        // For now, this for is across the z(y) direction
        cumulatestats();
    }

    return;
}

// The implementation of the fluid force on the grain. This rotine should be used for each grain that have a fluid interaction. One important parameter to this calcul is the zb.
void fluidforce(long i){
    // This function have 2 diferent origins: a drag component and an Archimed component

    // Calculating the fluid velocity
    double long fdx, fdy, fax, fay, vmodulus, fluidvelocity = 0.0L; //u; // u (velocity of the fluid) can be function of x component, but also from z (or y), in the future.
    double long fluid_density = grain_density/density_ratio;

    long index_top = LISTAT+ceil(ry[i]/dz), index_bottom = LISTAT+floor(ry[i]/dz);
    fluidvelocity = index_top == index_bottom ? u[index_top] : u[index_bottom]+(u[index_top]-u[index_bottom])*(ry[i]-(index_bottom-LISTAT)*dz)/((index_top-index_bottom)*dz);
    vmodulus = sqrt((fluidvelocity-vx[i])*(fluidvelocity-vx[i]) + (vy[i])*(vy[i])); // Modulus of the diference of velocity between grain and fluid
    reynoldsg = 2.0L*vmodulus*radius[i]/viscosity; // Reynolds number of the grain interacting with the fluid
    // This if takes care of a divergency of dragc -> 0
    if (((reynoldsg >= 1e-18L) && !(isnan(reynoldsg) || isinf(reynoldsg))) && (vmodulus > 1e-18L)){
        dragc = (sqrt(dragd)+sqrt(reynoldsc/reynoldsg)); // Calculating (drag coefficient)^(1/2)
        dragc *= dragc; // Drag Coefficient

        /// Calculating drag force
//        fdrag = PI/2.0L*radius[i]*radius[i]*dragc*fluid_density*(vmoduls*vmodulus); // Modulus of the force
        fdx = PI/2.0L*radius[i]*radius[i]*dragc*fluid_density*(vmodulus*(fluidvelocity-vx[i])); // Cosine of the vector
        fdy = PI/2.0L*radius[i]*radius[i]*dragc*fluid_density*(vmodulus*(-vy[i])); // Sine of the vector
    } else {
//        fdrag = PI/2.0L*radius[i]*radius[i]*fluid_density*vmodulus;
        fdx = PI/2.0L*radius[i]*reynoldsc*viscosity*fluid_density*(fluidvelocity-vx[i]); // Cosine of the vector
        fdy = PI/2.0L*radius[i]*reynoldsc*viscosity*fluid_density*(-vy[i]); // Sine of the vector
    }
//printf("\n");
    /// Archimedes force
    // Archimedes force has 2 components: against gravity and in the shearing direction
    // For geometry reasons, these expression uses the density ratio and the mass of the particle, if geometry changes, this change will be automatically done here
    // Calculating the hydrodynamic force
    index_top = LISTAT+ceil(ry[i]/dz-0.5L);
    index_bottom = LISTAT+floor(ry[i]/dz-0.5L);
    fax = index_top == index_bottom ? (tau[index_top+1] -tau[index_bottom-1])/(2.0L*dz) : (tau[index_top] -tau[index_bottom])/(dz);
    fax *= (4.0L*PI*radius[i]*radius[i]*radius[i])/(3.0L);
    fay = index_top == index_bottom ? dpressuredz[index_top] : dpressuredz[index_bottom]+(dpressuredz[index_top]-dpressuredz[index_bottom])*(ry[i]-(index_bottom-LISTAT)*dz)/((index_top-index_bottom)*dz);
    fay *= -(4.0L*PI*radius[i]*radius[i]*radius[i])/(3.0L);
    // Calculating the hydrostatic force
//    fay = gravity/(ivmass[i] * density_ratio); // Modulus of the force, oriented agains gravity

    // Force caused by grains on the fluid, only a drag contribution:
    fxfluid[i] = (fdx);
    fyfluid[i] = (fdy);

    // An important hint: The fluid it self have a time constant to ajust the drag force. If the fluid time is greater than the dt, divernces may happen.
    // Force caused by fluid on grains:
    fpx[i] += fxfluid[i] +fax;
    fpy[i] += fyfluid[i] +fay;

    return;
}
#endif

/*********** Predictor (require ax, ay, domedt)  ***********/
void predictor()
{
long i;

if (menu==15)
{
	hp=height+dt*dhdt+halfdt2*ddhdtt;
	dhdtp=dhdt+dt*ddhdtt;
	ddhdttp=ddhdtt;
	wallx+=dt*wallv;
	if (wallx>width)
		wallx-=width;
	for (i=1;i<=nwall;i++)
	{
		rpx[nfg+nwall+i]=distperio(xwall[i]+wallx);
		rpy[nfg+nwall+i]=ywall[i]+hp;
		vpy[nfg+nwall+i]=dhdtp;
	}
}
else if ((menu>15)&&(menu<29))
{
	if (menu==27)
	{
		wallv=wallvzero*exp(-(double long)it*timewall);
		for (i=1;i<=nwall;i++)
			vpx[nfg+nwall+i]=wallv;
	}
	wallx+=dt*wallv;
	if (wallx>width)
		wallx-=width;
	for (i=1;i<=nwall;i++)
	{
		rpx[nfg+nwall+i]=distperio(xwall[i]+wallx);
		rpy[nfg+nwall+i]=ywall[i]+height;//useless?
	}
}
else if (menu==29)
{
	hp=height+dt*dhdt+halfdt2*ddhdtt;
	dhdtp=dhdt+dt*ddhdtt;
	ddhdttp=ddhdtt;
	wallxp=wallx+dt*wallv+halfdt2*walla;
	if (wallxp>0.5L*width)
		wallxp-=width;
	wallvp=wallv+dt*walla;
	wallap=walla;
	for (i=1;i<=nwall;i++)
	{
		rpx[nfg+nwall+i]=distperio(xwall[i]+wallxp);
		rpy[nfg+nwall+i]=ywall[i]+hp;
		vpx[nfg+nwall+i]=wallvp;
		vpy[nfg+nwall+i]=dhdtp;
	}
}

/*********** Predictions sur les grains libres ***********/
for (i=1;i<=nfg;i++)
{
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
void corrector()
{
long i;

if (menu==15)
{
    ddhdtt=fyup-1.0L-gn*dhdtp;
    dhdt = dhdtp + halfdt*(ddhdtt-ddhdttp);
    height = hp;
    for (i=1;i<=nwall;i++)
    {
        rx[nfg+nwall+i]=distperio(xwall[i]+wallx);
        ry[nfg+nwall+i]=ywall[i]+height;
    }
}
else if ((menu>15)&&(menu<29))
{
	if (menu<28)
		height+=(fyup-1.0L)/springyy;
	for (i=1;i<=nwall;i++)
	{
		rx[nfg+nwall+i]=distperio(xwall[i]+wallx);
		ry[nfg+nwall+i]=ywall[i]+height;
	}
}
else if (menu==29)
{
	ddhdtt=fyup-1.0L-gn*dhdtp;
	dhdt = dhdtp + halfdt*(ddhdtt-ddhdttp);
	height = hp;
	walla=shear+fxup;
	wallv = wallvp + halfdt*(walla-wallap);
	wallx = wallxp;
	for (i=1;i<=nwall;i++)
	{
		rx[nfg+nwall+i]=distperio(xwall[i]+wallx);
		ry[nfg+nwall+i]=ywall[i]+height;
	}
}

/*********** Correction on free grains ***********/
for (i = 1;i<=nfg;i++)
{
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

/*********** Correction on free grains ***********/
void cumulatestats()
{
long i,j,k,il;
double long ratio,fx,fy;  // weight function
// bin size =0.1d;

/*********** initialize*********/
for (i =1;i<=nfg;i++)
    Zglis[i] = 0;
#ifdef FLUID
    for(k = 0; k < NSTAT; k++){
        instPhi[k] = 0.0L;
        instVx[k] = 0.0L;
        instVy[k] = 0.0L;
        instP[k] = 0.0L;
        instZcont[k] = 0.0L;
        instZglis[k] = 0.0L;
        instFx[k] = 0.0L;
        instFy[k] = 0.0L;
        instFn[k] = 0.0L;
        instFt[k] = 0.0L;
        instInterp[k] = 0.0L;
        packfraction[k] = 0.0L;
        fxF[k] = 0.0L;
        fyF[k] = 0.0L;
    }
#endif

/*********** contact between free grains ***********/
for (il = 1;il<= ncontfree;il++)
    if (fabsl(react[il])>=(0.99L*reacn[il]))
        {
            Zglis[ior[il]]+=1;
            Zglis[iex[il]]+=1;
        }
for (il=ncontfree+1;il<=ncontot;il++)
    if (fabsl(react[il])>=(0.99L*reacn[il]))
    {
        Zglis[ior[il]]+=1;
    }

for( il=1;il<=ncontot;il++)
{
    fx = reacn[il]*xnij[il] - react[il]*ynij[il];
    fy = reacn[il]*ynij[il] + react[il]*xnij[il];
    i = ior[il];
    j = iex[il];
    if (ry[i]>ry[j])
        for (k=ceil(ry[j]/dz); k<=floor(ry[i]/dz); k++)
        {
            statFx[LISTAT+k]+=fx;
            statFy[LISTAT+k]-=fy;
#ifdef FLUID
            instFx[LISTAT+k] += fx;
            instFy[LISTAT+k] -= fy;
            instFn[LISTAT+k] += reacn[il];
            instFt[LISTAT+k] += fabs(react[il]);
            instInterp[LISTAT+k] -= eij[il];
#endif
        }
    else
        for (k=ceil(ry[i]/dz); k<=floor(ry[j]/dz); k++)
        {
            statFx[LISTAT+k]-=fx;
            statFy[LISTAT+k]+=fy;
#ifdef FLUID
            instFx[LISTAT+k] -= fx;
            instFy[LISTAT+k] += fy;
            instFn[LISTAT+k] += reacn[il];
            instFt[LISTAT+k] += fabs(react[il]);
            instInterp[LISTAT+k] -= eij[il];
#endif
        }
}
ncumul+=width;
    for (i=1;i<=ngtot;i++)
    {
        for (k=floor((ry[i]-radius[i])/dz); k<=ceil((ry[i]+radius[i])/dz)-1; k++){
            if ((radius[i] > fabsl(ry[i]-(double long)(k+1)*dz)) && (radius[i] > fabsl(ry[i]-(double long)k*dz))){
                ratio = (sqrt(radius[i]*radius[i]-pow((ry[i]-(double long)(k+1)*dz),2.0L))*(dz*(double long)(k+1)-ry[i])-sqrt(radius[i]*radius[i]-pow((ry[i]-(double long)k*dz),2.0L))*(dz*(double long)k-ry[i]))+radius[i]*radius[i]*(asin(((double long)(k+1)*dz-ry[i])/radius[i])-asin(((double long)k*dz-ry[i])/radius[i]));
            } else if (radius[i] > fabsl(ry[i]-(double long)(k+1)*dz)){
                ratio = (sqrt(radius[i]*radius[i]-pow((ry[i]-(double long)(k+1)*dz),2.0L))*(dz*(double long)(k+1)-ry[i]))+radius[i]*radius[i]*(asin(((double long)(k+1)*dz-ry[i])/radius[i])+M_PI*0.5L);
            } else if (radius[i] > fabsl(ry[i]-(double long)k*dz)){
                ratio = (-sqrt(radius[i]*radius[i]-pow((ry[i]-(double long)k*dz),2.0L))*(dz*(double long)k-ry[i]))-radius[i]*radius[i]*(asin(((double long)k*dz-ry[i])/radius[i])-M_PI*0.5L);
            } else {
                ratio = 0.0L;
            }
#ifdef FLUID
            // To be sure that k is not invalidating any data:
            if ((k >= -LISTAT) && (k < NSTAT-LISTAT)){
#endif
                statPhi[LISTAT+k] += ratio;
                statVx[LISTAT+k] += ratio*vx[i];
                statVy[LISTAT+k] += ratio*vy[i];
                statP[LISTAT+k]+=ratio*sumfn[i]/(2.0L*M_PI*radius[i]);
                statZcont[LISTAT+k] += ratio*Zcont[i];
                statZglis[LISTAT+k] += ratio*Zglis[i];
#ifdef FLUID
                instPhi[LISTAT+k] += ratio;
                instVx[LISTAT+k] += ratio*vx[i];
                instVy[LISTAT+k] += ratio*vy[i];
                instP[LISTAT+k] += ratio*sumfn[i]/(2.0*PI*radius[i]);
                instZcont[LISTAT+k] += ratio*Zcont[i];
                instZglis[LISTAT+k] += ratio*Zglis[i];
if(isinf(ratio)){
DEBUGPRINT("ERROR, it %li ratio of particle[%li] on region[%li] is INF\n", it, i, k);
getchar();
}
                packfraction[LISTAT+k] += ratio/(width*dz);
if(isinf(packfraction[LISTAT+k])){
DEBUGPRINT("ERROR, it %li phi[%li] is INF, width %Le and dz %Le\n", it, LISTAT+k, width, dz);
getchar();
}
if((packfraction[LISTAT+k] >= 1.0L)&&(k>LISTAT+1.0L/dz)){
DEBUGPRINT("ERROR, it %li phi[%li] is %Le\n", it, LISTAT+k, packfraction[LISTAT+k]);
getchar();
}
                fxF[LISTAT+k] += ratio/(width*dz)*fxfluid[i]/(4.0L*PI*radius[i]*radius[i]);
                fyF[LISTAT+k] += ratio/(width*dz)*fyfluid[i]/(4.0L*PI*radius[i]*radius[i]);
            } else {
                fprintf(stderr, "Error, invalid profile of particle %li, z: %li %Le, it:%li\n", i, k, ry[i]-radius[i], it);
                getchar();
            }
#endif
        }
    }
}

/*********** timemeasurement ***********/
void timemeasurement()
{
    double long cineti,rotational,elastic,gravitational;
    long i,nsi,nsi1;
    long il;

    /*********** Test des violations ( = profondeurs de penetration) ***********/
    /*********** Maximal interpenetration distance ***********/
    viol = 0.0;
    for( il = 1;il<= ncontot;il++)
    if (eij[il]<viol)
    viol=eij[il];
    viol = -0.5L*viol;

    /*********** Energies cinetique et elastictielle totales ***********/
    cineti = 0.0L;
    rotational = 0.0L;
    elastic = 0.0L;
    gravitational = 0.0L;
    for (i = 1;i<=nfg;i++)
    cineti += 0.5L*(vx[i]*vx[i] + vy[i]*vy[i])/ivmass[i];
    for (i = 1; i <= nfg; i++)
    rotational += 0.5L*(ome[i]*ome[i])/ivmomine[i];
    for (i = 1; i <= nfg; i++)
    gravitational += ry[i]/ivmass[i];
    for (il = 1;il<=ncontot;il++)
    elastic += 0.5L*kn*eij[il]*eij[il];

    nsi = 0;
    nsi1 = 0;
    for( i=1;i<=nfg;i++)
    {
        if (Zcont[i]>0)
        nsi = nsi + 1;
        if (Zcont[i]<=1)
        nsi1 = nsi1 + 1;
    }
    fprintf(out,"%Le\t%ld\t%ld\t%ld\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%ld\n", (double long)it*dt,ncontot,nsi,nsi1,viol,cineti,elastic,rotational,gravitational,-fydown,fxdown,fyup,-fxup,height,wallx,ncontglis);
    fprintf(out1,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %ld %ld %ld\n", (double long)it*dt,cineti,rotational,elastic,gravitational,-fydown,fxdown,fyup,-fxup,height,wallx,ncontglis,ncontot,nlt);

    return;
}

/****************period save***********************/
void periodsave()
{
    FILE *out2;
#ifdef FLUID
    savestat();
    savefluidconfiguration();
#endif
    sprintf(periodsaveout, "save%ld.bin",nperiodsave);
    out2=fopen(periodsaveout,"w");
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
    fwrite(&ncontfree,sizeof(ncontfree),1,out2);
    fwrite(&ncontot,sizeof(ncontot),1,out2);
    fwrite(ior,sizeof(long),ncontot+1,out2);
    fwrite(iex,sizeof(long),ncontot+1,out2);
    fwrite(react,sizeof(double long),ncontot+1,out2);
    fflush(out2);
    fclose(out2);
    nperiodsave++;
    return;
}

#ifdef FLUID
void savestat(){
    long i;
    char statname[200];
    FILE *out2;

    sprintf(statname, "stat%li.igr", nperiodsave);
    if ((out2 = fopen (statname,"w")) == NULL)
    {
        fprintf(stderr,"Probleme d'ouverture de %s\n",statname);
        exit (1);
    }
    fprintf(out2,"IGOR\n\n");
    fprintf(out2,"WAVES/O/R/S sz\tsphi\tsvx\tsyy\tsxy\tpre\tsZt\tsZg\n");
    fprintf(out2,"BEGIN\n");
    for (i = 0; i < NSTAT; i++)
        fprintf(out2,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", (double long)(i-LISTAT)*dz,statPhi[i]/((double long)ncumul*dz),statVx[i]/statPhi[i],statFy[i]/(double long)ncumul,statFx[i]/(double long)ncumul,statP[i]/((double long)ncumul*dz), statZcont[i]/statPhi[i],statZglis[i]/statPhi[i]);

        //Used to be: statP[i]/statPhi[i]
    fprintf(out2,"END\n\n");

    fflush(out2);
    fclose(out2);

    sprintf(statname, "inst%li.igr", nperiodsave);
    if ((out2 = fopen (statname,"w")) == NULL)
    {
        fprintf(stderr,"Probleme d'ouverture de %s\n",statname);
        exit (1);
    }
    fprintf(out2,"#sz\tsphi\tsvx\tpre\tsZt\tsZg\tsn\tst\tInterp\n");
    for (i = 0; i < NSTAT; i++)
        fprintf(out2,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", (double long)(i-LISTAT)*dz,instPhi[i]/(width*dz),instVx[i]/instPhi[i],instP[i]/(width*dz), instZcont[i]/instPhi[i],instZglis[i]/instPhi[i],instFn[i]/width,instFt[i]/width,instInterp[i]);

    fflush(out2);
    fclose(out2);
}
#endif

/*********** Save Configuration ***********/
void saveconfiguration()
{
long i;
FILE *out2;
FILE *out3;

if ((out2 = fopen (nomout2,"w")) == NULL)
{
    fprintf(stderr,"Probleme d'ouverture de %s\n",nomout2);
    exit (1);
}
fprintf(out2,"IGOR\n\n");
fprintf(out2,"WAVES/O/R/S sz\tsphi\tsvx\tsyy\tsxy\tpre\tsZt\tsZg\n");
fprintf(out2,"BEGIN\n");
for (i = 0; i < NSTAT; i++)
    fprintf(out2,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", dz*(double long)(i-LISTAT),statPhi[i]/((double long)ncumul*dz),statVx[i]/statPhi[i],statFy[i]/(double long)ncumul,statFx[i]/(double long)ncumul,statP[i]/((double long)ncumul*dz), statZcont[i]/statPhi[i],statZglis[i]/statPhi[i]);

	//Used to be: statP[i]/statPhi[i]
fprintf(out2,"END\n\n");

fflush(out2);
fclose(out2);

out3=fopen(datanameout,"w");

	fwrite(&ngtot,sizeof(ngtot),1,out3);
	fwrite(radius,sizeof(double long),ngtot+1,out3);
	fwrite(rx,sizeof(double long),ngtot+1,out3);
	fwrite(ry,sizeof(double long),ngtot+1,out3);
	fwrite(vx,sizeof(double long),ngtot+1,out3);
	fwrite(vy,sizeof(double long),ngtot+1,out3);
	fwrite(ome,sizeof(double long),ngtot+1,out3);
	fwrite(ax,sizeof(double long),ngtot+1,out3);
	fwrite(ay,sizeof(double long),ngtot+1,out3);
	fwrite(domedt,sizeof(double long),ngtot+1,out3);
	fwrite(&ncontfree,sizeof(ncontfree),1,out3);
	fwrite(&ncontot,sizeof(ncontot),1,out3);
	fwrite(ior,sizeof(double long),ncontot+1,out3);
	fwrite(iex,sizeof(double long),ncontot+1,out3);
	fwrite(react,sizeof(double long),ncontot+1,out3);

fflush(out3);
fclose(out3);

return;
}

void savefluidconfiguration(){
    long j = 0;
    char statname[200];

    sprintf(statname, "fluid%li.igr", nperiodsave);
    FILE *out_fluid1;
    if ((out_fluid1 = fopen (statname,"w")) == NULL) {
        fprintf(stderr,"Probleme d'ouverture de %s\n",statname);
        exit (1);
    }
#ifdef FLUID
    fprintf(out_fluid1, "#z u(z) tau(z) dpress(z)/dz fxF(z)/phi(z) fyF(z)/phi(z)\n");
    for (j = 0; j < NSTAT; j++){
        fprintf(out_fluid1, "%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", (j-LISTAT)*dz, u[j], tau[j+1], dpressuredz[j+1], fxF[j]/(1.0L-packfraction[j]), fyF[j]/(1.0L-packfraction[j]));
    }
    fclose(out_fluid1);
#endif
}

/*********** Close and save final Igor file ***********/
void closeall()
{
long i;

fprintf(out,"END\n\n");
fprintf(out,"WAVES/O/R/S rr\txx\tyy\tvx\tvy\tome\n");
fprintf(out,"BEGIN\n");
for (i = 1; i<=ngtot; i++)
    fprintf(out,"%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n", radius[i],rx[i], ry[i],vx[i], vy[i],ome[i]);
fprintf(out,"END\n\n");
fflush(out);
fclose(out);

/** TEST TEST**/
#ifdef debug
	{
		fprintf(outest,"END\n\n");
		fflush(outest);
		fclose(outest);
	}
#endif
/** TEST TEST**/

fprintf(out1,"\n\nEn sortie, ngtot=%ld\n",ngtot);
fprintf(out1," ncontfree=%ld  et    %ld \n",ncontfree,ncontot);

#ifdef FLUID
    if ((menu >= 30) && (menu < 45))
        fclose(out_fluid);
#endif

return;
}

/*********** Main ***********/
/*********** Main ***********/
/*********** Main ***********/
/*********** Main ***********/
int main ()
{
init();
findneighbours2();
predictor();
detectcontacts();
forcecalculation();
timemeasurement();
periodsave();
fprintf(out1,"\n\nEn entrée test, ncontfree=%ld  et    %ld \n",ncontfree,ncontot);

it=0;
endreached=1;
while(endreached!=0)
{
    it++;
    if ((it%ifreqv)==0)
        findneighbours2();
    predictor();
    detectcontacts();
    forcecalculation();
    corrector();
    if ((it%itimemeasurements)==0)
        timemeasurement();
    if ((it%iperiodsave)==0)
        periodsave();
#ifdef FLUID
    if ((menu > 30) && (menu < 45))
        solvefluidvelocity();
#endif
#ifndef FLUID
    // This rotine will happen each timestep in the fluid code
    if ((it%icustats)==0)
        cumulatestats();
#endif
    if ((it%isave)==0){
        saveconfiguration(); //All state variables are saved, including tangential forces
    }
    if (it>nbiteration){         //nbr pas de temps (en unité de dt)
        long i = 0;
        long double vel = 0.0L;
        for (i = 1; i < nfg; i++){
            vel += sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
        }
        if (vel < 0.001L){
            endreached=0;
        }
    }
}
if (extract==1)
{
	fprintf(out1,"extract wall \n");
	extractwall();
	findneighbours2();
	predictor();
	detectcontacts();
	forcecalculation();
	fprintf(out1,"free grain number=%ld \n",nfg);
	fprintf(out1,"number of grains in each wall=%ld \n",nwall);
}
fprintf(out1,"'end reached !'\n");
findneighbours2();
predictor();
detectcontacts();
forcecalculation();
timemeasurement();
periodsave();
saveconfiguration();
#ifdef FLUID
    if((menu > 30) && (menu < 45))
        savefluidconfiguration(); //Save all variables of the fluid
#endif
fprintf(out1,"\n\nEn sortie test, ncontfree=%ld  et    %ld \n",ncontfree,ncontot);

closeall();
return 0;
}

/*********** Complementary functions ***********/
/*********** Complementary functions ***********/
/*********** Complementary functions ***********/
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
    if(idum < 0) idum+=IM1;
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

double gaussian2(long *idum)
{
    double long magni1,phiga1,phiga2;
    double long fxT,fzT;
    double long x1,x2;
    long islu;

    x1 = ran2(&islu);
    x2 = ran2(&islu);

    magni1 = sqrt(-2.0L*log(x1));
    phiga1 = sin(2.0L*M_PI*x2);
    phiga2 = cos(2.0L*M_PI*x2);
    fxT = magni1*phiga1;
    fzT = magni1*phiga2;
    return 0.0L;
}

//
////----------------Compute Local Gamma--------------------------//
//double
//alpha,denumnonplastic,denumplastic,numplastic,numnonplastic,denum,num;
//double distancez;
//double gauss,mean,mean2,varianceGamma;
//denum=0.0;
//num=0.0;
//alpha=0.25;
//k=0;
//mean=0.0;
//mean2=0.0;
//
//for(int i=0;i<Nparticules; i++)
//{
//	for(int j=0;j<Nparticules; j++)
//	{
//		if (Ep[j]!=-16 && Ep[j]!=-15) // to avoid particles on the walls
//		{
//			distancez=(z[i]-z[j])*(z[i]-z[j]); // Compute the distance z
//			component
//
//			distance=distancez+(x[i]-x[j])*(x[i]-x[j]); // Compute the
//			distance
//
//			gauss=exp(-distance*alpha); // Compute Gaussienne
//
//			gama[i]+=(u[i]-u[j])*(z[i]-z[j])*gauss; // Numerator of Gamma
//
//			norm[i]+=distancez*gauss;   // Denomenator: normalisation
//
//		}else
//		{
//			gama[i]=0.0;
//			norm[i]=1.0;
//		}
//	}
//	//------- Compute Mean of Gamma in the center of the cell between 5
//	and 45-----------//
//	if (z[i] >=5 && z[i] <= 45)
//	{
//		k+=1;
//		moyenne+=-gama[i]/(norm[i]*dt) ;
//
//		moyenne2+=(((gama[i]/(norm[i]*dt)))*((gama[i]/(norm[i]*dt))));
//
//		moyenne4+= pow((gama[i]/(norm[i]*dt)),4.);
//
//
//	}
//	//------- End Mean of Gamma over all the system-----------//
//
//
//	if (Ep[i]!=-16 && Ep[i]!=-15){
//		if (int(t)%dt==0)
//		{
//			dataout  << gama[i]/(norm[i]*dt*meangama)  << " " <<
//			(gama[i]/(norm[i]*dt*meangama)*gama[i]/(norm[i]*dt*meangama))
//			<< endl;
//
//
//
//			norm[i]=0.0;
//			gama[i]=0.0;
//
//		}}
//}
//
////———Calcul of non gaussian parameter -----------//
//if (int(t)%dt==0)
//{
//	moyennecumul +=moyenne/k;
//	moyenne2cumul=moyenne2/k;
//	moyenne4cumul=moyenne4/k;
//	varianceGamma= moyenne2cumul-(moyenne/k)*(moyenne/k);
//	nongaussian=3./5.*moyenne4cumul/(moyenne2cumul*moyenne2cumul)-1.0;
//	nongaussiancumul+=nongaussian;
//	gamaout << tindex << " " << moyennecumul/tindex << " " << moyenne/k
//	<<  " " << varianceGamma  << " " << nongaussian << " " <<
//	nongaussiancumul/tindex << endl;
//
//
//
//
//	k=0;
//	moyenne=0.0;
//	moyenne2=0.0;
//	moyenne4=0.0;
//	denum=0.0;
//	num=0.0;
//	varianceGamma=0;
//
//}
