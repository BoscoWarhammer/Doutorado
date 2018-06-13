/**************************************************************************************************
    Runge-Kutta of 4th order to solve equations 
**************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double RvD = 26.0, k = 0.4, ustar = 0.1, nu = 0.01;

// Write the differential equation in function of dy/dx = z: Called F
double F(double x, double y, double z){
    double r = 0.0;
    r = z;
    return r;
}

// Write the differential equation in function of dz/dx, with z = dy/dx: Called G
double G(double x, double y, double z){
    double r = 0.0;
    //r = ;
    return r;
}

// Runge-Kutta method for 4th order in 2nd ODE
void RK4(double x, double &y, double &z, double h){
    double f = 0.0, g = 0.0;
    double k1 = 0.0, k2 = 0.0, k3 = 0.0, k4 = 0.0;
    double l1 = 0.0, l2 = 0.0, l3 = 0.0, l4 = 0.0;
    k1 = F(x, y, z);
    l1 = G(x, y, z);
    k2 = F(x +h/2.0, y +k1*h/2.0, z +l1*h/2.0);
    l2 = G(x +h/2.0, y +k1*h/2.0, z +l1*h/2.0);
    k3 = F(x +h/2.0, y +k2*h/2.0, z +l2*h/2.0);
    l3 = G(x +h/2.0, y +k2*h/2.0, z +l2*h/2.0);
    k4 = F(x +h, y +k3*h, z +l3*h);
    l4 = G(x +h, y +k3*h, z +l3*h);
    l4 = G(x, y +k3*h, z +l3*h);
    f = y + (k1 +2.0*k2 +2.0*k3 +k4)*h/6.0;
    g = z + (l1 +2.0*l2 +2.0*l3 +l4)*h/6.0;
//printf("x=%e f=%e %e %e %e %e x=%e\n", x, f, k1, k2, k3, k4, x);
//printf("x=%e g=%e %e %e %e %e z=%e", x, g, l1, l2, l3, l4, z);
//getchar();
    y = f;
    z = g;
}

// Write the differential equation in function of dy/dx: Called F
double F(double x, double y){
    double r = 0.0, l = k*x*(1.0-exp(-(x*ustar)/(RvD*nu)));
    r = (-nu +sqrt(nu*nu +4.0*l*l*ustar*ustar))/(2.0*l*l);
    return r;
}

// Runge-Kutta method for 4th order in 1st ODE
double RK4(double x, double y, double h){
    double f = 0.0;
    double k1 = 0.0, k2 = 0.0, k3 = 0.0, k4 = 0.0;
    k1 = F(x, y);
    k2 = F(x +h/2.0, y +k1*h/2.0);
    k3 = F(x +h/2.0, y +k2*h/2.0);
    k4 = F(x +h, y +k3*h);
    f = y + (k1 +2.0*k2 +2.0*k3 +k4)*h/6.0;
//printf("x=%e f=%e %e %e %e %e x=%e\n", x, f, k1, k2, k3, k4, x);
//getchar();
    return f;
}

int main(){
    double x = 0.0, y = 0.0, h = 0.001;
    FILE *RK = fopen("RK4.txt","w");
    fprintf(RK, "#z u(z) tau(z) dpress(z)/dz fxF(z)/phi(z) fyF(z)/phi(z)\n");
    fprintf(RK, "%e\t%e\t%e\t%e\t%e\t%e\n", x, y, ustar*ustar, 0.0, 0.0, 0.0);
    for (int i = 1; i < 5000000; i++){
        x = i*h;
        y = RK4(x, y, h);
        if (x/50 > 1){
            if (x/500 > 1){
                if (x/5000 > 1){
                    if (i%1000 == 0)
                        fprintf(RK, "%e\t%e\t%e\t%e\t%e\t%e\n", x*ustar/nu, y/ustar, 1.0, 0.0, 0.0, 0.0);
                } else {
                    if (i%100 == 0)
                        fprintf(RK, "%e\t%e\t%e\t%e\t%e\t%e\n", x*ustar/nu, y/ustar, 1.0, 0.0, 0.0, 0.0);
                }
            } else {
                if (i%10 == 0)
                    fprintf(RK, "%e\t%e\t%e\t%e\t%e\t%e\n", x*ustar/nu, y/ustar, 1.0, 0.0, 0.0, 0.0);
            }
        } else {
            fprintf(RK, "%e\t%e\t%e\t%e\t%e\t%e\n", x*ustar/nu, y/ustar, 1.0, 0.0, 0.0, 0.0);
        }
    }
    fclose(RK);
    return 0;
}
