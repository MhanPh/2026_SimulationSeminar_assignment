#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <iostream>


void eom(double tau){
    char filename[128];
    sprintf(filename, "./output/v_x_%.4f.dat", tau);
    std::ofstream file(filename);
    int i, N = 2000;
    double dt = 0.01, tout = 0.01;
    double xold = 1.0, vold = 0.0, xnew, vnew;
    for(i=1;i<N;i++){
        vnew = vold - tau*vold*dt - xold *dt;
        xnew = xold + vnew*dt;
        xold = xnew;
        vold = vnew;
        if(i*dt>tout){
            file << i*dt << "\t" << xold << "\t" << vold << std::endl;
            //tout = tout*1.5;
        }
    }
    file.close();
}

int main(void){
    double tau = 0.5;
    for(tau;tau<10;tau=tau*4){
        eom(tau);
    }
    return 1;
}
