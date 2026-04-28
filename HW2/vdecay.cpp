#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <cfloat>
#include <iomanip>

double v(double dt){
    char filename[128];
    std::ofstream file;
    sprintf(filename, "./output/v%.6f.dat", dt);
    file.open(filename);
    int i, nstep, mstep;
    double v=10., tmax = 1.e4, ttemp = 5.0,vout;
    nstep = tmax/(1e4*dt);
    mstep = (int)ttemp/(10*dt);
    for(i=0;i*dt<tmax;i++){
        if(v<DBL_EPSILON){
            std::cout << "out of precision" << std::endl;
            break;
        }
        v = v-v*dt;
        if(i*dt < (double)5.0){
            if(i % mstep ==0){
            file << std::setprecision(7)<< (double)i*dt << " " << std::setprecision(17)<< v << std::endl;
            }
        } else if(i*dt >= (double)5.0){
            if(i%nstep==0){
            file << std::setprecision(7)<< (double)i*dt << " " << std::setprecision(17)<< v << std::endl;
            }
        }
        if(i*dt == double(30)){
            vout = v;
        }
    }
    std::cout << "finished" << std::endl;
    file.close();
    return vout;
}


int main(void){
    double dt;
    char filename2[128]="./output/dt/dt.dat";
    std::ofstream file;
    file.open(filename2);
    for(dt=1.e-1; dt>=1.e-6; dt=dt/10.){
        file << dt << " " << std::setprecision(17)<< v(dt) << std::endl;
        //v(dt);
    }
    file.close();
    return 0;
}