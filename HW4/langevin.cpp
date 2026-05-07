#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "BM.h"
#include <iostream>
#include <math.h>

#define dim 3
#define dt 0.01
#define M (int)1e+3 //time step
#define N (int)1e+3 //number of ensemble

void eom(double *x, double *v,double T){
    int i;
    for(i=0;i<dim;i++){
        v[i] = v[i]-v[i]*dt+sqrt(2.*T*dt)*gaussian_rand();
        x[i] = x[i]+v[i]*dt;
    }
}

void init_array(double *x, double *v){
    int i;
    for(i=0;i<dim;i++){
        x[i]=0.;
        v[i]=0.;
    }
}

double msd(double *x, double *xOld){
    double r=0.;
    int i;
    for(i=0;i<dim;i++){
        r = r+(x[i]-xOld[i])*(x[i]-xOld[i]);
    }
    return r;
}

double corre(double *v, double *vOld){
    double r=0.;
    int i;
    for(i=0;i<dim;i++){
        r = r+v[i]*vOld[i];
    }
    return r;
}

int main(void){
    double msdArray[M], cArray[M];
    double x[dim], v[dim], Temp,xOld[dim],vOld[dim];
    for(Temp = 0.0001;Temp<2;Temp*=10.){
        char filename[128];
        sprintf(filename, "./output/result/temp_%.4f", Temp);
        //char filename1[128] = "./output/test_coor_vel.dat";
        std::ofstream file(filename);
        int i,j,k;

        for(i=0;i<M;i++){
            msdArray[i]=0.;
            cArray[i]=0.;
        }
        
        init_array(x,v); //For some reason, we cannot do ensemble of same initial value, so we initialise it once here.
        for(i=0;i<N;i++){//this loop is ensemble
            //init_array(x,v);
            bool check = true;
            for(j=0;j<M;j++){ //This loop is time
                //file1 << j*dt << "\t" << x[0] << " " << x[1] << " "<< x[2] << " "<< v[0] << " "<< v[1] << " "<< v[2] <<std::endl;
                eom(x,v,Temp);
                if(check){
                    for(k=0;k<dim;k++){
                        xOld[k]=x[k];
                        vOld[k]=v[k];
                    }
                    check = false;
                }
                msdArray[j]=msdArray[j]+msd(x,xOld);
                cArray[j]=cArray[j]+corre(v,vOld);
            }
        }

        for(j=0;j<M;j++){
            file << j*dt <<"\t"<< msdArray[j]/N <<"\t"<< cArray[j]/N<< std::endl;
        }
        //file1.close();
        file.close();

    }

    
}
