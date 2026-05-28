#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include "./BM.h"
#include <queue>



const int N = 1024; //number of disk
const int dim =   2;
const double L = 40., aCut = 3.0, dt = 0.01, tmax = 100.;



void square_lattice(double *x){
    double a = 1.1;
    double a1[2]={a,0.0},a2[2]={0.0,a};
    int i,j,d, n = 1;
    for(i=0;i<dim;i++){
        x[i]=2.1;
    }
    i = 1;
    j = 0;
    while(n<N){
        if(abs(x[int(n-1)*dim]+a1[0]+a2[0]) > L || abs(x[(int)(n-1)*dim+1]+a1[1]+a2[1]) > L){
            j++;
            i=0;
        }
        for(d=0;d<dim;d++){
            x[n*dim+d]=x[d]+(double)i*a1[d]+(double)j*a2[d];
        }
        
        i++;
        n++;
    }
}


void square_lattice_adaptive(double *x){
    int num_x = (int)sqrt(N)+1;
    int num_y = (int)sqrt(N)+1;
    int i,j,k=0;
    for(j=0;j<num_y;j++){
        for(i=0;i<num_x;i++){
            x[k*dim] = i*L/(double)num_x;
            x[k*dim+1] = j*L/(double)num_y;
            k++;
            if(k==N){
                break;
            }
        }
        if(k==N){
            break;
        }
    }
}

void hexagonal_lattice(double* x){
    double a = 1.5;
    double a1[2]={a,0.0},a2[2]={a/2.0,a*sqrt(3.0)/2};
    int i,j,d, n = 1;
    for(i=0;i<dim;i++){
        x[i]=2.1;
    }
    i = 1;
    j = 0;
    while(n<N){
        if(abs(x[(int)(n-1)*dim]+a1[0]+a2[0]) > L || abs(x[(int)(n-1)*dim+1]+a1[1]+a2[1]) > L){
            j++;
            i= -(int)j/2;
        }
        for(d=0;d<dim;d++){
            x[n*dim+d]=x[d]+(double)i*a1[d]+(double)j*a2[d];
        }
        
        i++;
        n++;
    }
}

void particle_size(double *a){ //set diameter of particle
    int i;
    double aVar = 0.0;
    for(i=0;i<N;i++){
        a[i]=1.0+aVar*gaussian_rand();
    }
}

void ini_array(double *x, int n){
    int i;
    for(i=0;i<n;i++){
        x[i]=0;
    }
}

bool interaction(double *x, double *a, double *f){
    double dx,dy,dr2, aij, w2, w6, w12,F;
    double dxt, dyt;
    int i, j;
    bool found = false;
    ini_array(f,(int)N*dim);
    for(i=0;i<N;i++){
        if(found) break;
        for(j=0;j<N;j++){
            if(i<j){
                dx = x[i*dim]-x[j*dim];
                dy = x[i*dim+1]-x[j*dim+1];
                dx  -= L*floor((dx+0.5*L)/L);
                dy  -= L*floor((dy+0.5*L)/L);
                dr2 = dx*dx+dy*dy;
                
                if(dr2==0){
                    std::cout <<"=============================\n";
                    std::cout << "particle i (x,y) = " << x[(int)i*dim] << ", " << x[(int)(i*dim+1)] << std::endl;
                    std::cout << "particle j (x,y) = " << x[j*dim] << ", " << x[j*dim+1] << std::endl;
                    std::cout <<"============Before Periodic Boundary=================\n";
                    std::cout << "dx = " << dxt << " dy= " << dyt << std::endl;
                    std::cout << "dx = " << dx << " dy = " << dy << std::endl;
                    std::cout<<"dr2 = 0 " << dr2 <<  std::endl;
                    std::cout << "i = " << i << " j = " << j << std::endl;
                    //found = true;
                    break;
                }
                
                if(dr2<aCut*aCut){
                    aij = 0.5*(a[i]+a[j]);
                    w2 = aij*aij/dr2;
                    w6 = w2*w2*w2;
                    w12 = w6*w6;
                    F = -12.*w12/dr2;
                    f[i*dim] -=F*dx;
                    f[i*dim+1] -= F*dy;
                    f[j*dim] += F*dx;
                    f[j*dim +1 ] += F*dy;
                }
            }
        }
    }
    return found;
}

void boundary(double *x){
    for(int i =0; i<(int)N*dim;i++){
        x[i] -= L*floor(x[i]/L);
    }
}


void eom_langevin(double *x, double *v, double *f, double T){
    double zeta=1.0;
    double fluc = sqrt(2.*zeta*T*dt);
    for(int i = 0; i<N; i++){
        v[i*dim] = v[i*dim] - zeta*dt*v[i*dim] + dt*f[i*dim] + fluc*gaussian_rand();
        x[i*dim] = x[i*dim] + dt*v[i*dim];
        v[i*dim+1] = v[i*dim+1] - zeta*dt*v[i*dim+1] + dt*f[i*dim+1] + fluc*gaussian_rand();
        x[i*dim+1] = x[i*dim+1] + dt*v[i*dim+1];   
    }
    boundary(x);
}

void eom_md(double *x, double *v, double *f, double *a){
    interaction(x,a,f);
    for(int n = 0; n < N; n++){
        for(int i =0; i<dim ; i++){
            x[n*dim+i]=x[n*dim+i]+v[n*dim+i]*dt+0.5*dt*dt*f[n*dim+i];
            v[n*dim+i]=v[n*dim+i]+0.5*dt*f[n*dim+i];
        }
    }
    interaction(x,a,f);
    for(int n =0; n<N ;n++){
        for(int i =0; i<dim ;i++){
            v[n*dim+i]=v[n*dim+i]+0.5*dt*f[n*dim+i];
        }
    }
}

double potential(double *x, double* a){
    double U = 0.0;
    double dx, dy, dr2, aij, w2,w6,w12;
    for(int i = 0; i<N ;i++){
        for(int j = 0; j<N ;j++){
            if(i<j){
                dx = x[i*dim]-x[j*dim];
                dy = x[i*dim+1]-x[j*dim+1];
                dx  -= L*floor((dx+0.5*L)/L);
                dy  -= L*floor((dy+0.5*L)/L);
                dr2 = dx*dx+dy*dy;
                if(dr2<aCut*aCut){
                    aij = 0.5*(a[i]+a[j]);
                    w2 = aij*aij/dr2;
                    w6 = w2*w2*w2;
                    w12 = w6*w6;
                    U = U +w12;
                }
            }
        }
    }
    return U;
}

double kinetic(double *v){
    double K =0.0;
    for(int n =0; n<N;n++){
        for(int i =0; i<dim; i++){
            K = K + 0.5*v[n*dim+i]*v[n*dim+i];
        }
    }
    return K;
}



void v_rms(double *v, double *vrms_array){
    ini_array(vrms_array, N);
    double vTemp;
    for(int i =0 ; i< N ; i++){
       vTemp = 0.0;
       for(int j = 0 ; j< dim ; j++){
            vTemp = vTemp + v[i*dim+j]*v[i*dim+j];
       }
       vrms_array[i] = sqrt(vTemp);
    }
}

void outputV(double *vrms_array,double t){
    char filename[512];
    sprintf(filename, "./output/MBDist/vDist_%05.2f.dat", t);
    std::ofstream file(filename);
    for(int n =0; n<N; n++){
        file << vrms_array[n] << std::endl;
    }
}

void outputVNew(double *vrms_array,double T){
    char filename[512];
    sprintf(filename, "./output/vDist_%05.2f.dat", T);
    std::ofstream file(filename);
    for(int n =0; n<N; n++){
        file << vrms_array[n] << std::endl;
    }
}

void outputEnergy(double t, double *x, double *v, double *a){
    char filename1[512];
    sprintf(filename1, "./output/Energy.dat");
    std::ofstream file1(filename1, std::ios::app);
    int i,j,k,n=0;
    i = 0;
    double U = potential(x,a);
    double K = kinetic(v);
    file1 << t << "\t" << U << "\t" << K << "\t" << U+K <<std::endl;
    file1.close();
}

void outputOne(double *x, double *v, double *a, double t, double T, const char *message){
    char filename1[512];
    sprintf(filename1, "./output/%s/xva_t%05.2fT%.2f.dat", message, t, T);
    std::ofstream file1(filename1);
    int i,j,k,n=0;
    while(n<N){
        file1 << n;
        for(i = 0; i<dim ; i++){
            file1 << "\t" << x[n*dim+i];
        }
        file1 << "\t" << a[n] <<std::endl;
        
        n++;
    }
    file1.close();
    
}

void outputDebug(double *x, double *v, double *a, const char *message){
    char filename1[512];
    sprintf(filename1, "./debug/%s.dat", message);
    std::ofstream file1(filename1);
    int i,j,k,n=0;
    double speedArray[N];
    v_rms(v,speedArray);
    while(n<N){
        file1 << n;
        for(i = 0; i<dim ; i++){
            file1 << "\t" << x[n*dim+i];
        }
        file1 << "\t" << a[n] << "\t" << speedArray[n] << std::endl;
        
        n++;
    }
    file1.close();
    
}

int main(){
    //double x1[N][dim];
    int i=0;
    double T0 = 5.0, tout=0.0;
    double *x1 = new double[N*dim];
    double *v1 = new double[N*dim];
    double* a = new double[N];
    double* f = new double[N*dim];
    double* vrmsArray = new double[N];
    bool check = false;
    double temp;
    double tempArray[4] = {0.1,0.2,0.5,1.0};
    std::queue <double> E;
    for(int k = 0; k < 4 ; k++){
        i = 0;
        temp = tempArray[k];
        ini_array(v1,(int)N*dim);
        particle_size(a);
        square_lattice(x1);
        //outputDebug(x1,v1,a,"ini_coor");
        while(dt*i<10){ //preparation of high temp equilibrium
            check = interaction(x1,a,f);
            eom_langevin(x1,v1,f,5.0);
            i++;
        }
        //outputDebug(x1,v1,a, "high_temp");
        i = 0;
        while (dt*i<10){ //preparation at T =0.9
            check = interaction(x1,a,f);
            eom_langevin(x1,v1,f,0.9);
            i++;
        }
        //outputDebug(x1,v1,a,"T = 0.9");
        i=0;
        while (dt*i<50){ //preparation at T =0.9
            check = interaction(x1,a,f);
            eom_langevin(x1,v1,f,0.9);
            i++;
        }
        //outputDebug(x1,v1,a,"T = 0.9 long");
        v_rms(v1,vrmsArray);
        outputVNew(vrmsArray, 0.9);
        i=0;
        while (dt*i<51){
            eom_md(x1,v1,f,a);
            
            if(dt*i>=tout){
                tout = tout + 5.;
                outputEnergy(dt*i,x1,v1,a);
                v_rms(v1,vrmsArray);
                outputV(vrmsArray,dt*i);
            }
            i++;
        }
    }
    
    

    delete[] x1;
    delete[] v1;
    delete[] a;
    delete[] f;
    return 0;
}
