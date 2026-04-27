#include <stdio.h> // for printf, etc                                                                                                                                       
#include <stdlib.h> // for rand(), etc                                                                                                                                      
#include <math.h> // for sin(),cos(), etc                                                                                                                                   
#include <iostream>// for cout, etc 
#include <iomanip>// for setprecision()
#include <fstream> // for ifstream/ofstream                                                                                                                                 
#include <time.h>// for time(NULL), etc                                                                                                                                     

int main(void){
  int i, count = 0, max = 1e+9;
  int step = 1e+2;
  double x,y,z,pi;
  char fname[128];
  char hwname[128];
  std::ofstream file;
  srand(time(NULL));  // "time(NULL)" as a seed of ramdom number
  //  sprintnf(fname,"coord%d.dat",max); // Define the file name for fname[128]
  snprintf(fname, sizeof(fname), "./output/coord%d.dat", max);
  

  
  /*
  file.open(fname); // "file" with the name of fname[128]
  for(i=0;i<max;i++){
    x = (double)rand()/RAND_MAX;
    y = (double)rand()/RAND_MAX;
    z = x*x + y*y;
    if(z<1){
      count++;
      file << x <<" "<<y <<std::endl;
    }
  }
  
  file.close(); //close fname file
  pi = (double)count / max * 4.;
  */
  
  //Add the result to hwname
  bool check = false;
  snprintf(hwname, sizeof(hwname), "output/coord%d_hw_rand.dat", max);
  file.open(hwname);
  count = 0;
  for(i=0;i<max;i++){
    x = (double)rand()/RAND_MAX;
    y = (double)rand()/RAND_MAX;
    z = x*x + y*y;
    if(z<1){
      count++;
    }
    check = false;
    
    if(1<=i && i <1e+2){
      if(i%10==0){
        check = true;
      }
    } else if(1e+2<=i && i <1e+3){
      if(i%100 == 0){
        check = true;
      } 
    } else if(1e+3<=i && i <1e+4){
      if(i%1000 == 0){
        check = true;
      } 
    } else if(1e+4<=i && i <1e+5){
      if(i%10000 == 0){
        check = true;
      } 
    } else if(1e+5<=i && i <1e+6){
      if(i%100000 == 0){
        check = true;
      }
    } else if(1e+6<=i && i <1e+7){
      if(i%1000000 == 0){
        check = true;
      }
    } else if(1e+7<=i && i <1e+8){
      if(i%1000000 == 0){
        check = true;
      } 
    } else if(1e+8<=i && i <1e+9){
      if(i%100000000 == 0){
        check = true;
      } 
    } else if(1e+9<=i && i <1e+10){
      if(i%(int)1000000000 == 0){
        check = true;
      }
    }
    if(check){
      pi = (double)count / i *4;
      file << i << " "  << std::setprecision(21)<< pi << std::endl;
      //printf("here");
      //count << i;
    }
  }
  file.close();
  std::cout << "finished" << std::endl;
  //printf("%.20f\n",pi); //by C, %.20f -- Displaying with 20 decimal precision
  //std::cout<< std::setprecision(21) <<  pi  << std::endl; // by C++
  return 0;
}
