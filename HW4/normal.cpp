#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "BM.h"
#include <iostream>
#include <math.h>

int main(void){
    char filename[128];
    int m = 5000;
    sprintf(filename, "./output/normal_%d.dat", m);
    std::ofstream file;
    file.open(filename);
    int i;
    for(i=0;i<=m;i++){
        file << gaussian_rand() << std::endl;
    }
    file.close();
    return 0;
}
