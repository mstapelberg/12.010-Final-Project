#include <iostream>
#include "stdio.h"
#include "TH1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TF1.h"

#define nSamples 10
#define nCollisions 10
#define nParticles 1

int MCSimulator(){
/*Here we will include the Random samples*/

    const int N = 10;
    //double theta[N];
    TMatrixD *theta = new TMatrixD[nSamples];
    double phi[nSamples];
    double r[nSamples];
    for (int i = 0; i <nSamples; i++){
        theta->RndmArray(nSamples, theta);
    }
        theta.print();
    //for (int i = 0; i = nSamples; i++)
        //theta[i] = RndmArray(nSamples, *theta);




return 0;
}
