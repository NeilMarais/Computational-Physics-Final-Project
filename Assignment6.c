#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
double Bh, BJ;
int L;
double transitionprobabilities[18];
double antitransitionprobabilities[18];
int S;

//SECTION 1.1: Importing the ran1 function from Numerical Recipes
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

char *outputfilename1; 
FILE *outputfile1; 
char *outputfilename2; 
FILE *outputfile2; 
int S,N;
double H; //the energy of the system after a full sweep
double EaveragedivJ, E2averagedivJ2;
long r;
double percentrejected;
int L;
int BoundaryConds;
float ran1(long *idum)
/*“Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
successive deviates in a sequence. RNMX should approximate the largest floating value that is
less than 1.*/
{
int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;
if (*idum <= 0 || !iy) { //Initialize.
if (-(*idum) < 1) *idum=1; //Be sure to prevent idum = 0.
else *idum = -(*idum);
for (j=NTAB+7;j>=0;j--) { //Load the shuffle table (after 8 warm-ups).
k=(*idum)/IQ;
*idum=IA*(*idum-k*IQ)-IR*k;
if (*idum < 0) *idum += IM;
if (j < NTAB) iv[j] = *idum;
}
iy=iv[0];
}
k=(*idum)/IQ; //Start here when not initializing.
*idum=IA*(*idum-k*IQ)-IR*k; //Compute idum=(IA*idum) % IM without overif (*idum < 0) *idum += IM; flows by Schrage’s method.
j=iy/NDIV;// Will be in the range 0..NTAB-1.
iy=iv[j]; //Output previously stored value and refill the
iv[j] = *idum; //shuffle table.
if ((temp=AM*iy) > RNMX) return RNMX; //Because users don’t expect endpoint values.
else return temp;
}

void calculatetransitionprobabilities(){
int i = 0;
int j = 0;
for(i = 0; i<9; i++)
{transitionprobabilities[i] = exp(-2*BJ*(-4+i) -Bh);
//printf("Transition probability1 : %f\n", transitionprobabilities[i]);
}
for(j =0; j<9; j++)
{transitionprobabilities[j+9] = exp(2*BJ*(-4+j) -Bh);
//printf("Transition probability2 : %f\n", transitionprobabilities[j+9]);
}  
}

void calculateantiferromagneticprobabilities(){
int i = 0;
int j = 0;
for(i = 0; i<9; i++)
{antitransitionprobabilities[i] = exp(2*BJ*(-4+i) -Bh);
//printf("Transition probability1 : %f\n", transitionprobabilities[i]);
}
for(j =0; j<9; j++)
{antitransitionprobabilities[j+9] = exp(-2*BJ*(-4+j) -Bh);
//printf("Transition probability2 : %f\n", transitionprobabilities[j+9]);
}  
}

// Periodic Boundary Conditions
int periodic(int lattice[L][L],int S){// Simulate S sweeps for a lattice of size LxL
int numsweeps = 0;
H = 0;
double Mgenerated;
double Maverage, MAbsoluteAverage, M2Average;
double MTotal = 0;
Maverage = 0;
MAbsoluteAverage = 0;
M2Average = 0;
// Now we want to adjust the spins sequentially according to the acceptance probability
int x = 0;
int y = 0;
float acceptanceprobability;
int sumofneighbourspins =0;
for(numsweeps = 0; numsweeps <S;numsweeps++){
     //printf("First for loop ran\n");
    for (x=0; x<L; x=x+1){
   // printf("2nd for loop ran\n");
        for (y= 0; y<L; y++){
       // printf("3rd for loop ran\n");
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}

        //printf("Sum of neighbour spins first:%d\n",sumofneighbourspins);

        if (lattice[x][y] ==1) {
             //printf("First if loop ran\n");
            if(transitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = transitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = -1;}
        }
        else{
        if(transitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = transitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = 1;}}
        }
        }
    int magnetisation =0;
    x = 0;
    y = 0; 
    for (x = 0; x<L; x++){
        for(y=0; y<L; y++){
            magnetisation = magnetisation + lattice[x][y];}}

    //printf("Sum over neighbouring spins:%d\n",Hsumofneighbourspins); - potential issue here
    //printf("M:%d\n",magnetisation);
    if (numsweeps >= 5000){
        int Hsumofneighbourspins = 0;

            // Now to calculate H for the lattice after a full sweep
    for (x=0; x<L; x++){
    for (y= 0; y<L; y++){
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}
        Hsumofneighbourspins += sumofneighbourspins;}}
    Maverage =  Maverage + (1/(double)(S - 5000))*magnetisation;
    //printf("<M>:%d\n",Maverage);
    MAbsoluteAverage =  MAbsoluteAverage + (1/(double)(S - 5000))*abs(magnetisation);
   // printf("<|M|>:%d\n",Maverage);
    M2Average =  M2Average + (1/(double)(S - 5000))*magnetisation*magnetisation;
   // printf("<M2>:%d\n",M2Average);
    EaveragedivJ = EaveragedivJ + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation);
    E2averagedivJ2 = E2averagedivJ2 + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation)*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation) ;
    }   
    fprintf(outputfile1,"%d,",magnetisation);  
    fprintf(outputfile1,"%d\n",abs(magnetisation));
    }
  fprintf(outputfile1,"%f,",Maverage); 
  fprintf(outputfile1,"%f\n",MAbsoluteAverage); 
  fprintf(outputfile1,"%f,",M2Average); 
  fprintf(outputfile1,"%f\n",EaveragedivJ);  
  fprintf(outputfile1,"%f,0",E2averagedivJ2);  
return 0;
}
// END OF PERIODIC BOUNDARY CONDITIONS

// Polarised Boundary Conditions
int polarised(int lattice[L][L],int S){// Simulate S sweeps for a lattice of size LxL

// initialise the polarised boundary conditions
int rowcount = 0;
for (rowcount = 0; rowcount < L; rowcount++)
{lattice[rowcount][0] = 1; // left hand side of lattice is all spin up
 lattice[rowcount][L-1] = -1; // right hand side of lattice is all spin down
 }
int numsweeps = 0;
H = 0;
double Mgenerated;
double Maverage, MAbsoluteAverage, M2Average;
double MTotal = 0;
Maverage = 0;
MAbsoluteAverage = 0;
M2Average = 0;
// Now we want to adjust the spins sequentially according to the acceptance probability
int x = 0;
int y = 0;
float acceptanceprobability;
int sumofneighbourspins =0;
for(numsweeps = 0; numsweeps <S;numsweeps++){
     //printf("First for loop ran\n");
    for (x=0; x<L; x=x+1){
   // printf("2nd for loop ran\n");
        for (y= 1; y<L-1; y++){
       // printf("3rd for loop ran\n");
       // printf("3rd for loop ran\n");
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}

        //printf("Sum of neighbour spins first:%d\n",sumofneighbourspins);

        if (lattice[x][y] ==1) {
             //printf("First if loop ran\n");
            if(transitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = transitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = -1;}
        }
        else{
        if(transitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = transitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = 1;}}
        }
        }
    int magnetisation =0;
    x = 0;
    y = 0; 
    for (x = 0; x<L; x++){
        for(y=0; y<L; y++){
            magnetisation = magnetisation + lattice[x][y];}}

    //printf("Sum over neighbouring spins:%d\n",Hsumofneighbourspins); - potential issue here
    //printf("M:%d\n",magnetisation);
    if (numsweeps >= 5000){
        int Hsumofneighbourspins = 0;

            // Now to calculate H for the lattice after a full sweep
    for (x=0; x<L; x++){
    for (y= 0; y<L; y++){
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}
        Hsumofneighbourspins += sumofneighbourspins;}}
    Maverage =  Maverage + (1/(double)(S - 5000))*magnetisation;
    //printf("<M>:%d\n",Maverage);
    MAbsoluteAverage =  MAbsoluteAverage + (1/(double)(S - 5000))*abs(magnetisation);
   // printf("<|M|>:%d\n",Maverage);
    M2Average =  M2Average + (1/(double)(S - 5000))*magnetisation*magnetisation;
   // printf("<M2>:%d\n",M2Average);
    EaveragedivJ = EaveragedivJ + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation);
    E2averagedivJ2 = E2averagedivJ2 + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation)*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation) ;
    }   
    fprintf(outputfile1,"%d,",magnetisation);  
    fprintf(outputfile1,"%d\n",abs(magnetisation));
    }
  fprintf(outputfile1,"%f,",Maverage); 
  fprintf(outputfile1,"%f\n",MAbsoluteAverage); 
  fprintf(outputfile1,"%f,",M2Average); 
  fprintf(outputfile1,"%f\n",EaveragedivJ);  
  fprintf(outputfile1,"%f,0",E2averagedivJ2);

   int totalringspinplus = 0;
   int totalringspinminus = 0;
  int count = 0;
  for (count =0; count <L; count ++)
  {totalringspinplus = totalringspinplus + lattice[count][0];
  totalringspinminus = totalringspinminus + lattice[count][L-1];}
  printf("Total spin of left edge =%d\n",totalringspinplus);   
  printf("Total spin of right edge =%d\n",totalringspinminus);  
return 0;
}
// END OF POLARISED BOUNDARY CONDITIONS


// SEMI POLARISED BOUNDARY CONDITIONS
int semipolarised(int lattice[L][L],int S){
// initialise the semi-polarised boundary conditions
int rowcount = 0;
for (rowcount = 0; rowcount < L; rowcount++)
{
    if(BoundaryConds ==2)
    {lattice[rowcount][0] = 1;}// left hand side of lattice is all spin up
    else if (BoundaryConds ==3)
    {lattice[rowcount][0] = -1;} // left hand side of lattice is all spin down
}
int numsweeps = 0;
H = 0;
double Mgenerated;
double Maverage, MAbsoluteAverage, M2Average;
double MTotal = 0;
Maverage = 0;
MAbsoluteAverage = 0;
M2Average = 0;
// Now we want to adjust the spins sequentially according to the acceptance probability
int x = 0;
int y = 0;
float acceptanceprobability;
int sumofneighbourspins =0;
for(numsweeps = 0; numsweeps <S;numsweeps++){
     //printf("First for loop ran\n");
    for (x=0; x<L; x=x+1){
   // printf("2nd for loop ran\n");
        for (y= 1; y<L; y++){
       // printf("3rd for loop ran\n");
       // printf("3rd for loop ran\n");
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}

        //printf("Sum of neighbour spins first:%d\n",sumofneighbourspins);

        if (lattice[x][y] ==1) {
             //printf("First if loop ran\n");
            if(transitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = transitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = -1;}
        }
        else{
        if(transitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = transitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = 1;}}
        }
        }
    int magnetisation =0;
    x = 0;
    y = 0; 
    for (x = 0; x<L; x++){
        for(y=0; y<L; y++){
            magnetisation = magnetisation + lattice[x][y];}}

    //printf("Sum over neighbouring spins:%d\n",Hsumofneighbourspins); - potential issue here
    //printf("M:%d\n",magnetisation);
    if (numsweeps >= 5000){
        int Hsumofneighbourspins = 0;

            // Now to calculate H for the lattice after a full sweep
    for (x=0; x<L; x++){
    for (y= 0; y<L; y++){
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}
        Hsumofneighbourspins += sumofneighbourspins;}}
    Maverage =  Maverage + (1/(double)(S - 5000))*magnetisation;
    //printf("<M>:%d\n",Maverage);
    MAbsoluteAverage =  MAbsoluteAverage + (1/(double)(S - 5000))*abs(magnetisation);
   // printf("<|M|>:%d\n",Maverage);
    M2Average =  M2Average + (1/(double)(S - 5000))*magnetisation*magnetisation;
   // printf("<M2>:%d\n",M2Average);
    EaveragedivJ = EaveragedivJ + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation);
    E2averagedivJ2 = E2averagedivJ2 + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation)*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation) ;
    }   
    fprintf(outputfile1,"%d,",magnetisation);  
    fprintf(outputfile1,"%d\n",abs(magnetisation));
    }
  fprintf(outputfile1,"%f,",Maverage); 
  fprintf(outputfile1,"%f\n",MAbsoluteAverage); 
  fprintf(outputfile1,"%f,",M2Average); 
  fprintf(outputfile1,"%f\n",EaveragedivJ);  
  fprintf(outputfile1,"%f,0",E2averagedivJ2); 

  int totalringspin = 0;
  int count = 0;
  for (count =0; count <L; count ++)
  {totalringspin = totalringspin + lattice[count][0];}
  printf("Total spin of polarised edge =%d\n",totalringspin); 
return 0;}
// END OF SEMI POLARISED BOUNDARY CONDITIONS

// NON-INTERACTION BOUNDARY CONDITIONS
int noninteraction(int lattice[L][L],int S){
// initialise the non-interaction boundary conditions
int rowcount = 0;
int colcount =1;
for (rowcount = 0; rowcount < L; rowcount++)
{   if(BoundaryConds == 4)
    {lattice[rowcount][0] = 1;
    lattice[rowcount][L-1] =1;
     for (colcount = 1; colcount <L-1; colcount ++)
     {lattice[0][colcount] =1;
     lattice[L-1][colcount] = 1;}}
    else if (BoundaryConds == 5)
    {lattice[rowcount][0] = -1;
    lattice[rowcount][L-1] = -1;
     for (colcount = 1; colcount <L-1; colcount ++)
     {lattice[0][colcount] = -1;
     lattice[L-1][colcount] = -1;}} 
}
int numsweeps = 0;
H = 0;
double Mgenerated;
double Maverage, MAbsoluteAverage, M2Average;
double MTotal = 0;
Maverage = 0;
MAbsoluteAverage = 0;
M2Average = 0;
// Now we want to adjust the spins sequentially according to the acceptance probability
int x = 0;
int y = 0;
float acceptanceprobability;
int sumofneighbourspins =0;
for(numsweeps = 0; numsweeps <S;numsweeps++){
     //printf("First for loop ran\n");
    for (x=1; x<L-1; x=x+1){
   // printf("2nd for loop ran\n");
        for (y= 1; y<L-1; y++){
       // printf("3rd for loop ran\n");
       // printf("3rd for loop ran\n");
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}

        //printf("Sum of neighbour spins first:%d\n",sumofneighbourspins);

        if (lattice[x][y] ==1) {
             //printf("First if loop ran\n");
            if(transitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = transitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = -1;}
        }
        else{
        if(transitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = transitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = 1;}}
        }
        }
    int magnetisation =0;
    x = 0;
    y = 0; 
    for (x = 0; x<L; x++){
        for(y=0; y<L; y++){
            magnetisation = magnetisation + lattice[x][y];}}

    //printf("Sum over neighbouring spins:%d\n",Hsumofneighbourspins); - potential issue here
    //printf("M:%d\n",magnetisation);
    if (numsweeps >= 5000){
        int Hsumofneighbourspins = 0;

            // Now to calculate H for the lattice after a full sweep
    for (x=0; x<L; x++){
    for (y= 0; y<L; y++){
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}
        Hsumofneighbourspins += sumofneighbourspins;}}
    Maverage =  Maverage + (1/(double)(S - 5000))*magnetisation;
    //printf("<M>:%d\n",Maverage);
    MAbsoluteAverage =  MAbsoluteAverage + (1/(double)(S - 5000))*abs(magnetisation);
   // printf("<|M|>:%d\n",Maverage);
    M2Average =  M2Average + (1/(double)(S - 5000))*magnetisation*magnetisation;
   // printf("<M2>:%d\n",M2Average);
    EaveragedivJ = EaveragedivJ + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation);
    E2averagedivJ2 = E2averagedivJ2 + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation)*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation) ;
    }   
    fprintf(outputfile1,"%d,",magnetisation);  
    fprintf(outputfile1,"%d\n",abs(magnetisation));
    }
  fprintf(outputfile1,"%f,",Maverage); 
  fprintf(outputfile1,"%f\n",MAbsoluteAverage); 
  fprintf(outputfile1,"%f,",M2Average); 
  fprintf(outputfile1,"%f\n",EaveragedivJ);  
  fprintf(outputfile1,"%f,0",E2averagedivJ2);  
  int totalringspin =0;
  int index =  0;
  int innercount  = 0;
  for (index = 0; index<L; index++)
  {totalringspin = totalringspin+lattice[index][0];
   totalringspin = totalringspin + lattice[index][L-1];}
    for (innercount = 1; innercount <L-1; innercount++)
    {totalringspin = totalringspin + lattice[0][innercount] + lattice[L-1][innercount];}

printf("Spin of outer spins = %d\n",totalringspin);
return 0;

}
// END OF NON-INTERACTION BOUNDARY CONDITIONS
int ferromagnetic(int lattice[L][L],int S){
int numsweeps = 0;
H = 0;
double Mgenerated;
double Maverage, MAbsoluteAverage, M2Average;
double MTotal = 0;
Maverage = 0;
MAbsoluteAverage = 0;
M2Average = 0;
// Now we want to adjust the spins sequentially according to the acceptance probability
int x = 0;
int y = 0;
float acceptanceprobability;
int sumofneighbourspins =0;
for(numsweeps = 0; numsweeps <S;numsweeps++){
     //printf("First for loop ran\n");
    for (x=0; x<L; x=x+1){
   // printf("2nd for loop ran\n");
        for (y= 0; y<L; y++){
       // printf("3rd for loop ran\n");
        int sumofneighbourspins = 0;
        int edge = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];
        edge =1;}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];
         edge =1;}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];
         edge =1;}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];
         edge =1;}

        //printf("Sum of neighbour spins first:%d\n",sumofneighbourspins);

        if (lattice[x][y] ==1) {
             //printf("First if loop ran\n");
             if (edge == 0){
            if(antitransitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = antitransitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}}
        else if (edge == 1)
        {if(transitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = transitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = -1;}
        }
        else{
            if (edge ==0){
        if(antitransitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = antitransitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}}
            if (edge == 1){if(transitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = transitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = 1;}}
        }
        }
    int magnetisation =0;
    x = 0;
    y = 0; 
    for (x = 0; x<L; x++){
        for(y=0; y<L; y++){
            magnetisation = magnetisation + lattice[x][y];}}

    //printf("Sum over neighbouring spins:%d\n",Hsumofneighbourspins); - potential issue here
    //printf("M:%d\n",magnetisation);
    if (numsweeps >= 5000){
        int Hsumofneighbourspins = 0;

            // Now to calculate H for the lattice after a full sweep
    for (x=0; x<L; x++){
    for (y= 0; y<L; y++){
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}
        Hsumofneighbourspins += sumofneighbourspins;}}
    Maverage =  Maverage + (1/(double)(S - 5000))*magnetisation;
    //printf("<M>:%d\n",Maverage);
    MAbsoluteAverage =  MAbsoluteAverage + (1/(double)(S - 5000))*abs(magnetisation);
   // printf("<|M|>:%d\n",Maverage);
    M2Average =  M2Average + (1/(double)(S - 5000))*magnetisation*magnetisation;
   // printf("<M2>:%d\n",M2Average);
    EaveragedivJ = EaveragedivJ + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation);
    E2averagedivJ2 = E2averagedivJ2 + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation)*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation) ;
    }   
    fprintf(outputfile1,"%d,",magnetisation);  
    fprintf(outputfile1,"%d\n",abs(magnetisation));
    }
  fprintf(outputfile1,"%f,",Maverage); 
  fprintf(outputfile1,"%f\n",MAbsoluteAverage); 
  fprintf(outputfile1,"%f,",M2Average); 
  fprintf(outputfile1,"%f\n",EaveragedivJ);  
  fprintf(outputfile1,"%f,0",E2averagedivJ2);    
    return 0;}
int antiferromagnetic(int lattice[L][L],int S){
    int numsweeps = 0;
H = 0;
double Mgenerated;
double Maverage, MAbsoluteAverage, M2Average;
double MTotal = 0;
Maverage = 0;
MAbsoluteAverage = 0;
M2Average = 0;
// Now we want to adjust the spins sequentially according to the acceptance probability
int x = 0;
int y = 0;
float acceptanceprobability;
int sumofneighbourspins =0;
for(numsweeps = 0; numsweeps <S;numsweeps++){
     //printf("First for loop ran\n");
    for (x=0; x<L; x=x+1){
   // printf("2nd for loop ran\n");
        for (y= 0; y<L; y++){
       // printf("3rd for loop ran\n");
        int sumofneighbourspins = 0;
        int edge = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];
        edge =1;}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];
         edge =1;}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];
         edge =1;}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];
         edge =1;}

        //printf("Sum of neighbour spins first:%d\n",sumofneighbourspins);

        if (lattice[x][y] ==1) {
             //printf("First if loop ran\n");
             if (edge == 0){
            if(transitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = transitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}}
        else if (edge == 1)
        {if(antitransitionprobabilities[sumofneighbourspins + 4] < 1)
        {acceptanceprobability = antitransitionprobabilities[sumofneighbourspins + 4];}
        else {acceptanceprobability =1;}}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = -1;}
        }
        else{
            if (edge ==0){
        if(transitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = transitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}}
            if (edge == 1){if(antitransitionprobabilities[sumofneighbourspins + 13] < 1){
            acceptanceprobability = antitransitionprobabilities[sumofneighbourspins + 13];}
        else {acceptanceprobability =1;}}
        //printf("Acceptance prob:%d\n",acceptanceprobability);
        float choosingnum = ran1(&r);
        if (choosingnum < acceptanceprobability)
        {lattice[x][y] = 1;}}
        }
        }
    int magnetisation =0;
    x = 0;
    y = 0; 
    for (x = 0; x<L; x++){
        for(y=0; y<L; y++){
            magnetisation = magnetisation + lattice[x][y];}}

    //printf("Sum over neighbouring spins:%d\n",Hsumofneighbourspins); - potential issue here
    //printf("M:%d\n",magnetisation);
    if (numsweeps >= 5000){
        int Hsumofneighbourspins = 0;

            // Now to calculate H for the lattice after a full sweep
    for (x=0; x<L; x++){
    for (y= 0; y<L; y++){
        int sumofneighbourspins = 0;
        if (x == 0){sumofneighbourspins += lattice[L-1][y];}
        else{sumofneighbourspins += lattice[x-1][y];}

        if (y == 0){sumofneighbourspins += lattice[x][L-1];}
        else{sumofneighbourspins += lattice[x][y-1];}

        if (x == L-1){sumofneighbourspins += lattice[0][y];}
        else{sumofneighbourspins += lattice[x+1][y];}

        if (y == L -1){sumofneighbourspins += lattice[x][0];}
        else{sumofneighbourspins += lattice[x][y+1];}
        Hsumofneighbourspins += sumofneighbourspins;}}
    Maverage =  Maverage + (1/(double)(S - 5000))*magnetisation;
    //printf("<M>:%d\n",Maverage);
    MAbsoluteAverage =  MAbsoluteAverage + (1/(double)(S - 5000))*abs(magnetisation);
   // printf("<|M|>:%d\n",Maverage);
    M2Average =  M2Average + (1/(double)(S - 5000))*magnetisation*magnetisation;
   // printf("<M2>:%d\n",M2Average);
    EaveragedivJ = EaveragedivJ + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation);
    E2averagedivJ2 = E2averagedivJ2 + (1/(double)(S - 5000))*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation)*(-Hsumofneighbourspins + (Bh/BJ)*magnetisation) ;
    }   
    fprintf(outputfile1,"%d,",magnetisation);  
    fprintf(outputfile1,"%d\n",abs(magnetisation));
    }
  fprintf(outputfile1,"%f,",Maverage); 
  fprintf(outputfile1,"%f\n",MAbsoluteAverage); 
  fprintf(outputfile1,"%f,",M2Average); 
  fprintf(outputfile1,"%f\n",EaveragedivJ);  
  fprintf(outputfile1,"%f,0",E2averagedivJ2); 
    return 0;}

int main(int argc, char *argv[]){
    if (argc<7) {
		printf("Usage: Assignment6 Bh Bj\n");
		printf("Bh= product of 1/kT and external field strength h - BJ = product of 1/kT and spin-spin interaction strength J - L = size of lattice - S = number of sweeps -filename1 = name of file where averages output will be written - B = number which determines the boundary condition\n");
		return 1;
	}
Bh = atof(argv[1]); 
BJ = atof(argv[2]);
L = atoi(argv[3]); 
S = atoi(argv[4]);
outputfilename1=argv[5];
BoundaryConds=atoi(argv[6]);	
outputfile1=fopen(outputfilename1,"w"); // Try to open the output file for writing. Return a pointer to a file object.
//outputfile2=fopen(outputfilename2,"w"); // Try to open the output file for writing. Return a pointer to a file object.

// We will get outputfile=0 if fopen failed. In C, any non-zero integer is seen as True in an if statement, while 0 is seen as False. So !0=True.
if (!outputfile1) {
    printf("Cannot create %s.\n",outputfilename1);
	return 1;}

// We will get outputfile=0 if fopen failed. In C, any non-zero integer is seen as True in an if statement, while 0 is seen as False. So !0=True.
//if (!outputfile2) {
    //printf("Cannot create %s.\n",outputfilename2);
	//return 1;}
// Initialise the PRNG using the current time as a seed
r=(long) time(NULL);
r=-r;
ran1(&r);
printf("Bh = %f\n", Bh);
printf("BJ = %f\n", BJ);
printf("L = %d\n", L);
printf("S = %d\n", S);
calculatetransitionprobabilities();
calculateantiferromagneticprobabilities();
	printf("Generating S=%d samples on lattice with L=%d.\n",S,L);
	printf("Will write data to %s.\n",outputfilename1);
		
	
   
	// Set the averages to 0
	EaveragedivJ = 0;
    E2averagedivJ2 = 0;
        // Create a 2d array to represent the lattice - initialise the lattice with all the spins in the spin up state
        // The spin up state is represented by a value of 1 and the spin down state is represented by a value of -1
    int Nup = L*L; // number of spins in the spin up state
    int lattice[L][L];
    int i=0; int j=0;
    for (int i =0; i<L; i++){
        for (int j = 0; j<L;j++){
            lattice[i][j] =1;
        }

    }

    if (BoundaryConds == 0)
	{periodic(lattice,S);}
	else if (BoundaryConds == 1)
	{polarised(lattice,S);}
    else if (BoundaryConds == 2)
	{semipolarised(lattice,S);}
    else if (BoundaryConds == 3)
	{semipolarised(lattice,S);}
    else if (BoundaryConds == 4)
	{noninteraction(lattice,S);}
    else if (BoundaryConds == 5)
	{noninteraction(lattice,S);}
    else if (BoundaryConds == 6)
    {ferromagnetic(lattice,S);}
    else
    {antiferromagnetic(lattice,S);}
    // Force the OS to write any buffered data to file. This is important if we plan to read from the file directly after the program ends. Maybe fclose does this as well, but it is better to be sure...
    fflush(outputfile1);
    // Close the file
	fclose(outputfile1);
	 // Force the OS to write any buffered data to file. This is important if we plan to read from the file directly after the program ends. Maybe fclose does this as well, but it is better to be sure...
   // fflush(outputfile2);
    // Close the file
	//fclose(outputfile2);



	printf("Done!\n");
    printf("Test done!\n");
	return 0;
}