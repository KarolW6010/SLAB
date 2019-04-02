#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "music_func.h"
#include "Eigen/Dense"

//Global Variables
/*
const int ANTS = 4;					//Number of antennas
const int GRID_RES = 200;			//Resolution of Music Spectra		
const int NUM_SAMPLES = 1024;		//Samples per autocorrelation matrix
const int FRAC_BITS = 43;			//Fractional bits for auto corr numbers
const int TAGS = 1;					//Number of tags
const float FREQ = 915000000;		//Tag freq
const float DOL = 0.5;				//Distance over lambda
const float LAMBDA = 300000000/FREQ;	//Wavelength
const float D = DOL*LAMBDA;				//Distance
*/

//Function Prototypes

int main(){
	setConst(200,1);
	float freq = 915*pow(10,6);			//Frequency of RFID Tags
	float lambda = 3*pow(10,8)/freq;	//Wavelength

	float array1[ANTS][3];		//ANTS antennas, (x,y,z) coordinates
	float array2[ANTS][3];		//^ but for array 2
	
	//Array 1
	array1[0][0] = 0.0*D;
	array1[0][1] = 0.0*D;
	array1[0][2] = 0.0;

	array1[1][0] = 1.0*D;
	array1[1][1] = 0.0*D;
	array1[1][2] = 0.0;

	array1[2][0] = 0.0*D;
	array1[2][1] = 1.0*D;
	array1[2][2] = 0.0;

	array1[3][0] = 1.0*D;
	array1[3][1] = 1.0*D;
	array1[3][2] = 0.0;
	
	//Array 2
	array2[0][0] = 2.0;
	array2[0][1] = 0.0;
	array2[0][2] = 0.0;

	array2[1][0] = 2.0+1.0*D;
	array2[1][1] = 0.0;
	array2[1][2] = 0.0;

	array2[2][0] = 2.0;
	array2[2][1] = 1.0*D;
	array2[2][2] = 0.0;

	array2[3][0] = 2.0+1.0*D;
	array2[3][1] = 1.0*D;
	array2[3][2] = 0.0;

/*
	//Scale appropriately
	for(int i=0; i<ANTS; i++){
		for(int j=0; j<3; j++){
		//	array1[i][j] = array1[i][j]*D;
		//	array2[i][j] = array2[i][j]*D;
			cout << "Array 1: " << array1[i][j] << ", Array 2: " << array2[i][j] << endl;
		}
		cout << endl;
	}
*/
	float *a1, *a2;
	a1 = new float[ANTS*3];
	a2 = new float[ANTS*3];
	a1 = &array1[0][0];
	a2 = &array2[0][0];

	//Find centers of array 1 and 2
	float *center1, *center2;
	center1 = new float[3];
	center2 = new float[3];

	float _Complex *R1, *R2;
	R1 = new float _Complex[ANTS*ANTS];
	R2 = new float _Complex[ANTS*ANTS];

	//Autocorrelation for array 1
	float valsR1[] = {.3650, .0364, .2522, .0207,.3533,.0427, .2515,.3450, .0305,.3421};
	float valsC1[] = {.0000,-.2565,-.0124,-.2541,.0000,.2467,-.0100,.0000,-.2384,.0000};
	
	//Autocorrelation for array 2
	float valsR2[] = {.3499,-.1586, .2511,-.1612, .3476,-.1562, .2573, .3524,-.1621, .3621};
	float valsC2[] = {.0000,-.1904,-.0009,-.1932, .0000, .1940,-.0066, .0000,-.1969, .0000};

	vec2autocorr(valsR1,valsC1,R1);
	vec2autocorr(valsR2,valsC2,R2);

	MatrixXcf *eigvecs1, *eigvecs2;
	MatrixXf *eigvals1, *eigvals2;
	eigvecs1 = new MatrixXcf;
	eigvecs2 = new MatrixXcf;
	eigvals1 = new MatrixXf;
	eigvals2 = new MatrixXf;

	autocorr2eig(R1, eigvecs1, eigvals1);
	autocorr2eig(R2, eigvecs2, eigvals2);
	
	MatrixXcf *subspace1, *subspace2;
	subspace1 = new MatrixXcf;
	subspace2 = new MatrixXcf;

	subspaceMat(eigvals1, eigvecs1, subspace1);
	subspaceMat(eigvals2, eigvecs2, subspace2);

	MatrixXf *S_music1, *thetas1, *phis1;
	S_music1 = new MatrixXf;		//Music Spectrum values
	thetas1 = new MatrixXf;		//Theta grid
	phis1 = new MatrixXf;		//Phi grid
	
	MatrixXf *S_music2, *thetas2, *phis2;
	S_music2 = new MatrixXf;		//Music Spectrum values
	thetas2 = new MatrixXf;		//Theta grid
	phis2 = new MatrixXf;		//Phi grid

	musicSpectrum(subspace1, a1, S_music1, thetas1, phis1);

	musicSpectrum(subspace2, a2, S_music2, thetas2, phis2);
	
	float *thLocs1, *phLocs1;
	thLocs1 = new float[TAGS];	//Theta values corresponding to peaks
	phLocs1 = new float[TAGS]; 	//Phi values corresponding to peaks 
	findPeaks(S_music1, thetas1, phis1, thLocs1, phLocs1);

	float *thLocs2, *phLocs2;
	thLocs2 = new float[TAGS];	//Theta values corresponding to peaks
	phLocs2 = new float[TAGS]; 	//Phi values corresponding to peaks 
	findPeaks(S_music2, thetas2, phis2, thLocs2, phLocs2);

	float *dist;
	dist = new float;
	Vector3f *midpoint;
	midpoint = new Vector3f;

//	cout << "\nBefore Func: Th1 = " << thLocs1[0] << ", Ph1 = " << phLocs1[0]
//		 << "\nBefore Func: Th2 = " << thLocs2[0] << ", Ph2 = " << phLocs2[0] << endl;

//	ang2loc(a1,*thLocs1,*phLocs1,a2,*thLocs2,*phLocs2,ANTS,dist,midpoint);

//	cout << "Midpoint:\n" << *midpoint << endl;
//	cout << "Distance = " << *dist << endl;

	Vector3f *locations;
	locations = new Vector3f[TAGS];
	bestLocal(a1,thLocs1,phLocs1,a2,thLocs2,phLocs2,locations);

	cout << "\n\nLocations:\n" << *locations << endl;

	return 0;
}













