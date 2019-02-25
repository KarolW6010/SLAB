#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "music_func.h"
#include "Eigen/Dense"

//Global Variables
const int ANTS = 4;					//Number of antennas
const int GRID_RES = 200;			//Resolution of Music Spectra		
const int NUM_SAMPLES = 1024;		//Samples per autocorrelation matrix
const int FRAC_BITS = 43;			//Fractional bits for auto corr numbers
const int TAGS = 2;					//Number of tags

//Function Prototypes

int main(){
	float freq = 915*pow(10,6);			//Frequency of RFID Tags
	float lambda = 3*pow(10,8)/freq;	//Wavelength

	float array1[ANTS][3];		//ANTS antennas, (x,y,z) coordinates
	float array2[ANTS][3];		//^ but for array 2
	array1[0][0] = 0.0;
	array1[0][1] = 0.0;
	array1[0][2] = 0.0;

	array1[1][0] = 1.0;
	array1[1][1] = 0.0;
	array1[1][2] = 0.0;

	array1[2][0] = 0.0;
	array1[2][1] = 1.0;
	array1[2][2] = 0.0;

	array1[3][0] = 1.0;
	array1[3][1] = 1.0;
	array1[3][2] = 0.0;

	float *a1, *a2;
	a1 = new float[ANTS*3];
	a1 = &array1[0][0];
	a2 = &array2[0][0];

	cout << "Ant_locs: " << &a1 << endl;
	cout << "ant_locs: " << &a1 << endl;

	//Find centers of array 1 and 2
	float *center1, *center2;
	center1 = new float[3];
	center2 = new float[3];

//	centroid(a1,ANTS,center1);
//	centroid(a2,ANTS,center2);

	float _Complex *R;
	R =  new float _Complex[ANTS*ANTS];
	float valsR[] = {.6079,-.174,.0195,.3498,.6065,.1417,-.0273,.5921,-.2394,.6018};
	float valsC[] = {0,0.0893,-.2822,.3478,0,-.2195,-.0814,0,.2838,0};
//	float valsR[] = {1,2,3,4,5,6,7,8,9,10};
//	float valsC[] = {0,2,3,4,0,6,7,0,9,0};

	vec2autocorr(valsR,valsC,ANTS,R);

/*
	for (int i=0; i<ANTS*ANTS; i++){
		printf("%f %fj,  ", creal(*(R+i)), cimag(*(R+i)));
		if ((i+1)%ANTS == 0){
			printf("\n");
		}
	}
*/

	printf("\nEigen Stuff: \n");

	MatrixXcf *eigvecs;
	MatrixXf *eigvals;
	eigvecs = new MatrixXcf;
	eigvals = new MatrixXf;


	autocorr2eig(R, ANTS, eigvecs, eigvals);

	MatrixXcf *subspace;
	subspace = new MatrixXcf;

	subspaceMat(eigvals, eigvecs, TAGS, ANTS, subspace);

	MatrixXf *S_music, *thetas, *phis;
	S_music = new MatrixXf;
	thetas = new MatrixXf;
	phis = new MatrixXf;

	musicSpectrum(subspace, ANTS, a1, GRID_RES, S_music, thetas, phis);

	return 0;
}













