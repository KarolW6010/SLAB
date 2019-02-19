#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "music_func.h"

//Global Variables
const int ANTS = 4;			//Number of antennas
const int GRID_RES = 200;		//Resolution of Music Spectra		
const int NUM_SAMPLES = 1024;		//Samples per autocorrelation matrix

//Function Prototypes
//void centroid(float points[ANTS][3], float *center);	//Finds centroid of the points

int main(){
	float freq = 915*pow(10,6);		//Frequency of RFID Tags
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
	a1 = &array1[0][0];
	a2 = &array2[0][0];

	//Find centers of array 1 and 2
	float *center1, *center2;
	center1 = malloc(3*sizeof(float));
	center2 = malloc(3*sizeof(float));
//	centroid(array1, center1);
//	centroid(array2, center2);

	centroid(a1,ANTS,center1);
	centroid(a2,ANTS,center2);

	for (int i=0; i<3; i++){
		printf("%f\n",*(center1+i));
	}

	return 0;
}

//void centroid(float points[ANTS][3], float *center){
//	for (int i=0; i<ANTS; i++){
//		for (int j=0; j<3; j++){
//			*(center+j) += points[i][j]/(double)ANTS;
//		}
//	}
//
//	printf("X: %f, Y: %f, Z: %f\n", *(center), *(center+1), *(center+2));
//}
