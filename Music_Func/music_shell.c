#include <stdio.h>
#include <math.h>

void centroid(float points[][3]);	//Finds centroid of the points

int ANTS = 4;				//Number of antennas

int main(){
	float freq = 915*pow(10,6);		//Frequency of RFID Tags
	float lambda = 3*pow(10,8)/freq;	//Wavelength

	float array1[ANTS][3];		//ants antennas, (x,y,z) coordinates

	int grid_res = 200;		//Resolution of Music Spectra
	int num_samples = 1024;		//Samples per autocorrelation matrix	
	
	float center1[3] = centroid(array1);
	return 0;
}

void centroid(float points[][3]){

	float ans[2];
	return ans;
}
