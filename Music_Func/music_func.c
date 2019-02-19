#include "music_func.h"

int bitVal(char bit){
//Takes a '0' or '1' char and returns as int.
	int val;
	if(bit == '0'){
		val = 0;
	}
	else{
		val = 1;
	}
	return val;
}


void centroid(float *points, int ants, float *center){
/* Finds the centroid of the points.
 *
 * Inputs:
 * 	*points : Pointer to a 2D array of size [ants][3] where each row contains (x,y,z) coordinate
 * 	ants	: Number of antennas to find the centroid of
 *
 * Outputs:
 * 	*center : Pointer to a 1D array of size [3] containing (x,y,z) coordinates of centroid
 */
	for (int i=0; i<ants; i++){
		for (int j=0; j<3; j++){
				*(center+j) += *(points + 3*i + j)/(double)ants;
		}
	}
	
	printf("X: %f, Y: %f, Z: %f\n", *(center), *(center+1), *(center+2));
}

float str2val(char bin[], int fracBits){
//Takes a string representing a 2s complement binary number and converts to a float.
	float val = 0;		//Value of binary string
	int len = strlen(bin);	//How many bits to convert
	int bit;		//Bit value: either 0 or 1

	if(bin[0] == '0'){	//Positive value
		for (int i=len-1; i>=1; i--){
			bit = bitVal(bin[i]);
			val += bit*pow(2,len-1-i);
		}
	}
	else{			//Negative value
		for (int i=len-1; i>=1; i--){
			bit = -(bitVal(bin[i])-1);	//2s complement flip bits
			val += bit*pow(2,len-1-i);
		}
		val++;		//2s complement increment
		val *= -1;	//Negate
	}

	val /= pow(2,fracBits); //Scale appropriately
	return val;
}







