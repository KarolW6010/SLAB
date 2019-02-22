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

void autocorr2eig(float *R, int ants, MatrixXcf *eigmat, VectorXcf *eigvals){
	MatrixXcf Rmat(ants,ants);

	//Place R values into appropriate data type
	for (int i=0; i<ants*ants; i++){
		row = i/ants;
		col = i%ants;
		Rmat(row,col) = *(R+i);
	}

	Eigen::SelfAdjointEigenSolver<MatrixXcf> es(Rmat);
	*eigvals = es.eigenvalues();
	*eigmat = es.eigenvectors();

	std::cout << "Eigen Values: " << *eigvals << std::endl;
	std::cout << "Eigen Vectors: " << *eigmat << std::endl;
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

void vec2autocorr(float valsReal[], float valsComp[],  int ants, float complex *R){
/* Returns autocorrelation matrix formed from the values in vals.
 *
 * Inputs:
 * 	valsReal: Array of the real portion of values containing the upper traingular portion
 * 		  [r11, r12, r13, r22, r23, r33] are the indexes for the 3x3 R case.
 * 	valsComp: Array of the complex portion of values containing the upper traingular portion
 * 	ants	: Number of antennas (aka size of square matrix)
 *
 * Outputs: 
 * 	*R	: Pointer to a 2D array of complex values
 */	
	float complex temp;
	for (int i=0; i<ants; i++){
		for (int j=0; j<ants; j++){
			if (j>=i){	//Upper triangular
				int ind = i*ants - i*(i+1)/2 + j;
				temp = valsReal[ind] + I*valsComp[ind];	//Form the complex number. I is imag unit
				*(R + ants*i + j) = temp;
			}	
			else{		//Lower triangular 
				int ind = j*ants - j*(j+1)/2 + i;
				temp = valsReal[ind] - I*valsComp[ind];	//Form the complex number
				*(R + ants*i + j) = temp;
			}
		}
   	}
}




