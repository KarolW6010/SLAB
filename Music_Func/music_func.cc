#include "music_func.h"

int GRID_RES;
int TAGS;

void ang2loc(float *antloc1, float th1, float ph1, float *antloc2, float th2,
             float ph2, float *dist, Vector3f *midpoint){
/* Reutrns the minimum distance between the two rays pertruding from the center
   of the antenna arrays specified by the angles.

Inputs:
	*antloc	: [ANTS][3] Contains the coordiantes of an antenna array
	th		: theta value corresponding to array
	ph 		: phi value corresponding to array

Outputs:
	*dist	: Minimum distance between the two rays
	*midpoint	: Midpoint between the two lines
*/

	float *c1, *c2;
	c1 = new float[3];	//Center of array 1
	c2 = new float[3];	//Center of array 2

	centroid(antloc1, ANTS, c1);
	centroid(antloc2, ANTS, c2);

	//Store centers in vector objects
	Vector3f center1, center2;
	for(int i=0; i<3; i++){
		center1(i) = *(c1+i);
		center2(i) = *(c2+i);
	}

//	cout << "\nAngles1: Th = " << th1 << ", Ph = " << ph1 << endl;
//	cout << "\nAngles2: Th = " << th2 << ", Ph = " << ph2 << endl;

	Vector3f dir1, dir2; 	//Direction of ray
	dir1(0) = cos(th1);		//x component
	dir1(1) = sin(th1);		//y component
	dir1(2) = cos(ph1)/sin(ph1);	//z component
	dir2(0) = cos(th2);		//x component
	dir2(1) = sin(th2);		//y component
	dir2(2) = cos(ph2)/sin(ph2);	//z component

//	cout << "\nDirection 1:\n" << dir1 << endl;
//	cout << "\nDirection 2:\n" << dir2 << endl;

	Vector3f nor1, nor2, distNor;		//Normal vectors
	nor1 = dir1.cross(dir2.cross(dir1));		//For center calc
	nor2 = dir2.cross(dir1.cross(dir2));		//For center calc
	distNor = (dir1.cross(dir2)).normalized();	//For distance calc
	
//	cout << "\nNormal 1:\n" << nor1 << endl;
//	cout << "\nNormal 2:\n" << nor2 << endl;
//	cout << "\nDistance Normal:\n" << distNor << endl;

	float distance = abs(distNor.dot(center1-center2));		//Distance between 2 skew lines
	*dist = distance;

	Vector3f point1,  point2; 	//Points on line1 and line2 closest to line 2 and line 1 respectively
	point1 = center1 + ((center2-center1).dot(nor2)/dir1.dot(nor2))*dir1;
	point2 = center2 + ((center1-center2).dot(nor1)/dir2.dot(nor1))*dir2;

//	cout << "\nPoint 1:\n" << point1 << endl;
//	cout << "\nPoint 2:\n" << point2 << endl;

	Vector3f mid;
	mid = (point1+point2)/(float)(2);
	*midpoint = mid;
}

void autocorr2eig(float _Complex *R, MatrixXcf *eigmat, MatrixXf *eigvals){
/*Takes an autocorrelation matrix and returns the eigenvectors and eigenvalues.

Inputs:
	*R 		: Pointer to R matrix values
	*eigmat	: Pointer to eigenvectors stored in matrix
	*eigvals: Pointer to eigenvalues sotred in matrix (row vector)
*/
	MatrixXcf Rmat(ANTS,ANTS);

	//Place R values into appropriate data type
	for (int i=0; i<ANTS*ANTS; i++){
		int row = i/ANTS;
		int col = i%ANTS;
		Rmat(row,col) = *(R+i);
	}

	Eigen::SelfAdjointEigenSolver<MatrixXcf> es(Rmat);
	*eigvals = es.eigenvalues();
	*eigmat = es.eigenvectors();

//	std::cout << "Eigen Values: \n" << *eigvals << std::endl;
//	std::cout << "\nEigen Vectors: \n" << *eigmat << std::endl;
//	std::cout << "\nLambda*v: \n" << (*eigvals)(0)*((*eigmat).col(0)) << std::endl;
//	std::cout << "\nR*eig: \n" << (Rmat)*((*eigmat).col(0)) << std::endl;
}

void bestLocal(float *antloc1, float *thLocs1, float *phLocs1, float *antloc2,
			   float *thLocs2, float *phLocs2, Vector3f *locations){
//Takes the thetas and phis from each array and finds the best matching to 
//localize the tags.
	vector <int> avail;		//Unpaired AOA indices
	pair <int,int> matches[TAGS];

	for(int i=0; i<TAGS; i++){
		matches[i].first = i;
		avail.push_back(i);
	}

	float minDist = FLT_MAX;	//Smallest distance found so far
	float *dist;				//Current distance found
	dist = new float;		
	Vector3f *midpoint;			//Current localized point
	midpoint = new Vector3f;
	Vector3f *loc;				//Best localized point
	loc = new Vector3f;
	int bestInd;				//Index of best match

	Vector3f locals[TAGS];		//Array of tags locations

	for(int i=0; i<TAGS; i++){
		for(int j=0; j<TAGS-i; j++){
			ang2loc(antloc1,*(thLocs1+i),*(phLocs1+i),antloc2,
				*(thLocs2+avail.at(j)),*(phLocs2+avail.at(j)),dist,midpoint);
			if(*dist < minDist){
				minDist = *dist;
				bestInd = j;
				loc = midpoint;

			//	cout << "MIDPOINT " << *midpoint << endl;
			};
		}
		matches[i].second = avail.at(bestInd);
		avail.erase(avail.begin()+bestInd);
		locals[i] = *loc;
	}

	*locations = locals[0];

	for(int i=0; i<TAGS; i++){
		cout << "Matches: " << matches[i].first << " and "
		     << matches[i].second << endl;
	}
}

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
 	*center = 0.0;
	*(center+1) = 0.0;
	*(center+2) = 0.0;
// 	cout << "AAAA: " << &points << endl;
	for (int i=0; i<ants; i++){
		for (int j=0; j<3; j++){
				*(center+j) += *(points + 3*i + j)/(double)ants;
		}
	}
	
//	printf("X: %f, Y: %f, Z: %f\n", *(center), *(center+1), *(center+2));
}

bool comparator(pair <float, int> p1, pair <float, int> p2){
	return (p1.first > p2.first);
}

void findPeaks(MatrixXf *S_music, MatrixXf *th, MatrixXf *ph, float *thetas,
			   float *phis){
/* Find the theta and phi values of the peaks of the Music Spectrum.
 * 
 * Inputs: 
 * 	*S_music: Music spectrum
 *  *th		: Theta value grid
 *  *ph 	: Phi value grid
 * 
 * Outputs:
 * 	*thetas	: theta values of the peaks
 *  *phis	: phi values of the peals  
 */
	if(TAGS == 1){
		float thMax;
		float phMax;
		float maxPeak = FLT_MIN;

		for(int i=0; i<GRID_RES; i++){
			for(int j=0; j<GRID_RES; j++){
				if((*S_music)(i,j) > maxPeak){
					thMax = (*th)(i,j);
					phMax = (*ph)(i,j);
					maxPeak = (*S_music)(i,j);
				}
			}
		}
		*thetas = thMax;
		*phis = phMax;
	
		cout << "Theta: " << thMax << ", Phi: " << phMax << endl;
	}
	else{
	
		bool isPeak = false;
		int im, ip, jm, jp, k;
		im = 0;		//Previous Row
		ip = 0;		//Next Row
		jm = 0;		//Previous Column
		jp = 0;		//Next Column
		k = 0;

		int gg = GRID_RES*GRID_RES;
		float peaks[gg];
		float thTemp[gg];
		float phTemp[gg];

		for(int i=0; i<GRID_RES; i++){		//Theta loop (rows constant)		
			//Handle edge cases for rows
			if(i==0){
				im = GRID_RES-1;		//Wrap around on the theta
				ip = 1;
			}
			else if(i == (GRID_RES-1)){
				im = i-1;			
				ip = 0;				//Wrap around on the theta
			}
			else{
				im = i-1;
				ip = i+1;
			}
		
			for(int j=0; j<GRID_RES; j++){	//Phi loop (cols constant)
				//Handle edge cases for rows
				if(j==0){
					jm = 0;
					jp = 1;
				}
				else if(j == (GRID_RES-1)){
					jm = j-1;
					jp = GRID_RES-1;
				}
				else{
					jm = j-1;
					jp = j+1;
				}

				//Check if current value is larger or equal to neighbors
				isPeak = ((*S_music)(im,j)<=(*S_music)(i,j))
					   & ((*S_music)(i,jm)<=(*S_music)(i,j))
					   & ((*S_music)(ip,j)<=(*S_music)(i,j))
					   & ((*S_music)(i,jp)<=(*S_music)(i,j));

				//Store peak value and corresponding theta and phi
				bool isok = true;
				float differ;
				if(isPeak){
					for(int n=0; n<k; n++){
	//					float tempth = thTemp[n];
	//					if(thTemp[n] > 0){
	//						tempth = thTemp[n] - M_PI;
	//					}
	//					else{
	//						tempth = thTemp[n] + M_PI;
	//					}
						differ = abs(thTemp[n]-(*th)(i,j)) + abs(phTemp[n]-(*ph)(i,j));
				//		cout << "Differ = " << differ << ", k = " << k << endl;
						isok = isok && (differ > oneNorm);
					}

					if((*th)(i,j) < (-M_PI + oneNorm)){
						for(int n=0; n<k; n++){
							differ = abs(thTemp[n] + (float)(2.0*M_PI) - (*th)(i,j)
										+ abs(phTemp[n] - (*ph)(i,j)));
							isok = isok && (differ > oneNorm);	
						}
					}
					else if((*th)(i,j) < (M_PI - oneNorm)){
						for(int n=0; n<k; n++){
							differ = abs(thTemp[n] - (float)(2.0*M_PI) - (*th)(i,j)
										+ abs(phTemp[n] - (*ph)(i,j)));
							isok = isok && (differ > oneNorm);
						}
					}

					if(isok){
				//		cout << "End my misery!\n";
						peaks[k] = (*S_music)(i,j);
						thTemp[k] = (*th)(i,j);
						phTemp[k] = (*ph)(i,j);

	//					if(thTemp[k] > 0){
	//						thTemp[k] = thTemp[k] - M_PI;
	//					}
	//					else{
	//						thTemp[k] = thTemp[k] + M_PI;
	//					}
						k++;
					}
				}

		//				cout << "Peak " << k <<" found! Val: " << peaks[k] << ", th: " << thTemp[k] << ", ph: " << phTemp[k] << endl;
			}
		}

		//Sort the peaks so that highest peak is first (Stronget Signal)
		pair <float, int> locs[k];
		for(int i=0; i<k; i++){
			locs[i].first = peaks[i];
			locs[i].second = i;
		}
		int n = sizeof(locs)/sizeof(locs[0]);
		std::sort(locs,locs+n,comparator);

		float allOrdTh[k], allOrdPh[k];
		for(int i=0; i<k; i++){
			allOrdTh[i] = thTemp[locs[i].second];
			allOrdPh[i] = phTemp[locs[i].second];
		}

		//The above two for loops organize all the peaks and corresponding in angles in
		//order of highest peaks to lowest peaks.


		/*
		for(int i=0; i<k; i++){
			cout << "value: " << locs[i].first << " , Position: " << locs[i].second << endl;
		}
		*/
		//Take the n largest peaks where n is the number of tags;
		float ordTh[TAGS];		//Ordered thetas
		float ordPh[TAGS];		//Ordered phis


		for(int i=0; i<TAGS; i++){
//			cout << "Location index: " << locs[i].second << endl;
//			cout <<"\t TH: " << thTemp[locs[i].second]*(180/M_PI) << endl;
//			cout <<"\t PH: " << phTemp[locs[i].second]*(180/M_PI) << endl;
			ordTh[i] = thTemp[locs[i].second];
			ordPh[i] = phTemp[locs[i].second];
		}

		*thetas = ordTh[0];
		*phis = ordPh[0];

		/*
		for(int i=0; i<tags; i++){
			cout << "ThPh" << i << ": " << ordTh[i] << ", " << ordPh[i] << endl;
		}
		*/
	}	
}

void musicSpectrum(MatrixXcf *subspace, float *antPos, MatrixXf *S_music, MatrixXf *thetas, MatrixXf *phis){
/*Calculates the music spectrum based on the noise subspace matrix

Inputs:	
	*subspace	: Noise subspace matrix
	*antPos		: [ants][3] pointer to antenna coordinates

Outputs:
	*S_music	: [gridRes][girdRes] matrix containing values of the spectrum
					theta in first dim, phi in second dim
	*thetas		: [gridRes][gridRes] matrix containing theta values
	*phis		: [gridRes][gridRes] matrix containing phis values
*/

	MatrixXcf music_spec(GRID_RES,GRID_RES);	//To be filled in and then sent back

	float *center;
	center = new float[3];
//	cout << "Antpos: " <<  &antPos << endl;
//	printf("Center %f, %f, %f\n", *center, *(center+1), *(center+2));
	centroid(antPos, ANTS, center);

//	printf("Center %f, %f, %f\n", *center, *(center+1), *(center+2));

	MatrixXf centered(3,ANTS);	//Shifted antPos by center
	for(int i=0; i<ANTS; i++){
		for(int j=0; j<3; j++){
			centered(j,i) = *(antPos+3*i+j) - *(center+j);
		}
	}

//	cout << "Ant_locs\n" <<  centered << endl;
	
	MatrixXf dir_vec(1,3);		//Directional Vector
	MatrixXcf steerRow(1,ANTS);	//Steering vector

	float th, ph;
	MatrixXf theta(GRID_RES,GRID_RES), phi(GRID_RES,GRID_RES);

	float kwav = 2.0*M_PI/LAMBDA;
	MatrixXcf temp1(1,1);
	for(int i=0; i<GRID_RES; i++){
		for(int j=0; j<GRID_RES; j++){
			th = -M_PI + i*(2*M_PI)/(float)GRID_RES;
			ph = j*(M_PI)/(float)(GRID_RES*2);
			theta(i,j) = th;
			phi(i,j) = ph;
			dir_vec(0,0) = sin(ph)*cos(th);		//ak x component
			dir_vec(0,1) = sin(ph)*sin(th);		//ak y component
			dir_vec(0,2) = cos(ph);				//ak z component

			steerRow = -((float)(kwav))*If*dir_vec*centered;
			for(int k=0; k<ANTS; k++){
				steerRow(0,k) = exp(steerRow(0,k))/((float)(sqrt(ANTS)));
			}

			temp1(0,0) =((steerRow.conjugate())*(*subspace)*(steerRow.transpose()))(0,0);
			music_spec(i,j) = (float)1.0/temp1(0,0);

		//	cout << "th: " << th << ", ph: " << ph << ", vect: " << steerRow << endl;
		}
	}

	*thetas = theta;
	*phis = phi;
	MatrixXf temp(GRID_RES,GRID_RES);
	
	//Take abs of music spec and print it
	for(int i=0; i<GRID_RES; i++){
		for(int j=0; j<GRID_RES; j++){
			temp(i,j) = ((float)20)*log10(abs(music_spec(i,j)));
			//temp(i,j) = abs(music_spec(i,j));
//			printf("%f\n",temp(i,j));
		}
	}
	*S_music = temp;
}

void setConst(int grid_res, int tags){
//Sets the grid resolution and number of tags that are being used.
	GRID_RES = grid_res;
	TAGS = tags;
}

float str2val(char bin[], int fracBits){
//Takes a string representing a 2s complement binary number and converts to a float.
	float val = 0;		//Value of binary string
	int len = strlen(bin)-1;	//How many bits to convert
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

void subspaceMat(MatrixXf *eigvals, MatrixXcf *eigvecs, MatrixXcf *subspace){
/*Takes eigenvectors and eigenvalues and returns the appropriate Subspace matrix.

Inputs:
	*eigvals	: Pointer to eigenvalues
	*eigvecs	: Pointer to eigenvector matrix

Ouputs:
	*subspace	: Pointer to noise subspace matrix
*/
	pair <float, int> lambs[ANTS];
	for(int i=0; i<ANTS; i++){
		lambs[i].first = (*eigvals)(i);
		lambs[i].second = i;
	}

	int n = sizeof(lambs)/sizeof(lambs[0]);
	std::sort(lambs,lambs+n,comparator);

/*
	for(int i=0; i<ants; i++){
		cout << "Lambs at " << i << ": First  " << lambs[i].first
			 << ", Second " << lambs[i].second << endl;
	}
*/
	MatrixXcf subs(ANTS, ANTS-TAGS);

	for(int i=0; i<ANTS-TAGS; i++){
		int ind = lambs[i+TAGS].second;
		subs.col(i) = (*eigvecs).col(ind);
	}

	(*subspace) = subs*subs.adjoint();

//	std::cout << "\nSubs :\n" << *subspace << std::endl;
}

void vec2autocorr(float valsReal[], float valsComp[], float _Complex *R){
/* Returns autocorrelation matrix formed from the values in vals.
 *
 * Inputs:
 * 	valsReal: Array of the real portion of values containing the upper
 			  traingular portion
 * 		  [r11, r12, r13, r22, r23, r33] are the indexes for the 3x3 R case.
 * 	valsComp: Array of the complex portion of values containing the upper
 			  traingular portion
 * 	ants	: Number of antennas (aka size of square matrix)
 *
 * Outputs: 
 * 	*R	: Pointer to a 2D array of complex values
 */	
	float _Complex temp;
	for (int i=0; i<ANTS; i++){
		for (int j=0; j<ANTS; j++){
			if (j>=i){	//Upper triangular
				int ind = i*ANTS - i*(i+1)/2 + j;
				temp = valsReal[ind] + I*valsComp[ind];	//Form the complex number. I is imag unit
				*(R + ANTS*i + j) = temp;
			}	
			else{		//Lower triangular 
				int ind = j*ANTS - j*(j+1)/2 + i;
				temp = valsReal[ind] - I*valsComp[ind];	//Form the complex number
				*(R + ANTS*i + j) = temp;
			}
		}
   	}
}


void	Rs2Loc(int tags, int grid_res,
			   char r1_0[], char r1_1[], char r1_2[], char r1_3[], char r1_4[],
			   char r1_5[], char r1_6[], char r1_7[], char r1_8[], char r1_9[], 
			   char r2_0[], char r2_1[], char r2_2[], char r2_3[], char r2_4[],
			   char r2_5[], char r2_6[], char r2_7[], char r2_8[], char r2_9[],
			   Vector3f *locs){
	setConst(grid_res, tags);

	int len = strlen(r1_0);

	//Break up char arrays into real and imaginary arrays
	char rr1_0[len/2+1], rc1_0[len/2+1];
	char rr1_1[len/2+1], rc1_1[len/2+1];
	char rr1_2[len/2+1], rc1_2[len/2+1];
	char rr1_3[len/2+1], rc1_3[len/2+1];
	char rr1_4[len/2+1], rc1_4[len/2+1];
	char rr1_5[len/2+1], rc1_5[len/2+1];
	char rr1_6[len/2+1], rc1_6[len/2+1];
	char rr1_7[len/2+1], rc1_7[len/2+1];
	char rr1_8[len/2+1], rc1_8[len/2+1];
	char rr1_9[len/2+1], rc1_9[len/2+1];

	char rr2_0[len/2+1], rc2_0[len/2+1];
	char rr2_1[len/2+1], rc2_1[len/2+1];
	char rr2_2[len/2+1], rc2_2[len/2+1];
	char rr2_3[len/2+1], rc2_3[len/2+1];
	char rr2_4[len/2+1], rc2_4[len/2+1];
	char rr2_5[len/2+1], rc2_5[len/2+1];
	char rr2_6[len/2+1], rc2_6[len/2+1];
	char rr2_7[len/2+1], rc2_7[len/2+1];
	char rr2_8[len/2+1], rc2_8[len/2+1];
	char rr2_9[len/2+1], rc2_9[len/2+1];
	
	//Fill the arrays
	for(int i=0; i<len/2; i++){
		rr1_0[i] = r1_0[i];		rc1_0[i] = r1_0[i+len/2];
		rr1_1[i] = r1_1[i];		rc1_1[i] = r1_1[i+len/2];
		rr1_2[i] = r1_2[i];		rc1_2[i] = r1_2[i+len/2];
		rr1_3[i] = r1_3[i];		rc1_3[i] = r1_3[i+len/2];
		rr1_4[i] = r1_4[i];		rc1_4[i] = r1_4[i+len/2];
		rr1_5[i] = r1_5[i];		rc1_5[i] = r1_5[i+len/2];
		rr1_6[i] = r1_6[i];		rc1_6[i] = r1_6[i+len/2];
		rr1_7[i] = r1_7[i];		rc1_7[i] = r1_7[i+len/2];
		rr1_8[i] = r1_8[i];		rc1_8[i] = r1_8[i+len/2];
		rr1_9[i] = r1_9[i];		rc1_9[i] = r1_9[i+len/2];

		rr2_0[i] = r2_0[i];		rc2_0[i] = r2_0[i+len/2];
		rr2_1[i] = r2_1[i];		rc2_1[i] = r2_1[i+len/2];
		rr2_2[i] = r2_2[i];		rc2_2[i] = r2_2[i+len/2];
		rr2_3[i] = r2_3[i];		rc2_3[i] = r2_3[i+len/2];
		rr2_4[i] = r2_4[i];		rc2_4[i] = r2_4[i+len/2];
		rr2_5[i] = r2_5[i];		rc2_5[i] = r2_5[i+len/2];
		rr2_6[i] = r2_6[i];		rc2_6[i] = r2_6[i+len/2];
		rr2_7[i] = r2_7[i];		rc2_7[i] = r2_7[i+len/2];
		rr2_8[i] = r2_8[i];		rc2_8[i] = r2_8[i+len/2];
		rr2_9[i] = r2_9[i];		rc2_9[i] = r2_9[i+len/2];
	}

	int t = len/2 + 1;
	rr1_0[t] = '\0';	rc1_0[t] = '\0';
	rr1_1[t] = '\0';	rc1_1[t] = '\0';
	rr1_2[t] = '\0';	rc1_2[t] = '\0';
	rr1_3[t] = '\0';	rc1_3[t] = '\0';
	rr1_4[t] = '\0';	rc1_4[t] = '\0';
	rr1_5[t] = '\0';	rc1_5[t] = '\0';
	rr1_6[t] = '\0';	rc1_6[t] = '\0';
	rr1_7[t] = '\0';	rc1_7[t] = '\0';
	rr1_8[t] = '\0';	rc1_8[t] = '\0';
	rr1_9[t] = '\0';	rc1_9[t] = '\0';

	rr2_0[t] = '\0';	rc2_0[t] = '\0';
	rr2_1[t] = '\0';	rc2_1[t] = '\0';
	rr2_2[t] = '\0';	rc2_2[t] = '\0';
	rr2_3[t] = '\0';	rc2_3[t] = '\0';
	rr2_4[t] = '\0';	rc2_4[t] = '\0';
	rr2_5[t] = '\0';	rc2_5[t] = '\0';
	rr2_6[t] = '\0';	rc2_6[t] = '\0';
	rr2_7[t] = '\0';	rc2_7[t] = '\0';
	rr2_8[t] = '\0';	rc2_8[t] = '\0';
	rr2_9[t] = '\0';	rc2_9[t] = '\0';

//	cout << "RR1_0 : " << rr1_0[0] << rr1_0[1] << rr1_0[2] << rr1_0[3] <<
//						  rr1_0[4] << endl;
//
//	cout << "RC1_0 : " << rc1_0[0] << rc1_0[1] << rc1_0[2] << rc1_0[3] <<
//						  rc1_0[4] << endl;

	//Convert character arrays to numbers
		//R1 real part
	float Rr1_0 = str2val(rr1_0, FRAC_BITS);
	float Rr1_1 = str2val(rr1_1, FRAC_BITS);
	float Rr1_2 = str2val(rr1_2, FRAC_BITS);
	float Rr1_3 = str2val(rr1_3, FRAC_BITS);
	float Rr1_4 = str2val(rr1_4, FRAC_BITS);
	float Rr1_5 = str2val(rr1_5, FRAC_BITS);
	float Rr1_6 = str2val(rr1_6, FRAC_BITS);
	float Rr1_7 = str2val(rr1_7, FRAC_BITS);
	float Rr1_8 = str2val(rr1_8, FRAC_BITS);
	float Rr1_9 = str2val(rr1_9, FRAC_BITS);
		//R1 imagianry part
	float Rc1_0 = str2val(rc1_0, FRAC_BITS);
	float Rc1_1 = str2val(rc1_1, FRAC_BITS);
	float Rc1_2 = str2val(rc1_2, FRAC_BITS);
	float Rc1_3 = str2val(rc1_3, FRAC_BITS);
	float Rc1_4 = str2val(rc1_4, FRAC_BITS);
	float Rc1_5 = str2val(rc1_5, FRAC_BITS);
	float Rc1_6 = str2val(rc1_6, FRAC_BITS);
	float Rc1_7 = str2val(rc1_7, FRAC_BITS);
	float Rc1_8 = str2val(rc1_8, FRAC_BITS);
	float Rc1_9 = str2val(rc1_9, FRAC_BITS);
		
	float valsR1[] = {Rr1_0, Rr1_1, Rr1_2, Rr1_3, Rr1_4,
					  Rr1_5, Rr1_6, Rr1_7, Rr1_8, Rr1_9};
	float valsC1[] = {Rc1_0, Rc1_1, Rc1_2, Rc1_3, Rc1_4,
					  Rc1_5, Rc1_6, Rc1_7, Rc1_8, Rc1_9};

		//R2 real part
	float Rr2_0 = str2val(rr2_0, FRAC_BITS);
	float Rr2_1 = str2val(rr2_1, FRAC_BITS);
	float Rr2_2 = str2val(rr2_2, FRAC_BITS);
	float Rr2_3 = str2val(rr2_3, FRAC_BITS);
	float Rr2_4 = str2val(rr2_4, FRAC_BITS);
	float Rr2_5 = str2val(rr2_5, FRAC_BITS);
	float Rr2_6 = str2val(rr2_6, FRAC_BITS);
	float Rr2_7 = str2val(rr2_7, FRAC_BITS);
	float Rr2_8 = str2val(rr2_8, FRAC_BITS);
	float Rr2_9 = str2val(rr2_9, FRAC_BITS);
		//R2 imaginary part
	float Rc2_0 = str2val(rc2_0, FRAC_BITS);
	float Rc2_1 = str2val(rc2_1, FRAC_BITS);
	float Rc2_2 = str2val(rc2_2, FRAC_BITS);
	float Rc2_3 = str2val(rc2_3, FRAC_BITS);
	float Rc2_4 = str2val(rc2_4, FRAC_BITS);
	float Rc2_5 = str2val(rc2_5, FRAC_BITS);
	float Rc2_6 = str2val(rc2_6, FRAC_BITS);
	float Rc2_7 = str2val(rc2_7, FRAC_BITS);
	float Rc2_8 = str2val(rc2_8, FRAC_BITS);
	float Rc2_9 = str2val(rc2_9, FRAC_BITS);
	
	float valsR2[] = {Rr2_0, Rr2_1, Rr2_2, Rr2_3, Rr2_4,
					  Rr2_5, Rr2_6, Rr2_7, Rr2_8, Rr2_9};
	float valsC2[] = {Rc2_0, Rc2_1, Rc2_2, Rc2_3, Rc2_4,
					  Rc2_5, Rc2_6, Rc2_7, Rc2_8, Rc2_9};
//	cout << valsR1[0] << " + j" << valsC1[0] << endl;

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
	thetas1 = new MatrixXf;			//Theta grid
	phis1 = new MatrixXf;			//Phi grid
	
	MatrixXf *S_music2, *thetas2, *phis2;
	S_music2 = new MatrixXf;		//Music Spectrum values
	thetas2 = new MatrixXf;			//Theta grid
	phis2 = new MatrixXf;			//Phi grid

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

	locs = locations;
}




