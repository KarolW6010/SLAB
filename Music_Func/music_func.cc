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

void ang2loc(float *antloc1, float th1, float ph1, float *antloc2, float th2, float ph2, int ants, float *dist, Vector3f *midpoint){
/* Reutrns the minimum distance between the two rays pertruding from the center of the antenna arrays specified by the angles.

Inputs:
	*antloc	: [ANTS][3] Contains the coordiantes of an antenna array
	th		: theta value corresponding to array
	ph 		: phi value corresponding to array
	ants	: number of antennas (ANTS)

Outputs:
	*dist	: Minimum distance between the two rays
	*midpoint	: Midpoint between the two lines
*/

	float *c1, *c2;
	c1 = new float[3];	//Center of array 1
	c2 = new float[3];	//Center of array 2

	centroid(antloc1, ants, c1);
	centroid(antloc2, ants, c2);

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

void autocorr2eig(float _Complex *R, int ants, MatrixXcf *eigmat, MatrixXf *eigvals){
/*Takes an autocorrelation matrix and returns the eigenvectors and eigenvalues.

Inputs:
	*R 		: Pointer to R matrix values
	ants	: Number of antennas
	*eigmat	: Pointer to eigenvectors stored in matrix
	*eigvals: Pointer to eigenvalues sotred in matrix (row vector)
*/
	MatrixXcf Rmat(ants,ants);

	//Place R values into appropriate data type
	for (int i=0; i<ants*ants; i++){
		int row = i/ants;
		int col = i%ants;
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

void bestLocal(float *antloc1, float *thLocs1, float *phLocs1, float *antloc2, float *thLocs2, float *phLocs2, int ants, int tags, Vector3f *locations){

	vector <int> avail;		//Unpaired AOA indices
	pair <int,int> matches[tags];

	for(int i=0; i<tags; i++){
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

	Vector3f locals[tags];		//Array of tags locations

	for(int i=0; i<tags; i++){
		for(int j=0; j<tags-i; j++){
			ang2loc(antloc1,*(thLocs1+i),*(phLocs1+i),antloc2,
				*(thLocs2+avail.at(j)),*(phLocs2+avail.at(j)),ants,dist,midpoint);
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

	for(int i=0; i<tags; i++){
		cout << "Matches: " << matches[i].first << " and " << matches[i].second << endl;
	}
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

void findPeaks(MatrixXf *S_music, MatrixXf *th, MatrixXf *ph, int gridRes, int tags, float *thetas, float *phis){
/* Find the theta and phi values of the peaks of the Music Spectrum.
 * 
 * Inputs: 
 * 	*S_music: Music spectrum
 *  *th		: Theta value grid
 *  *ph 	: Phi value grid
 *	gridRes	: Resolution of the grid
 *  tags	: Number of tags present (aka peaks to search for)
 * 
 * Outputs:
 * 	*thetas	: theta values of the peaks
 *  *phis	: phi values of the peals  
 */
	if(tags == 1){
		float thMax;
		float phMax;
		float maxPeak = FLT_MIN;

		for(int i=0; i<gridRes; i++){
			for(int j=0; j<gridRes; j++){
				if((*S_music)(i,j) > maxPeak){
					thMax = (*th)(i,j);
					phMax = (*ph)(i,j);
					maxPeak = (*S_music)(i,j);
				}
			}
		}
		*thetas = thMax;
		*phis = phMax;

	}
	else{
	
		bool isPeak = false;
		int im, ip, jm, jp, k;
		im = 0;		//Previous Row
		ip = 0;		//Next Row
		jm = 0;		//Previous Column
		jp = 0;		//Next Column
		k = 0;

		float peaks[gridRes*gridRes];
		float thTemp[gridRes*gridRes];
		float phTemp[gridRes*gridRes];

		for(int i=0; i<gridRes; i++){		//Theta loop (rows constant)		
			//Handle edge cases for rows
			if(i==0){
				im = gridRes-1;		//Wrap around on the theta
				ip = 1;
			}
			else if(i == (gridRes-1)){
				im = i-1;			
				ip = 0;				//Wrap around on the theta
			}
			else{
				im = i-1;
				ip = i+1;
			}
		
			for(int j=0; j<gridRes; j++){	//Phi loop (cols constant)
				//Handle edge cases for rows
				if(j==0){
					jm = 0;
					jp = 1;
				}
				else if(j == (gridRes-1)){
					jm = j-1;
					jp = gridRes-1;
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
		float ordTh[tags];		//Ordered thetas
		float ordPh[tags];		//Ordered phis


		for(int i=0; i<tags; i++){
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

void musicSpectrum(MatrixXcf *subspace, int ants, float *antPos, int gridRes, float lambda, MatrixXf *S_music, MatrixXf *thetas, MatrixXf *phis){
/*Calculates the music spectrum based on the noise subspace matrix

Inputs:	
	*subspace	: Noise subspace matrix
	ants		: Number of antennas
	*antPos		: [ants][3] pointer to antenna coordinates
	gridRes 	: Resolution of music spectrum
	lambda 		: Wavelength

Outputs:
	*S_music	: [gridRes][girdRes] matrix containing values of the spectrum
					theta in first dim, phi in second dim
	*thetas		: [gridRes][gridRes] matrix containing theta values
	*phis		: [gridRes][gridRes] matrix containing phis values
*/

	MatrixXcf music_spec(gridRes,gridRes);	//To be filled in and then sent back

	float *center;
	center = new float[3];
//	cout << "Antpos: " <<  &antPos << endl;
//	printf("Center %f, %f, %f\n", *center, *(center+1), *(center+2));
	centroid(antPos, ants, center);

//	printf("Center %f, %f, %f\n", *center, *(center+1), *(center+2));

	MatrixXf centered(3,ants);	//Shifted antPos by center
	for(int i=0; i<ants; i++){
		for(int j=0; j<3; j++){
			centered(j,i) = *(antPos+3*i+j) - *(center+j);
		}
	}

//	cout << "Ant_locs\n" <<  centered << endl;
	
	MatrixXf dir_vec(1,3);		//Directional Vector
	MatrixXcf steerRow(1,ants);	//Steering vector

	float th, ph;
	MatrixXf theta(gridRes,gridRes), phi(gridRes,gridRes);

	float kwav = 2.0*M_PI/lambda;
	MatrixXcf temp1(1,1);
	for(int i=0; i<gridRes; i++){
		for(int j=0; j<gridRes; j++){
			th = -M_PI + i*(2*M_PI)/(float)gridRes;
			ph = j*(M_PI)/(float)(gridRes*2);
			theta(i,j) = th;
			phi(i,j) = ph;
			dir_vec(0,0) = sin(ph)*cos(th);		//ak x component
			dir_vec(0,1) = sin(ph)*sin(th);		//ak y component
			dir_vec(0,2) = cos(ph);				//ak z component

			steerRow = -((float)(kwav))*If*dir_vec*centered;
			for(int k=0; k<ants; k++){
				steerRow(0,k) = exp(steerRow(0,k))/((float)(sqrt(ants)));
			}

			temp1(0,0) =((steerRow.conjugate())*(*subspace)*(steerRow.transpose()))(0,0);
			music_spec(i,j) = (float)1.0/temp1(0,0);

		//	cout << "th: " << th << ", ph: " << ph << ", vect: " << steerRow << endl;
		}
	}

	*thetas = theta;
	*phis = phi;
	MatrixXf temp(gridRes,gridRes);
	
	//Take abs of music spec and print it
	for(int i=0; i<gridRes; i++){
		for(int j=0; j<gridRes; j++){
			temp(i,j) = ((float)20)*log10(abs(music_spec(i,j)));
			//temp(i,j) = abs(music_spec(i,j));
//			printf("%f\n",temp(i,j));
		}
	}
	*S_music = temp;
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

void subspaceMat(MatrixXf *eigvals, MatrixXcf *eigvecs, int tags, int ants,  MatrixXcf *subspace){
/*Takes eigenvectors and eigenvalues and returns the appropriate Subspace matrix.

Inputs:
	*eigvals	: Pointer to eigenvalues
	*eigvecs	: Pointer to eigenvector matrix
	tags		: Number of tags
	ants		: Number of antennas

Ouputs:
	*subspace	: Pointer to noise subspace matrix
*/
	pair <float, int> lambs[ants];
	for(int i=0; i<ants; i++){
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
	MatrixXcf subs(ants, ants-tags);

	for(int i=0; i<ants-tags; i++){
		int ind = lambs[i+tags].second;
		subs.col(i) = (*eigvecs).col(ind);
	}

	(*subspace) = subs*subs.adjoint();

//	std::cout << "\nSubs :\n" << *subspace << std::endl;
}

void vec2autocorr(float valsReal[], float valsComp[],  int ants, float _Complex *R){
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
	float _Complex temp;
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




