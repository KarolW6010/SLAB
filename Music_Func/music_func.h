#ifndef __MUSIC_FUNC_H

#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "Eigen/Dense"
#include <iostream>
#include <utility>
#include "Eigen/Geometry"
#include <cfloat>
#include <vector>

using namespace std;
using Eigen::MatrixXcf;
using Eigen::VectorXcf;
using Eigen::VectorXd;
using Eigen::MatrixXf;
using Eigen::Vector3f;

//Global Variables
const int ANTS = 4;							//Number of antennas
extern int GRID_RES;						//Resolution of Music Spectra
const int NUM_SAMPLES = 1024;				//Samples per autocorrelation matrix
const int FRAC_BITS = 43;					//Fractional nits for autocorr numbers
extern int TAGS;							//Number of tags (can be changed)
const float FREQ = 915000000;				//Tag frequency
const float DOL = 0.5;						//Distance over lambda
const float LAMBDA = 300000000/FREQ;		//Wavelength
const float D = DOL*LAMBDA;					//Distance
const std::complex<float> If(0.0f,1.0f);	//Imaginary unit
const float oneNorm = 5.0*M_PI/180.0;		//Min distance between peaks

//Functions in alphabetical order
void 	ang2loc(float *antloc1, float th1, float ph1, float *antloc2, float th2,
			    float ph2, float *dist, Vector3f *midpoint);
void	autocorr2eig(float _Complex *R, MatrixXcf *eigmat, MatrixXf *eigvals);	
void	bestLocal(float *antloc1, float *thLocs1, float *phLocs1, float *antloc2,
				  float *thLocs2, float *phLocs2, Vector3f *locations);
int 	bitVal(char bit);
void	centroid(float *points, int ants, float *center);
bool 	comparator(pair <float, int> p1, pair <float, int> p2);
void 	findPeaks(MatrixXf *S_music, MatrixXf *th, MatrixXf *ph, float *thetas,
				  float *phis);
void 	musicSpectrum(MatrixXcf *subspace, float *antPos, MatrixXf *S_music,
					  MatrixXf *thetas, MatrixXf *phis);
void 	setConst(int grid_res, int tags);
float	str2val(char bin[], int fracBits);
void 	subspaceMat(MatrixXf *eigvals, MatrixXcf *eigvecs, MatrixXcf *subspace);
void 	vec2autocorr(float valsReal[], float valsComp[], float _Complex *R);

//Big func
void	Rs2Loc(int tags, int grid_res,
			   char r1_0[], char r1_1[], char r1_2[], char r1_3[], char r1_4[],
			   char r1_5[], char r1_6[], char r1_7[], char r1_8[], char r1_9[], 
			   char r2_0[], char r2_1[], char r2_2[], char r2_3[], char r2_4[],
			   char r2_5[], char r2_6[], char r2_7[], char r2_8[], char r2_9[],
			   Vector3f *locs); 

#define __MUSIC_FUNC_H
#endif
