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
	
const std::complex<float> If(0.0f,1.0f);
const float oneNorm = 5.0*M_PI/180.0;

void 	ang2loc(float *antloc1, float th1, float ph1, float *antloc2, float th2, float ph2, int ants, float *dist, Vector3f *midpoint);
void	autocorr2eig(float _Complex *R, int ants, MatrixXcf *eigmat, MatrixXf *eigvals);	
void	 bestLocal(float *antloc1, float *thLocs1, float *phLocs1, float *antloc2, float *thLocs2, float *phLocs2, int ants, int tags, Vector3f *locations);

int 	bitVal(char bit);
void	centroid(float *points, int ants, float *center);
bool 	comparator(pair <float, int> p1, pair <float, int> p2);
void 	findPeaks(MatrixXf *S_music, MatrixXf *th, MatrixXf *ph, int gridRes, int tags, float *thetas, float *phis);

void 	musicSpectrum(MatrixXcf *subspace, int ants, float *antPos, int gridRes, float lambda, MatrixXf *S_music, MatrixXf *thetas, MatrixXf *phis);
float	str2val(char bin[], int fracBits);
void 	subspaceMat(MatrixXf *eigvals, MatrixXcf *eigvecs, int tags, int ants, MatrixXcf *subspace);
void 	vec2autocorr(float valsReal[], float valsComp[], int ants, float _Complex *R);



#define __MUSIC_FUNC_H
#endif
