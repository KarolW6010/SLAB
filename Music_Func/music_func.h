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

using namespace std;
using Eigen::MatrixXcf;
using Eigen::VectorXcf;
using Eigen::VectorXd;
using Eigen::MatrixXf;
	
const std::complex<float> If(0.0f,1.0f);

void	autocorr2eig(float _Complex *R, int ants, MatrixXcf *eigmat, MatrixXf *eigvals);	
int 	bitVal(char bit);
void	centroid(float *points, int ants, float *center);
bool 	comparator(pair <float, int> p1, pair <float, int> p2);
int 	findMax(MatrixXf *S_music, int gridres, float &max);
float	findPeakRec(MatrixXf *S_music, int gridRes, int mid);

void 	musicSpectrum(MatrixXcf *subspace, int ants, float *antPos, int gridRes, MatrixXf *S_music, MatrixXf *thetas, MatrixXf *phis);
float	str2val(char bin[], int fracBits);
void 	subspaceMat(MatrixXf *eigvals, MatrixXcf *eigvecs, int tags, int ants, MatrixXcf *subspace);
void 	vec2autocorr(float valsReal[], float valsComp[], int ants, float _Complex *R);



#define __MUSIC_FUNC_H
#endif
