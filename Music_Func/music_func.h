#ifndef __MUSIC_FUNC_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
 	
int 	bitVal(char bit);
void	centroid(float *points, int ants, float *center);	
float	str2val(char bin[], int fracBits);
void 	vec2autocorr(float valsReal[], float valsComp[], int ants, float _Complex *R);

#define __MUSIC_FUNC_H
#endif
