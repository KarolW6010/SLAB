#include <stdio.h>
#include <string.h>
#include <math.h>

float str2val(char bin[], int fracBits);	//Convert binary strings to floats
int bitVal(char bit);				//Convert binary char to int

int main(){
	int fb = 2;		//Fractional bits
	char magic[10];		//Number to convert

	printf("Input 10 binary digits\n");
	scanf("%s",magic);
	float val = str2val(magic,fb);

	printf("%f\n",val);
	
	return 0;	
}

int bitVal(char bit){
	int val;
	if(bit == '0'){
		val = 0;
	}
	else{
		val = 1;
	}
	
	return val;
}

float str2val(char bin[], int fracBits){
	float val = 0;		//Value of binary string
	int len = strlen(bin);	//How many bits to convert
	int bit;		//bit value: either 0 or 1

	if(bin[0] == '0'){	//Positive Value
		for(int i = len-1; i >= 1; i--){
			bit = bitVal(bin[i]);
			val += bit*pow(2,len-1-i);
		}
	}
	else{			//Negative Value
		for(int i = len-1; i>=1; i--){
			bit = -(bitVal(bin[i])-1);	//2s complement flip bits
			val += bit*pow(2,len-1-i);
		}
		val++;		//2s complement increment
		val *= -1;	//Negate
	}

	val/= pow(2,fracBits);	//Scale appropriately
	
	return val;
}
