/*
 Phan Le Son
 plson03@gmail.com 
 */

//#pragma once

typedef struct 
{
    float Re;
    float Im;
} complex;

typedef float real;

#define PI 3.14159265358979323846264338327
#define NUMFILTERBANK 64
#define NUMBINHALF    512 
#define FS            16000


float NormalizationFactor(int NumFilters, int m);

float GetFilterParameter(float currVal, unsigned int filterBand , float * Mel_Fre);

float GetCenterFrequency(unsigned int filterBand, int fs, int numPointBank);

/* precalculation the filterbank */
void PreCalcFilterBank(float Array[NUMFILTERBANK][NUMBINHALF], float fNorm[NUMFILTERBANK], int numBin, int numBank);

void MFCC(complex *Sample, float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm,float *Mel);

void GetMFCC(float* spectralData,float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm, float *out);

void EnergyFac(float FilterBank[NUMFILTERBANK][NUMBINHALF], float fEnergyFac[NUMFILTERBANK]);

void fft( complex *v, int n, complex *tmp );

void ifft( complex *v, int n, complex *tmp );

void abs_complex(complex * In, float *Out, int Size);

void lifter(float * Cepstra, unsigned char Len);

short _FFT(short int dir,long m,float *x,float *y);

/* Hanning window */
void Window (float *FIRCoef);

void Preemphasis(complex *Data);

