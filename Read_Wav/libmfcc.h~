/*
 * libmfcc.h - Header for libMFCC
 * Copyright (c) 2010 Jeremy Sawruk
 *
 * This code is released under the MIT License. 
 * For conditions of distribution and use, see the license in LICENSE
 */
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

// Returns the specified (mth) MFCC
float GetCoefficient(float* spectralData, unsigned int samplingRate, unsigned int NumFilters, unsigned int binSize, unsigned int m);

/* 
 * Computes the Normalization Factor (Equation 6)
 * Used for internal computation only - not to be called directly
 */
float NormalizationFactor(int NumFilters, int m);

// Compute the filter parameter for the specified frequency and filter bands (For internal computation only - not the be called directly)
float GetFilterParameter(unsigned int samplingRate, unsigned int binSize, unsigned int frequencyBand, unsigned int filterBand);

// Compute the band-dependent magnitude factor for the given filter band (For internal computation only - not the be called directly)
float GetMagnitudeFactor(unsigned int filterBand);

// Compute the center frequency (fc) of the specified filter band (l) (For internal computation only - not the be called directly)
float GetCenterFrequency(unsigned int filterBand);


/* precalculation the filterbank */
void PreCalcFilterBank(float Array[NUMFILTERBANK][NUMBINHALF], float fNorm[NUMFILTERBANK], int numBin, int numBank);

void MFCC(complex *Sample, float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm,float *Mel);
//void MFCC(float *x, float *y, float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm, float *Mel);

void GetMFCC(float* spectralData,float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm, float *out);

void fft( complex *v, int n, complex *tmp );

void ifft( complex *v, int n, complex *tmp );

void abs_complex(complex * In, float *Out, int Size);

short _FFT(short int dir,long m,float *x,float *y);
