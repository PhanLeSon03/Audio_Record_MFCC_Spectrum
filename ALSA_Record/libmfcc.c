/*
 * libmfcc.c - Code implementation for libMFCC
 * Copyright (c) 2010 Jeremy Sawruk
 *
 * This code is released under the MIT License. 
 * For conditions of distribution and use, see the license in LICENSE
 */
/*
 Phan Le Son
 plson03@gmail.com 
 */

#include <math.h>
#include "libmfcc.h"
#include <stdio.h>
#include <stdlib.h>
 
float FilterBank[NUMFILTERBANK][NUMBINHALF];
float fNorm[NUMFILTERBANK];
float result[NUMFILTERBANK];
complex tmpBuff[2*NUMBINHALF];

/* 
 * Computes the specified (mth) MFCC
 *
 * spectralData - array of floats containing the results of FFT computation. This data is already assumed to be purely real
 * samplingRate - the rate that the original time-series data was sampled at (i.e 44100)
 * NumFilters - the number of filters to use in the computation. Recommended value = 48
 * binSize - the size of the spectralData array, usually a power of 2
 * m - The mth MFCC coefficient to compute
 *
 */

/* 
 * Computes the specified (mth) MFCC
 *
 * spectralData - array of floats containing the results of FFT computation. This data is already assumed to be purely real
 * samplingRate - the rate that the original time-series data was sampled at (i.e 44100)
 * NumFilters - the number of filters to use in the computation. Recommended value = 48
 * binSize - the size of the spectralData array, usually a power of 2
 * m - The mth MFCC coefficient to compute
 *
 */
float GetCoefficient(float* spectralData, unsigned int samplingRate, unsigned int NumFilters, unsigned int binSize, unsigned int m)
{
	float result = 0.0f;
	float outerSum = 0.0f;
	float innerSum = 0.0f;
	unsigned int k, l;

	// 0 <= m < L
	if(m >= NumFilters)
	{
		// This represents an error condition - the specified coefficient is greater than or equal to the number of filters. The behavior in this case is undefined.
		return 0.0f;
	}

	result = NormalizationFactor(NumFilters, m);

	
	for(l = 1; l <= NumFilters; l++) //48
	{
		// Compute inner sum
		innerSum = 0.0f;
		for(k = 0; k < binSize - 1; k++)  //1024
		{
			innerSum += fabs(spectralData[k] * GetFilterParameter(samplingRate, binSize, k, l));
		}

		if(innerSum > 0.0f)
		{
			innerSum = log(innerSum); // The log of 0 is undefined, so don't use it
		}

		innerSum = innerSum * cos(((m * PI) / NumFilters) * (l - 0.5f));

		outerSum += innerSum;
	}

	result *= outerSum;

	return result;
}


/* 
 * Computes the Normalization Factor (Equation 6)
 * Used for internal computation only - not to be called directly
 */
float NormalizationFactor(int NumFilters, int m)
{
	float normalizationFactor = 0.0f;

	if(m == 0)
	{
		normalizationFactor = sqrt(1.0f / NumFilters);
	}
	else 
	{
		normalizationFactor = sqrt(2.0f / NumFilters);
	}
	
	return normalizationFactor;
}

/* 
 * Compute the filter parameter for the specified frequency and filter bands (Eq. 2)
 * Used for internal computation only - not the be called directly
 */
float GetFilterParameter(unsigned int samplingRate, unsigned int binSize, unsigned int frequencyBand, unsigned int filterBand)
{
	float filterParameter = 0.0f;

	float boundary = (frequencyBand * samplingRate) / binSize;		        // k * Fs / N
	float prevCenterFrequency = GetCenterFrequency(filterBand - 1);		// fc(l - 1) etc.
	float thisCenterFrequency = GetCenterFrequency(filterBand);
	float nextCenterFrequency = GetCenterFrequency(filterBand + 1);

    //printf(" %f - %f - %f - %f \n",boundary, prevCenterFrequency, thisCenterFrequency, nextCenterFrequency);
	if(boundary >= 0 && boundary < prevCenterFrequency)
	{
		filterParameter = 0.0f;
	}
	else if(boundary >= prevCenterFrequency && boundary < thisCenterFrequency)
	{
		filterParameter = (boundary - prevCenterFrequency) / (thisCenterFrequency - prevCenterFrequency);
		//filterParameter *= GetMagnitudeFactor(filterBand);
	}
	else if(boundary >= thisCenterFrequency && boundary < nextCenterFrequency)
	{
		filterParameter = (boundary - nextCenterFrequency) / (thisCenterFrequency - nextCenterFrequency);
		//filterParameter *= GetMagnitudeFactor(filterBand);
	}
	else if(boundary >= nextCenterFrequency && boundary < samplingRate)
	{
		filterParameter = 0.0f;
	}

	return filterParameter;
}

/* 
 * Compute the band-dependent magnitude factor for the given filter band (Eq. 3)
 * Used for internal computation only - not the be called directly
 */
float GetMagnitudeFactor(unsigned int filterBand)
{
	float magnitudeFactor = 0.0f;
	
	if(filterBand >= 1 && filterBand <= 14)
	{
		magnitudeFactor = 0.015;
	}
	else if(filterBand >= 15 && filterBand <= 48) //48////////////////////////////////////////////////////////////////////////////
	{
		magnitudeFactor = 2.0f / (GetCenterFrequency(filterBand + 1) - GetCenterFrequency(filterBand -1));
	}

	return magnitudeFactor;
}

/*
 * Compute the center frequency (fc) of the specified filter band (l) (Eq. 4)
 * This where the mel-frequency scaling occurs. Filters are specified so that their
 * center frequencies are equally spaced on the mel scale
 * Used for internal computation only - not the be called directly
 */
float GetCenterFrequency(unsigned int filterBand)
{
	float centerFrequency = 0.0f;
	float exponent;
    float tmp;
#if 1
	if(filterBand == 0)
	{
		centerFrequency = 0;
	}
	else if(filterBand >= 1 && filterBand <= 14)
	{
		centerFrequency = (200.0f * filterBand) / 3.0f;
	}
	else
	{
		exponent = filterBand - 14.0f;
		centerFrequency = pow(1.0711703, exponent);
		centerFrequency *= 1073.4;
	}
#else /* http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/ */
    float Mel_Max, Mel_Val=0;
    
    if(filterBand == 0)
	{
		centerFrequency = 0;
	}
    else  
    {
        
        /* 1125*ln(1+f/700) */
        tmp = (float)(FS/(2.0*700.0)) + 1.0; 
        Mel_Max = (float)(1125.0*log(tmp));
        //printf("Mel(8000hz): %f ", Mel_Max);

        /* Mel value */
        Mel_Val = filterBand*Mel_Max/(NUMFILTERBANK+1);
        /* 700(e^(m/1125)-1) */ 
        tmp = (float)(Mel_Val/1125.0f);
        centerFrequency = 700*(exp(tmp)-1);                
    }
    
#endif
       	
	return centerFrequency;
}


/* precalculation the filterbank */
void PreCalcFilterBank(float Array[NUMFILTERBANK][NUMBINHALF], float fNorm[NUMFILTERBANK], int numBin, int numBank)
{
    int i,j;
    for(i=0; i< numBank; i++)
    {
        for (j=0; j<numBin; j++)
        {
            Array[i][j]= GetFilterParameter(FS,2*numBin, j, i+1);     //2*512
            //printf(" %f ",Array[i][j]);
        }
    
        fNorm[i] = NormalizationFactor(NUMFILTERBANK, i);   
        //printf(" %f ",fNorm[i]); 
    }    
}

/* Mel-Frequency Cepstal Coefficient (MFCC) update 
   Inputs: 
       + *Sample                       : samples in time domain
       + Filter[NUMFILTERBANK][NUMBINHALF]: FilterBank
       
       Note: Macro "NUMFILTERBANK": number of filter bank, e.x. NUMFILTERBANK=64
             Macro "NUMBINHALF": half of number of bin frequency    
   Output:
       + *Mel: MFCC according to input sample 

   Using: "PreCalcFilterBank" function should be called in advance to
   update date Filter[NUMFILTERBANK][NUMBINHALF]
*/
//void MFCC(float *x, float *y, float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm, float *Mel)
void MFCC(complex *Sample, float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm,float *Mel)
{
    float spectrum[NUMBINHALF];

    int i;
   
    //tmpBuff = (complex *)malloc(sizeof(complex)*2*NUMBINHALF);

    /* FFT calculation for input samples */
    fft( Sample, 2*NUMBINHALF, tmpBuff);

    //_FFT(1,10,x,y);
    

    //for (i=0; i < 20; i++)
    //{
    //    printf("%5.2f : %5.2f \n", Sample[i].Re,Sample[i].Im);
    //}
    /* Get Magnitude of FFT bin */
    abs_complex(Sample,spectrum,NUMBINHALF);

    //for (i=20; i < 30; i++)
    //{
    //    printf("%5.2f  \n", spectrum[i]);
    //}
    //for (i=50; i < 70; i++)
    //{
    //    printf("%5.2f  \n", spectrum[i]);
    //}
    // Compute coefficients
    GetMFCC(spectrum, Filter,fNorm, Mel); 
  

    //free(tmpBuff); 
}


void GetMFCC(float* spectralData,float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm, float *out)
{
    float LogEnergy[NUMFILTERBANK];
	float outerSum = 0.0f;
	float innerSum = 0.0f;
	unsigned int k, l, m;

    //LogEnergy = (float *)malloc(sizeof(float)*NUMFILTERBANK);

	/* log filterbank energies */
	for(l = 0; l < NUMFILTERBANK; l++) 
	{
		LogEnergy[l] = 0.0f;
		for(k = 0; k < NUMBINHALF ; k++)  
		{
			LogEnergy[l] += fabs(spectralData[k] *  Filter[l][k]);
		}

		if(LogEnergy[l] > 0.0f)
		{
			LogEnergy[l] = log(LogEnergy[l]); // The log of 0 is undefined, so don't use it
		}

         
    }

    /* Discrete Cosine Transform (DCT) */
    for (m=0; m < NUMFILTERBANK; m++)
    {
        out[m] = fNorm[m];
         
        outerSum = 0.0f; 
        for(l = 0; l <  NUMFILTERBANK; l++) 
        {            
		    innerSum = LogEnergy[l] * cos(((m * PI) / NUMFILTERBANK) * (l + 0.5f));

		    outerSum += innerSum;        
	    }        
        
        out[m] *= outerSum;
        //out[m] = LogEnergy[m];

        //printf("%5.2f  \n", result[m]); 
    } 

    //free(LogEnergy);
 	
}


/* 
   fft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute fft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute fft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = -sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
void fft( complex *v, int n, complex *tmp )
{
    if(n>1)  /* otherwise, do nothing and return */
    {			
          int k,m;    
          complex z, w;
          complex *vo, *ve;
        
          ve = tmp; 
          vo = tmp+n/2;
          for(k=0; k<n/2; k++) 
          {
            ve[k] = v[2*k];   //even sample
            vo[k] = v[2*k+1]; //odd sample
          }
          
          fft( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
          fft( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
          for(m=0; m<n/2; m++)
          {
              w.Re = cos(2*PI*m/(float)n);
              w.Im = -sin(2*PI*m/(float)n);
              z.Re = w.Re*vo[m].Re - w.Im*vo[m].Im;	/* Re(w*vo[m]) */
              z.Im = w.Re*vo[m].Im + w.Im*vo[m].Re;	/* Im(w*vo[m]) */
              v[  m  ].Re = ve[m].Re + z.Re;
              v[  m  ].Im = ve[m].Im + z.Im;
              v[m+n/2].Re = ve[m].Re - z.Re;
              v[m+n/2].Im = ve[m].Im - z.Im;
         }
    }
    return;
}

/* 
   ifft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute ifft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute ifft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
void ifft( complex *v, int n, complex *tmp )
{
  if(n>1) /* otherwise, do nothing and return */
  {			
    int k,m;    complex z, w, *vo, *ve;
    ve = tmp; vo = tmp+n/2;
    for(k=0; k<n/2; k++) {
      ve[k] = v[2*k];
      vo[k] = v[2*k+1];
    }
    ifft( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
    ifft( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
    for(m=0; m<n/2; m++) {
      w.Re = cos(2*PI*m/(float)n);
      w.Im = sin(2*PI*m/(float)n);
      z.Re = w.Re*vo[m].Re - w.Im*vo[m].Im;	/* Re(w*vo[m]) */
      z.Im = w.Re*vo[m].Im + w.Im*vo[m].Re;	/* Im(w*vo[m]) */
      v[  m  ].Re = ve[m].Re + z.Re;
      v[  m  ].Im = ve[m].Im + z.Im;
      v[m+n/2].Re = ve[m].Re - z.Re;
      v[m+n/2].Im = ve[m].Im - z.Im;
    }
  }
  return;
}


void abs_complex(complex * In, float *Out, int Size)
{
    int i;
    for (i=0; i < Size; i++)
    {
        Out[i] = sqrt(In[i].Re*In[i].Re +In[i].Im*In[i].Im)/(2*Size);    
    }
}

/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
short _FFT(short int dir,long m,float *x,float *y)
{
   long n,i,i1,j,k,i2,l,l1,l2;
   float c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++)
   {
      if (i < j) 
      {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j)
      {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++)
   {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++)
      {
         for (i=j;i<n;i+=l2)
         {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1)
   {
      for (i=0;i<n;i++) 
      {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return 1;
}

