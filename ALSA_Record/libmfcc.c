/*
 Phan Le Son
 plson03@gmail.com 
 */

#include <math.h>
#include "libmfcc.h"
#include <stdio.h>
#include <stdlib.h>
 
float FilterBank[NUMFILTERBANK-1][NUMBINHALF];
float fNorm[NUMFILTERBANK-1];
float result[NUMFILTERBANK];
float FIRCoef[2*NUMBINHALF];
complex tmpBuff[2*NUMBINHALF];


/* precalculation the filterband */
void PreCalcFilterBank(float Array[NUMFILTERBANK][NUMBINHALF], float fNorm[NUMFILTERBANK], int numBin, int numBand)
{
    float currVal;
    float * Mel_Fre;
    unsigned int i,j,iBand;

    Window(FIRCoef);

    Mel_Fre =(float *)malloc(sizeof(float)*(numBand+1));
    Mel_Fre[0] = 0;
    for (iBand=1; iBand <numBand+1; iBand++)
    {
        Mel_Fre[iBand] = GetCenterFrequency(iBand,FS, numBand +1);
    }
    
    for(i=0; i< numBand -1; i++)  //number of Triangle band
    {
        for (j=0; j<numBin; j++)
        {
            currVal  = (float)(j*FS/(2*NUMBINHALF));
            Array[i][j]= GetFilterParameter(currVal,i,Mel_Fre);  
            //printf(" %f ",Array[i][j]);
        }
    
        fNorm[i] = NormalizationFactor(NUMFILTERBANK, i);   
        //EnergyFac(Array, fNorm );
        //printf(" %f ",fNorm[i]); 
    }    

    free(Mel_Fre);
}


float GetFilterParameter(float currVal, unsigned int filterBand , float * Mel_Fre)
{
	float filterParameter = 0.0f;

	
	float Prev = Mel_Fre[filterBand];
	float Center = Mel_Fre[filterBand+1];
	float Next = Mel_Fre[filterBand+2];

    //printf(" %f - %f - %f - %f \n",boundary, prevCenterFrequency, thisCenterFrequency, nextCenterFrequency);
	if(currVal >= 0 && currVal < Prev)
	{
		filterParameter = 0.0f;
	}
	else if(currVal >= Prev && currVal < Center)
	{
		filterParameter = (currVal - Prev) / (Center - Prev);		
	}
	else if(currVal >= Center && currVal < Next)
	{
		filterParameter = (currVal - Center) / (Next - Center);		
	}
	else 
	{
		filterParameter = 0.0f;
	}

	return filterParameter;
}



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

/* Energy factors*/
void EnergyFac(float FilterBank[NUMFILTERBANK][NUMBINHALF], float fEnergyFac[NUMFILTERBANK])
{
    int i,j;
    for (i=0; i < NUMFILTERBANK; i++)
    {
        fEnergyFac[i] = 0.0;
        for (j=0; j < NUMBINHALF; j++) 
        {
            fEnergyFac[i] += FilterBank[i][j]*FilterBank[i][j];   
        }
        
        fEnergyFac[i] = 1.0/fEnergyFac[i];
    }

}


/* http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/ */
float GetCenterFrequency(unsigned int filterBand, int fs, int numPointBank)
{
	float centerFrequency = 0.0f;
    float tmp;
   
    float Mel_Max, Mel_Val=0;
    
    if (filterBand == 0)
	{
		centerFrequency = 0;
	}
    else  
    {
        
        /* 1125*ln(1+f/700) */
        tmp = (float)(fs/(2.0*700.0)) + 1.0; 
        Mel_Max = (float)(1125.0*log(tmp));
        //printf("Mel(8000hz): %f ", Mel_Max);

        /* Mel value */
        Mel_Val = filterBand*Mel_Max/(numPointBank);
        /* 700(e^(m/1125)-1) */ 
        tmp = (float)(Mel_Val/1125.0f);
        centerFrequency = 700*(exp(tmp)-1);                
    }    
       	
	return centerFrequency;
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
void MFCC(complex *Sample, float Filter[NUMFILTERBANK][NUMBINHALF],float *fNorm,float *Mel)
{
    float spectrum[NUMBINHALF];

    int i;

    /* Preemphasis filter*/
    //Preemphasis(Sample);
    
    //tmpBuff = (complex *)malloc(sizeof(complex)*2*NUMBINHALF);

    /* Hanning Window */
    for (i=0; i < 2*NUMBINHALF; i++)
    {
         Sample[i].Re = FIRCoef[i]*Sample[i].Re; 
    }
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
  
    //lifter(Mel,NUMFILTERBANK); 
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
	for(l = 0; l < NUMFILTERBANK-1; l++) 
	{
		LogEnergy[l] = 0.0f;
		for(k = 0; k < NUMBINHALF ; k++)  
		{
			LogEnergy[l] += fabs(spectralData[k] *  Filter[l][k]);
		}

        
		if(LogEnergy[l] < 0.001f) LogEnergy[l] = 0.01f;
		
		LogEnergy[l] = log(LogEnergy[l]); // The log of 0 is undefined, so don't use it
		
         
    }

    /* Discrete Cosine Transform (DCT) */
    for (m=0; m < NUMFILTERBANK - 1; m++)
    {
        out[m] = fNorm[m];
         
        outerSum = 0.0f; 
        for(l = 0; l <  NUMFILTERBANK-1; l++) 
        {            
		    innerSum = LogEnergy[l] * cos(((m * PI) / (NUMFILTERBANK-1)) * (l + 0.5f));

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

void lifter(float * Cepstra, unsigned char Len)
{
    int i=0;
    float Lift;
    /* L =22 */ 
    /* Lift = 1 + (L/2)*sin(PI*n/L) */
    for (i=0; i < Len; i++)
    {
        Lift = 1.0f + 32*sin(PI*i/64.0);
        Cepstra[i] = Lift*Cepstra[i];
    
    } 
    

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

/* Hanning window */
void Window (float *FIRCoef)
{
    int i;
    for(i=0;i<NUMBINHALF*2;i++)
    {
        FIRCoef[i] = (float) NUMBINHALF*2.0f;
        FIRCoef[i] *= 0.5f*(
                            1.0f - cos((2.0f*PI*i)/((float)NUMBINHALF*2.0f -1.0f))
                           );  
    }
}


void Preemphasis(complex *Data)
{
    int i;
    float fCoef = 0.95f;
    for(i=1; i < 2*NUMBINHALF ; i++)
    {
        Data[i].Re = Data[i].Re -0.95f*Data[i-1].Re; 
    }
}

















