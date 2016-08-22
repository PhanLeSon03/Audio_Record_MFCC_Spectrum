

#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <stdint.h>
#include "libmfcc.h"
#include "bmp.h"
#include "wav.h"

#define T 78
RGB_data BMP[NUMFILTERBANK-1][2*T];
complex V[2*NUMBINHALF];

extern float FilterBank[NUMFILTERBANK-1][NUMBINHALF];
extern float fNorm[NUMFILTERBANK-1];


main()
{
    float mfcc_result[NUMFILTERBANK-1];
    
    float Min=0.0f,Max=0.0f;
    int16_t *samples = NULL;
    uint16_t numFrame = 0;
    uint16_t i,j,k;
    int32_t Val_RGB;
    uint32_t idxCoeff;
    int tmp;

    memset(BMP, 0, sizeof(BMP));
    PreCalcFilterBank(FilterBank,fNorm, NUMBINHALF, NUMFILTERBANK);

    wavread("record_2.wav", &samples);
    printf("No. of channels: %d\n", header->num_channels);
    printf("Sample rate: %d\n", header->sample_rate);
    printf("Bit rate: %dkbps\n", header->byte_rate*8 / 1000);
    printf("Bits per sample: %d\n\n", header->bps);
    //printf("Sample 0: %d\n", samples[0]);
    //printf("Sample 1: %d\n", samples[1]);
    // Modify the header values & samples before writing the new file
    wavwrite("track2.wav", samples);

    numFrame = sizeof(samples)/(2*NUMBINHALF) ;     
    for (i = 0; i < 2*T-1; i++) 
    {
       for (j=0; j< 2*NUMBINHALF;j++)
       {
           V[j].Re = (float)(samples[i*NUMBINHALF + j]); 
           //printf(" %d", samples[i*NUMBINHALF + j]);
           V[j].Im = (float)0.0f;
       } 
       //printf("\n");
       MFCC(V, FilterBank ,fNorm, mfcc_result);

       //printf("MFCC:");
       for(idxCoeff = 0; idxCoeff < NUMFILTERBANK -1; idxCoeff++)
       {	     
           //printf(" %f ", mfcc_result[idxCoeff]);
           //Val_RGB = (int)(255/(63.164356+46.877232)*(mfcc_result[idxCoeff]+10));
           Val_RGB = (int)(255*(mfcc_result[idxCoeff]));  
           //Val_RGB = tmp*tmp;
           //Val_RGB = Val_RGB/255;
           if (Val_RGB < 0)  Val_RGB = 0;
           if (Val_RGB > 255)  Val_RGB = 255;

           BMP[NUMFILTERBANK - idxCoeff-2][i].g =  (BYTE)(Val_RGB); 
           BMP[NUMFILTERBANK - idxCoeff-2][i].b =  (BYTE)(Val_RGB);
           BMP[NUMFILTERBANK - idxCoeff-2][i].r =  (BYTE)(Val_RGB); 
           //BMP[NUMFILTERBANK - idxCoeff-1][i].b = (BYTE)((Val_RGB & 0x00FF00)>>8);
           //BMP[NUMFILTERBANK - idxCoeff-1][i].r = (BYTE)((Val_RGB & 0xFF0000)>>16); 
           //printf(" %d  ", Val_RGB); 
           //if (mfcc_result[idxCoeff] < Min) Min=mfcc_result[idxCoeff];
           //if (mfcc_result[idxCoeff] > Max) Max=mfcc_result[idxCoeff];  
       }
       //printf("\n");
    }

    //printf("Min: %f  Max: %f \n", Min, Max);

    bmp_generator("./record_2.bmp", 2*T, NUMFILTERBANK -1 ,(BYTE*) (BMP)); 

    //bmpread("record_21.bmp");

    free(header);
    free(samples);
    exit (0);
}
