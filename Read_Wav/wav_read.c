
#include <stdio.h>
#include <stdlib.h>
//#include <alsa/asoundlib.h>
#include <stdint.h>
#include "libmfcc.h"
#include "bmp.h"
#include "wav.h"
#include <string.h>
#define T 78 
RGB_data * BMP;
complex V[2*NUMBINHALF];
char FileOut[200];
char FileOutTxt[200];
extern float FilterBank[NUMFILTERBANK-1][NUMBINHALF];
extern float fNorm[NUMFILTERBANK-1];


int main (int argc, char **argv)
{
    float mfcc_result[NUMFILTERBANK-1];
    FILE *fptr;
   
    float Min=0.0f,Max=0.0f;
    int16_t *samples = NULL;
    uint16_t numFrame = 0;
    uint16_t i,j,k;
    int32_t Val_RGB;
    uint32_t idxCoeff;
    int tmp;
    char Len;
    
    Len = strlen(argv[1]);
    for ( i=0; i < Len-4; i++)
    {
        FileOut[i]= argv[1][i];
        FileOutTxt[i]= argv[1][i];
    }
    
    FileOut[Len-4] = '.';
    FileOut[Len-3] = 'b';
    FileOut[Len-2] = 'm';
    FileOut[Len-1] = 'p';
    FileOut[Len] = '\0';

    FileOutTxt[Len-4] = '.';
    FileOutTxt[Len-3] = 't';
    FileOutTxt[Len-2] = 'x';
    FileOutTxt[Len-1] = 't';
    FileOutTxt[Len] = '\0';
    
    /* Sotarage  MFCC coeffience */
    fptr=fopen(FileOutTxt,"w");
    
    PreCalcFilterBank(FilterBank,fNorm, NUMBINHALF, NUMFILTERBANK);
    //printf("Length: %d \n", strlen(argv[1]));
    wavread(argv[1], &samples);
    printf("No. of channels: %d\n", header->num_channels);
    printf("Sample rate: %d\n", header->sample_rate);
    printf("Bit rate: %dkbps\n", header->byte_rate*8 / 1000);
    printf("Bits per sample: %d\n\n", header->bps);
    //printf("Sample 0: %d\n", samples[0]);
    //printf("Sample 1: %d\n", samples[1]);
    // Modify the header values & samples before writing the new file
    wavwrite("track2.wav", samples);

    numFrame = header->datachunk_size/(2*2*NUMBINHALF) ;   
    //printf("Num Frame:  %d \n", numFrame );
    /*
    BMP = (RGB_data **)malloc((NUMFILTERBANK-1)*sizeof(RGB_data*));
    for (i=0; i < NUMFILTERBANK-1; i++)
    {
        BMP[i] = (RGB_data *)malloc(sizeof(RGB_data)*(2*numFrame)) ; 
        if(BMP[i] == NULL)
        {
            fprintf(stderr, "out of memory\n");
            exit (0);
        }        
    } 
    */
    BMP = (RGB_data *)malloc((NUMFILTERBANK-1)*sizeof(RGB_data*)*2*numFrame);     
    //memset(BMP, 0, sizeof(BMP));    
    for (i = 0; i < 2*numFrame-1; i++) 
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
           fprintf(fptr,"%f ", mfcc_result[idxCoeff]);
           //Val_RGB = (int)(255/(63.164356+46.877232)*(mfcc_result[idxCoeff]+10));
           Val_RGB = (int)(255*(mfcc_result[idxCoeff]));  
           //Val_RGB = tmp*tmp;
           //Val_RGB = Val_RGB/255;
           if (Val_RGB < 0)  Val_RGB = 0;
           if (Val_RGB > 255)  Val_RGB = 255;

           BMP[(NUMFILTERBANK - idxCoeff-2)*2*numFrame +i].g =  (BYTE)(Val_RGB); 
           BMP[(NUMFILTERBANK - idxCoeff-2)*2*numFrame +i].b =  (BYTE)(Val_RGB);
           BMP[(NUMFILTERBANK - idxCoeff-2)*2*numFrame +i].r =  (BYTE)(Val_RGB); 
           //BMP[NUMFILTERBANK - idxCoeff-1][i].b = (BYTE)((Val_RGB & 0x00FF00)>>8);
           //BMP[NUMFILTERBANK - idxCoeff-1][i].r = (BYTE)((Val_RGB & 0xFF0000)>>16); 
           //printf(" %d  ", Val_RGB); 
           //if (mfcc_result[idxCoeff] < Min) Min=mfcc_result[idxCoeff];
           //if (mfcc_result[idxCoeff] > Max) Max=mfcc_result[idxCoeff];  
       }
       fprintf(fptr,"\r\n");
    }

    //printf("Min: %f  Max: %f \n", Min, Max);

    bmp_generator(FileOut, 2*numFrame, NUMFILTERBANK -1 ,(BYTE*) (BMP)); 

    //bmpread("record_21.bmp");

    free(header);
    free(samples);
    free(BMP);
    fclose(fptr);
    exit (0);
}
