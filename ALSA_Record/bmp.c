#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmp.h"
#include <sys/stat.h>
#include <fcntl.h>

int bmp_generator(char *filename, int width, int height, unsigned char *data)
{
    BITMAPFILEHEADER bmp_head;
    BITMAPINFOHEADER bmp_info;
    int size = width * height * 3;

    bmp_head.bfType = 0x4D42; // 'BM'
    bmp_head.bfSize= size + sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER); // 24 + head + info no quad    
    bmp_head.bfReserved1 = bmp_head.bfReserved2 = 0;
    bmp_head.bfOffBits = bmp_head.bfSize - size;
    // finish the initial of head

    bmp_info.biSize = 40;
    bmp_info.biWidth = width;
    bmp_info.biHeight = height;
    bmp_info.biPlanes = 1;
    bmp_info.biBitCount = 24; // bit(s) per pixel, 24 is true color
    bmp_info.biCompress = 0;
    bmp_info.biSizeImage = size;
    bmp_info.biXPelsPerMeter = 0;
    bmp_info.biYPelsPerMeter = 0;
    bmp_info.biClrUsed = 256 ;//256
    bmp_info.biClrImportant = 0;
    // finish the initial of infohead;

    // copy the data
    FILE *fp;
    if (!(fp = fopen(filename,"wb"))) return 0;
    fwrite(&bmp_head, 1, sizeof(BITMAPFILEHEADER), fp);
    fwrite(&bmp_info, 1, sizeof(BITMAPINFOHEADER), fp);
    fwrite(data, 1, size, fp);
    fclose(fp);

    return 1;
}

void bmpread(char *file_name)
{
    BITMAPFILEHEADER bmp_head;
    BITMAPINFOHEADER bmp_info; 
    int fd;
    if (!file_name)
        errx(1, "Filename not specified");

    if ((fd = open(file_name, O_RDONLY)) < 1)
        errx(1, "Error opening file");

    //if (!bmp_header)
    //    bmp_head = (bmp_head*)malloc(sizeof(bmp_head));

    //if (!bmp_info)
    //    bmp_info = (bmp_info*)malloc(sizeof(bmp_info));

    if (read(fd, bmp_head, sizeof(bmp_head)) < sizeof(bmp_head))
    {
        errx(1, "File broken: header");
    }
    else
    {
         printf("bmp_head.bfType: %d\n",bmp_head.bfType);
         printf("bmp_head.bfSize: %d\n",bmp_head.bfSize);  
         printf("bmp_head.bfReserved1: %d\n",bmp_head.bfReserved1);
         printf("bmp_head.bfReserved2: %d\n",bmp_head.bfReserved2);  
         printf("bmp_head.bfOffBits: %d\n",bmp_head.bfOffBits);
    }

    if (read(fd, bmp_info, sizeof(bmp_info)) < sizeof(bmp_info))
    {
        errx(1, "File broken: header");
    } 
    else
    {
         printf("bmp_info.biSize: %d\n",bmp_info.biSize);
         printf("bmp_info.biWidth: %d\n",bmp_info.biWidth);  
         printf("bmp_info.biHeight: %d\n",bmp_info.biHeight); 
         printf("bmp_info.biPlanes: %d\n",bmp_info.biPlanes);
         printf("bmp_info.biBitCount: %d\n",bmp_info.biBitCount);  
         printf("bmp_info.biCompress: %d\n",bmp_info.biCompress);
         printf("bmp_info.biSizeImage: %d\n",bmp_info.biSizeImage);
         printf("bmp_info.biXPelsPerMeter: %d\n",bmp_info.biXPelsPerMeter);  
         printf("bmp_info.biYPelsPerMeter: %d\n",bmp_info.biYPelsPerMeter);
         printf("bmp_info.biClrUsed: %d\n",bmp_info.biClrUsed);  
         printf("bmp_info.biClrImportant : %d\n",bmp_info.biClrImportant );
    }
  

    //if (*samples) free(*samples);
    //    *samples = (BYTE*)malloc(bmp_header.bfSize);
    
    //if (!*samples)
    //    errx(1, "Error allocating memory");
    
    //if (read(fd, *samples, bmp_header.bfSize) < bmp_header.size)
    //    errx(1, "File broken: samples");
    close(fd);
    //free(bmp_head);
    //free(bmp_info);
    //free(*samples);
}
