gcc -o alsa_record_MFCC alsa_record_MFCC.c libmfcc.c bmp.c -lm -lasound
sudo ./alsa_record_MFCC hw:2,0
