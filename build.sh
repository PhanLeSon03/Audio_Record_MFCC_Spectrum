gcc -o alsa_record_MFCC alsa_record_MFCC.c libmfcc.c -lm -lasound
sudo ./alsa-record-example hw:2,0
