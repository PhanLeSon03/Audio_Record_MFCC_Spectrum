gcc -o alsa-record-example alsa-record-example.c libmfcc.c -lm -lasound
sudo ./alsa-record-example hw:2,0
