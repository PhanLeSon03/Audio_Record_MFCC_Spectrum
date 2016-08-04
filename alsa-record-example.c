/* 
  A Minimal Capture Program

  This program opens an audio interface for capture, configures it for
  stereo, 16 bit, 44.1kHz, interleaved conventional read/write
  access. Then its reads a chunk of random data from it, and exits. It
  isn't meant to be a real program.

  From on Paul David's tutorial : http://equalarea.com/paul/alsa-audio.html

  Fixes rate and buffer problems

  sudo apt-get install libasound2-dev
  gcc -o alsa-record-example -lasound alsa-record-example.c && ./alsa-record-example hw:0
*/

#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <stdint.h>

typedef struct
{
    char RIFF_marker[4];
    uint32_t file_size;
    char filetype_header[4];
    char format_marker[4];
    uint32_t data_header_length;
    uint16_t format_type;
    uint16_t number_of_channels;
    uint32_t sample_rate;
    uint32_t bytes_per_second;
    uint16_t bytes_per_frame;
    uint16_t bits_per_sample;
} WaveHeader;



WaveHeader *genericWAVHeader(uint32_t sample_rate, uint16_t bit_depth, uint16_t channels);
int writeWAVHeader(FILE *fd, WaveHeader *hdr);

	      
main (int argc, char *argv[])
{
  int i;
  int err;
  char *buffer;
  char *buffer_mic1;
  int buffer_frames = 1024;
  unsigned int rate = 16000;
  char Channel = 8;
  snd_pcm_t *capture_handle;
  snd_pcm_hw_params_t *hw_params;
  snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;

  WaveHeader *fileAudio;
  char fileOut[] = "record.wav";
  FILE *out = fopen(fileOut,"wb");

  fileAudio = genericWAVHeader(16000, 16, 1);

  uint32_t pcm_data_size = 1024*2*100;
  fileAudio->file_size = pcm_data_size +36;
  writeWAVHeader(out,fileAudio); 
 

  if ((err = snd_pcm_open (&capture_handle, argv[1], SND_PCM_STREAM_CAPTURE, 0)) < 0) {
    fprintf (stderr, "cannot open audio device %s (%s)\n", 
             argv[1],
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "audio interface opened\n");
		   
  if ((err = snd_pcm_hw_params_malloc (&hw_params)) < 0) {
    fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n",
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "hw_params allocated\n");
				 
  if ((err = snd_pcm_hw_params_any (capture_handle, hw_params)) < 0) {
    fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n",
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "hw_params initialized\n");
	
  if ((err = snd_pcm_hw_params_set_access (capture_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
    fprintf (stderr, "cannot set access type (%s)\n",
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "hw_params access setted\n");
	
  if ((err = snd_pcm_hw_params_set_format (capture_handle, hw_params, format)) < 0) {
    fprintf (stderr, "cannot set sample format (%s)\n",
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "hw_params format setted\n");
	
  if ((err = snd_pcm_hw_params_set_rate_near (capture_handle, hw_params, &rate, 0)) < 0) {
    fprintf (stderr, "cannot set sample rate (%s)\n",
             snd_strerror (err));
    exit (1);
  }
	
  fprintf(stdout, "hw_params rate setted\n");

  if ((err = snd_pcm_hw_params_set_channels (capture_handle, hw_params, 8)) < 0) {
    fprintf (stderr, "cannot set channel count (%s)\n",
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "hw_params channels setted\n");
	
  if ((err = snd_pcm_hw_params (capture_handle, hw_params)) < 0) {
    fprintf (stderr, "cannot set parameters (%s)\n",
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "hw_params setted\n");
	
  snd_pcm_hw_params_free (hw_params);

  fprintf(stdout, "hw_params freed\n");
	
  if ((err = snd_pcm_prepare (capture_handle)) < 0) {
    fprintf (stderr, "cannot prepare audio interface for use (%s)\n",
             snd_strerror (err));
    exit (1);
  }

  fprintf(stdout, "audio interface prepared\n");

  buffer = malloc(buffer_frames * snd_pcm_format_width(format) / 8 * Channel);
  buffer_mic1 = malloc(buffer_frames * snd_pcm_format_width(format) / 8);
  fprintf(stdout, "buffer allocated\n");

  for (i = 0; i < 100; ++i) 
  {
      if ((err = snd_pcm_readi (capture_handle, buffer, buffer_frames)) != buffer_frames) 
      {
          //fprintf (stderr, "read from audio interface failed (%s)\n",
          //     err, snd_strerror (err));
          exit (1);
      }
      fprintf(stdout, "read %d done\n", i);
      uint16_t iSample;
      for (iSample=0; iSample < 2*buffer_frames*Channel; iSample++)
      {
          if (iSample%(2*Channel)==0)
          {  
              buffer_mic1[iSample/(Channel)]=buffer[iSample];
              buffer_mic1[iSample/(Channel)+1]=buffer[iSample+1];       
          }
      }
      fwrite(buffer_mic1,1,2*buffer_frames,out);
  }

  free(buffer);
  free(buffer_mic1);
 
  fclose(out);
  fprintf(stdout, "buffer freed\n");
	
  snd_pcm_close (capture_handle);
  fprintf(stdout, "audio interface closed\n");
  /////////////////////////////////////////////////////////////////

  
  //recordWAV(fileOut, fileAudio, 1);
  
  exit (0);
}

WaveHeader *genericWAVHeader(uint32_t sample_rate, uint16_t bit_depth, uint16_t channels)
{
    WaveHeader *hdr;
    hdr = malloc(sizeof(*hdr));
    if (!hdr) return NULL;

    memcpy(&hdr->RIFF_marker, "RIFF", 4);
    memcpy(&hdr->filetype_header, "WAVE", 4);
    memcpy(&hdr->format_marker, "fmt ", 4);
    hdr->data_header_length = 16;
    hdr->format_type = 1; //PCM
    hdr->number_of_channels = channels;
    hdr->sample_rate = sample_rate;
    hdr->bytes_per_second = sample_rate * channels * bit_depth / 8;
    hdr->bytes_per_frame = channels * bit_depth / 8;
    hdr->bits_per_sample = bit_depth;

    return hdr;
}


int writeWAVHeader(FILE *fd, WaveHeader *hdr)
{
    if (!hdr)
        return -1;

    fwrite(&hdr->RIFF_marker, 1, 4,fd); //RIFF
    fwrite(&hdr->file_size, 1,4, fd); //size
    fwrite(&hdr->filetype_header,1, 4,fd); //WAVE
    fwrite(&hdr->format_marker, 1,4, fd); //fmt
    fwrite(&hdr->data_header_length, 1, 4, fd);//length of fmt
    fwrite(&hdr->format_type, 1,2, fd); //PCM
    fwrite(&hdr->number_of_channels, 1,2, fd); //channel
    fwrite(&hdr->sample_rate, 1, 4, fd); 
    fwrite(&hdr->bytes_per_second, 1, 4, fd);
    fwrite(&hdr->bytes_per_frame, 1, 2, fd);
    fwrite(&hdr->bits_per_sample, 1, 2, fd);
    fwrite("data", 1,  4, fd);
    fwrite(&hdr->file_size, 1, 4, fd);
    //fwrite(&hdr->file_size-36, 1, 4, fd);
    //int count;
    //char * Dummy;
    //for (count = 44; count < 512; count++)
    //{
    //   *Dummy = 0x80;
    //   fwrite(Dummy, 1, 1, fd);
    //}

    return 0;
}





