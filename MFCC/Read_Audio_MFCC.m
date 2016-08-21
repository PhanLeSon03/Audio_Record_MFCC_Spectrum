clc;
clear;
NUMBINHALF = 512;
NUMFILTERBANK = 64;
FS = 16000;
T = 78; % number of frames
[FilterBank, fNorm, FIRCoef, Mel_Fre] =PreCalcFilterBank(NUMBINHALF, NUMFILTERBANK,FS);
Audio = audioread('record_2.wav');
V = zeros(2*NUMBINHALF,1);
Spectra = zeros(NUMFILTERBANK-1,2*T-1);
for i=0:2*T-2
    V = Audio(i*NUMBINHALF+1:(i+2)*NUMBINHALF);
    mfcc_result = MFCC(V, FilterBank,fNorm,FIRCoef);
    Spectra(:,i+1) = mfcc_result;
end
    
imagesc(Spectra)
set(gca,'YDir','normal');

imwrite(Spectra,'isabella_mfcc.bmp');
