# Audio-Compression-using-Matlab-using-Haar-wavelets
This project consist of audio compression technology using matlab in which we use haar wavelets algorithm.

This code performs wavelet compression over a .wav file.The scheme is a simplified version of the one described on the paper "low bit rate" transparent audio compression using adapted wavelets.

NOTE: The file must be in the same folder where this file is located.If you want to try this scheme with other audio file, please change the name of the cariable "file".avoid using long audio files or long silences at the beginning of the file for computational constraints.


clear;clc;
file='Coltrane.wav';
wavelet='dB10';
level=5;
frame_size=2048;
psychoacoustic='on '; %if it is off it uses 8 bits/frame as default
wavelet_compression = 'on ';
heavy_compression='off';
compander='on ';
quantization ='on ';
%%%%%%%%%%%%%%%%%%%%%%%%%
% ENCODER %
%%%%%%%%%%%%%%%%%%%%%%%%%
[x,Fs,bits] = wavread(file);
xlen=length(x);
t=0:1/Fs:(length(x)-1)/Fs;
%decomposition using N equal frames
step=frame_size;
N=ceil(xlen/step);
%computational variables
Cchunks=0;
Lchunks=0;
Csize=0;
PERF0mean=0;
PERFL2mean=0;
n_avg=0;
n_max=0;
n_0=0;
n_vector=[];
for i=1:1:N
if (i==N)
frame=x([(step*(i-1)+1):length(x)]);
else
frame=x([(step*(i-1)+1):step*i]);
end
%wavelet decomposition of the frame
[C,L] = wavedec(frame,level,wavelet);
%wavelet compression scheme
if wavelet_compression=='on '
[thr,sorh,keepapp] = ddencmp('cmp','wv',frame);
if heavy_compression == 'on '
thr=thr*10^6;
end
[XC,CXC,LXC,PERF0,PERFL2] = wdencmp('gbl',C, L, wavelet,
level,thr,sorh,keepapp);
C=CXC;
L=LXC;
PERF0mean=PERF0mean + PERF0;
PERFL2mean=PERFL2mean+PERFL2;
end
%Psychoacoustic model
if psychoacoustic=='on '
P=10.*log10((abs(fft(frame,length(frame)))).^2);
Ptm=zeros(1,length(P));
%Inspect spectrum and find tones maskers
for k=1:1:length(P)
if ((k<=1) | (k>=250))
bool = 0;
elseif ((P(k)<P(k-1)) | (P(k)<P(k+1))),
bool = 0;
elseif ((k>2) & (k<63)),
bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)));
elseif ((k>=63) & (k<127)),
bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)) & (P(k)>(P(k-
3)+7)) & (P(k)>(P(k+3)+7)));
elseif ((k>=127) & (k<=256)),
bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)) & (P(k)>(P(k-
3)+7)) & (P(k)>(P(k+3)+7)) & (P(k)>(P(k-4)+7)) & (P(k)>(P(k+4)+7)) &
(P(k)>(P(k-5)+7)) & (P(k)>(P(k+5)+7)) & (P(k)>(P(k-6)+7)) &
(P(k)>(P(k+6)+7)));
else
bool = 0;
end
if bool==1
Ptm(k)=10*log10(10.^(0.1.*(P(k1)))+10.^(0.1.*(P(k)))+10.^(0.1.*P(k+1)));
end
end
sum_energy=0;%sum energy of the tone maskers
for k=1:1:length(Ptm)
sum_energy=10.^(0.1.*(Ptm(k)))+sum_energy;
end
E=10*log10(sum_energy/(length(Ptm)));
SNR=max(P)-E;
n=ceil(SNR/6.02);%number of bits required for quantization
if n<=3%to avoid distortion by error of my psychoacoustic model.
n=4;
n_0=n_0+1;
end
if n>n_max
n_max=n;
end
n_avg=n+n_avg;
n_vector=[n_vector n];
end
%Compander(compressor)
if compander=='on '
Mu=255;
C = compand(C,Mu,max(C),'mu/compressor');
end
%Quantization
if quantization=='on '
if psychoacoustic=='off'
n=8;%default number of bits for each frame - sounds better but
uses more bits
end
partition = [min(C):((max(C)-min(C))/2^n):max(C)];
codebook = [min(C):((max(C)-min(C))/2^n):max(C)];
[index,quant,distor] = quantiz(C,partition,codebook);
%find and correct offset
offset=0;
for j=1:1:N
if C(j)==0
offset=-quant(j);
break;
end
end
quant=quant+offset;
C=quant;
end
%Put together all the chunks
Cchunks=[Cchunks C]; %NOTE: if an error appears in this line just
modify the transpose of C
Lchunks=[Lchunks L'];
Csize=[Csize length(C)];
Encoder = round((i/N)*100) %indicator of progess
end
Cchunks=Cchunks(2:length(Cchunks));
Csize=[Csize(2) Csize(N+1)];
Lsize=length(L);
Lchunks=[Lchunks(2:Lsize+1) Lchunks((N-1)*Lsize+1:length(Lchunks))];
PERF0mean=PERF0mean/N %indicator
PERFL2mean=PERFL2mean/N %indicator
n_avg=n_avg/N%indicator
n_max%indicator
end_of_encoder='done'%indicator of progess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In this part the signal is stored with the new format
or transmitted by frames
This new format uses this parameters:
header: N, Lsize, Csize.
body: Lchunks (small), Cchunks(smaller signal because now it is
quantized with less bit and coded)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% DECODER %
%%%%%%%%%%%%%%%%%%%%%%%%%
%reconstruction using N equal frames of length step (except the last
one)
xdchunks=0;
for i=1:1:N
if i==N
Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(2)+(Csize(1)*(i-
1))]);
%Compander (expander)
if compander=='on '
if max(Cframe)==0
else
Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
end
end
xd = waverec(Cframe,Lchunks(Lsize+2:length(Lchunks)),wavelet);
else
Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(1)*i]);
%Compander (expander)
if compander=='on '
if max(Cframe)==0
else
Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
end
end
xd = waverec(Cframe,Lchunks(1:Lsize),wavelet);
end
xdchunks=[xdchunks xd];
Decoder = round((i/N)*100) %indicator of progess
end
xdchunks=xdchunks(2:length(xdchunks));
distorsion = sum((xdchunks-x').^2)/length(x)
end_of_decoder='done'
%creating audio files with compressed schemes
wavwrite(xdchunks,Fs,bits,'output.wav') %this does not represnet the
real compression achieved. It is only to hear the results
end_of_writing_file='done'%indicator of progess
