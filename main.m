clear all
close all
%% read wavs
[mix,fs] = audioread("data/1observation.wav");
[source,fs] = audioread('data/1source.wav');

%% stft
windowSize=32;
shiftCof=0.5;
fsResample = 16000;
fftSize = fsResample*windowSize/1000;   % Number of points for FFT
shiftSize = shiftCof*fftSize; % Number of points for window shift

%% set parameters
beta = 1; % shape parameter
p = 0.5; % domain parameter
refMic = 1; % reference channel
it = 40;  % iterations for updates
itNMF = 1;   % iterations for NMF parameters in each update
nb = 40; % the number of basis vectors in NMF
L = 17;  % prediction order
delta = 1; % prediction delay

%% perform dereverberation
seed=1; 
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));
[drb] = WPECGGNMF(mix(:,1:4), nb, L, delta, fftSize, shiftSize, it,itNMF, beta,p,refMic);
drb = drb(:,refMic);
audiowrite('data/dereverberatedSignal.wav',drb,fsResample);




