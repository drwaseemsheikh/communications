% MATLAB script to show analog Quadrature Amplitude Modulation (QAM) 
% Modulation and Demodulation, also known as Quadrature Multiplexing

clc
close all
fs = 4000;  % sampling frequency
ts = 1/fs;  % time sampling period
df = 1;     % frequency sampling period
N = fs/df;  % number of DFT smaples in the spectrum
T = 15e-2;  % time interval for plotting time domain signal 
fc = 250;   % carrier frequency in Hz
pw = 0.05;  % rectangular and triangular puse width in (s)
N1 = N*pw;   % number of samples in the positive pulse
N2 = N*pw;   % number of samples in the negative pulse
N3 = fs-N1-N2;             % number of samples in the 0 level of the pulse
t1 = 0:ts:(N1-1)*ts;         % time samples for the positive pulse
t2 = N1*ts:ts:(N1+N2-1)*ts;  % time samples for the negative pulse
t3 = (N1+N2)*ts:ts:(N-1)*ts; % time samples for the 0 level of the pulse
t = [t1 t2 t3];              % complete time vector
m1 = [rectangularPulse(0,0.05,t1) -2*rectangularPulse(0.05+ts,0.1,t2) zeros(1,length(t3))]; % message signal 1
m2 = [triangularPulse(0,0.05,t1) -2*triangularPulse(0.05+ts,0.1,t2) zeros(1,length(t3))]; % message signal 2

% plotting message signal 1
subplot(3,2,1);
plot(t,m1);
grid on;
axis([0 2*T -2.5 1.5]);
xlabel('time (s)')
ylabel('message signal 1');

% plotting message signal 1 spectrum
X1 = (1/fs) * fft(m1);           % computes the Fast Fourier Transform (FFT)  
f = 0 : df : df*(N-1);           % frequency vector for the FFT
fnew = f - fs/2;                 % shift the frequency vector to go from -fs/2 to fs/2
subplot(3,2,3);
plot(fnew, fftshift(abs(X1)));   % plots the magnitude spectrum of the message signal
grid on;
axis([-100 100 0 0.12]);
xlabel('frequency (Hz)');
ylabel('magnitude spectrum (message 1)');

subplot(3,2,5);
plot(fnew, fftshift(angle(X1))); % plots the magnitude spectrum of the message signal
grid on;
axis([-100 100 -pi-0.2 pi+0.2]);
xlabel('frequency (Hz)');
ylabel('phase spectrum (message 1)');

% plotting message signal 2
subplot(3,2,2);
plot(t,m2);
grid on;
axis([0 2*T -2.5 1.5]);
xlabel('time (s)')
ylabel('message signal 2');

% plotting message signal 2 spectrum
X2 = (1/fs) * fft(m2);         % computes the Fast Fourier Transform (FFT)  
f = 0 : df : df*(N-1);         % frequency vector for the FFT
fnew = f - fs/2;               % shift the frequency vector to go from -fs/2 to fs/2
subplot(3,2,4);
plot(fnew, fftshift(abs(X2))); % plots the magnitude spectrum of the message signal
grid on;
axis([-100 100 0 0.12]);
xlabel('frequency (Hz)');
ylabel('magnitude spectrum (message 2)');

subplot(3,2,6);
plot(fnew, fftshift(angle(X2))); % plots the magnitude spectrum of the message signal
grid on;
axis([-100 100 -pi-0.2 pi+0.2]);
xlabel('frequency (Hz)');
ylabel('phase spectrum (message 2)');

% plotting QAM signal and its spectrum
figure;
subplot(3,1,1);
sq = m1.*cos(2*pi*fc*t) + m2.*sin(2*pi*fc*t);
plot(t,sq);
grid on;
axis([0 2*T -4 4]);
xlabel('time (s)')
ylabel('QAM signal');

% plotting QAM signal spectrum
SQ = (1/fs) * fft(sq);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(3,1,2);
plot(fnew, fftshift(abs(SQ)));  % plots the magnitude spectrum of the modulated signal
grid on;
axis([-400 400 0 0.12]);
xticks([-400 -300 -fc -200 -100 0 100 200 fc 300 400]);
xlabel('frequency (Hz)');
ylabel('magnitude spectrum (QAM)');

subplot(3,1,3);
plot(fnew, fftshift(angle(SQ))); % plots the magnitude spectrum of the modulated signal
grid on;
axis([-400 400 -pi-0.2 pi+0.2]);
xticks([-400 -300 -fc -200 -100 0 100 200 fc 300 400]);
xlabel('frequency (Hz)');
ylabel('phase spectrum (QAM)');

% demodulation of the QAM signal using coherent detection
% plotting the in-phase input signal to the LPF
x1 = 2*sq.*cos(2*pi*fc*t);        % input signal to the the lowpass filter (LPF)
figure;
subplot(2,2,1);
plot(t,x1);
grid on;
axis([0 T -5 3]);
title('Received in-phase (I) signal after coherent signal multiplication and before LPF');
xlabel('time (s)')
ylabel('voltage (V)');

% plotting the spectrum of the received in-phase signal after coherent signal multiplication and before LPF
X1 = (1/fs) * fft(x1);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,2,3);
plot(fnew, fftshift(abs(X1)));
grid on;
axis([-800 800 0 0.12]);
title({'Spectrum of the in-phase (I) signal after coherent signal multiplication','and before LPF'});
xlabel('frequency (Hz)');
ylabel('magnitude spectrum');

% plotting the quadrature input signal to the LPF
x2 = 2*sq.*sin(2*pi*fc*t);        % input signal to the the lowpass filter (LPF)
subplot(2,2,2);
plot(t,x2);
grid on;
axis([0 T -5 3]);
title('Received quadrature (Q) signal after coherent signal multiplication and before LPF');
xlabel('time (s)')
ylabel('voltage (V)');

% plotting the spectrum of the received quadrature signal after coherent signal multiplication and before LPF
X2 = (1/fs) * fft(x2);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,2,4);
plot(fnew, fftshift(abs(X2)));
grid on;
axis([-800 800 0 0.12]);
title({'Spectrum of the quadrature (Q) signal after coherent signal multiplication','and before LPF'});
xlabel('frequency (Hz)');
ylabel('magnitude spectrum');

% plotting the in-phase signal after LPF
B = 200;            % bandwidth of the lowpass filter
H = [ones(1,B/df+1) zeros(1,N-2*B/df-1) ones(1,B/df)]; % DFT of an ideal lowpass filter with bandwidth 200 Hz
figure;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,1);
plot(fnew, fftshift(abs(H))); % plots the magnitude spectrum of the LPF
grid on;
axis([-800 800 0 1.2]);
title('Magnitude spectrum of the LPF');
xlabel('frequency (Hz)');
ylabel('magnitude spectrum');

%plotting the impulse response of the LPF
h = fs * real(ifft(H));
subplot(2,1,2);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(h));
grid on;
axis([-0.1 0.1 -200 500]);
title('Impulse response of the LPF');
xlabel('time (s)')
ylabel('impulse response');

% plotting the spectrum of the demodulated sigal i.e., the output of the LPF
figure;
D1 = H.*X1;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,2,1);
plot(fnew, fftshift(abs(D1))); % plots the magnitude spectrum of the LPF
grid on;
axis([-100 100 0 0.12]);
title('Magnitude spectrum of the demodulated in-phase signal at the output of LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude spectrum');

%plotting the in-phase demodulated signal
d1 = fs * real(ifft(D1));
subplot(2,2,3);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(d1));
grid on;
axis([0 2*T -2.5 1.5]);
xlabel('time (s)')
ylabel('demodulated I signal');

D2 = H.*X2;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,2,2);
plot(fnew, fftshift(abs(D2))); % plots the magnitude spectrum of the LPF
grid on;
axis([-100 100 0 0.12]);
title('Magnitude spectrum of the demodulated quadrature signal at the output of LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude spectrum');

%plotting the in-phase demodulated signal
d2 = fs * real(ifft(D2));
subplot(2,2,4);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(d2));
grid on;
axis([0 2*T -2.5 1.5]);
xlabel('time (s)')
ylabel('demodulated Q signal');



