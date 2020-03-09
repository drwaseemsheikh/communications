% MATLAB script to show Double-Sideband (DSB) Modulation and Demodulation

clc
close all
fs = 4000;  % sampling frequency
ts = 1/fs;  % time sampling period
df = 1;     % frequency sampling period
N = fs/df;  % number of DFT smaples in the spectrum
T = 15e-2;  % time interval for plotting time domain signal 
fc = 250;   % carrier frequency in Hz
N1 = 100;   % number of samples in the positive pulse
N2 = 100;   % number of samples in the negative pulse
N3 = fs-N1-N2;             % number of samples in the 0 level of the pulse
t1 = 0:ts:(N1-1)*ts;         % time samples for the positive pulse
t2 = N1*ts:ts:(N1+N2-1)*ts;  % time samples for the negative pulse
t3 = (N1+N2)*ts:ts:(N-1)*ts; % time samples for the 0 level of the pulse
t = [t1 t2 t3];              % complete time vector
x = [ones(1,length(t1)) -2*ones(1,length(t2)) zeros(1,length(t3))]; % message signal

%plotting message signal
subplot(2,1,1);
plot(t,x);
grid on;
axis([0 2*T -2.5 1.5]);
title('Message or Source Signal');
xlabel('time (s)')
ylabel('Voltage(V)');

%plotting message signal spectrum
X = (1/fs) * fft(x);          % computes the Fast Fourier Transform (FFT)  
f = 0 : df : df*(N-1);        % frequency vector for the FFT
fnew = f - fs/2;              % shift the frequency vector to go from -fs/2 to fs/2
subplot(2,1,2);
plot(fnew, fftshift(abs(X))); % plots the magnitude spectrum of the message signal
grid on;
axis([-100 100 0 0.06]);
title('Message Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

%plotting modulated signal
figure;
s = x .* cos(2*pi*fc*t);      % double-sideband (DSB) modulated signal
subplot(2,1,1);
plot(t,s);
grid on;
axis([0 2*T -2.5 2.5]);
title('Modulated Signal');
xlabel('time (s)')
ylabel('Voltage(V)');

% plotting modulated signal spectrum
S = (1/fs) * fft(s);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,2);
plot(fnew, fftshift(abs(S))); % plots the magnitude spectrum of the modulated signal
grid on;
axis([-400 400 0 0.06]);
xticks([-400 -300 -fc -200 -100 0 100 200 fc 300 400]);
title('Modulated Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% demodulation of the DSB signal using coherent detection
%plotting the input signal to the LPF
v = 2 * s .* cos(2*pi*fc*t);          % input signal to the the lowpass filter (LPF)
figure;
subplot(2,1,1);
plot(t,v);
grid on;
axis([0 2*T -4.5 2.5]);
title('Received signal after coherent signal multiplication and before LPF');
xlabel('time (s)')
ylabel('Voltage(V)');

% plotting the spectrum of the received signal after coherent signal multiplication and before LPF
V = (1/fs) * fft(v);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,2);
plot(fnew, fftshift(abs(V))); % plots the magnitude spectrum of the modulated signal
grid on;
axis([-800 800 0 0.06]);
xticks([-800 -600 -2*fc -400 -200 0 200 400 2*fc 600 800]);
title('Spectrum of the received signal after coherent signal multiplication and before LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% plotting the signal after LPF
Bl = 200;            % bandwidth of the lowpass filter
H = [ones(1,Bl/df+1) zeros(1,N-2*Bl/df-1) ones(1,Bl/df)]; % DFT of an ideal lowpass filter with bandwidth 200 Hz
figure;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(3,1,1);
plot(fnew, fftshift(abs(H))); % plots the magnitude spectrum of the LPF
grid on;
axis([-800 800 0 1.2]);
title('Magnitude Spectrum of the LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

%plotting the impulse response of the LPF
h = fs * real(ifft(H));
subplot(3,1,2);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(h));
grid on;
axis([-0.1 0.1 -200 500]);
title('Impulse response of the LPF');
xlabel('time (s)')
ylabel('Impulse Response');

% plotting the spectrum of the demodulated sigal i.e., the output of the LPF
D = H.*V;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(3,1,3);
plot(fnew, fftshift(abs(D))); % plots the magnitude spectrum of the LPF
grid on;
axis([- 1000 1000 0 0.06]);
title('Magnitude Spectrum of the demodulated signal at the output of LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

%plotting the demodulated signal
d = fs * real(ifft(D));
figure;
subplot(2,1,1);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(d));
grid on;
axis([0 2*T -2.5 1.5]);
title('Demodulated Signal');
xlabel('time (s)')
ylabel('Voltage(V)');
