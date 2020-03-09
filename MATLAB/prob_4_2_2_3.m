% MATLAB script to plot signals and their spectra for textbook problem 4.2-2(iii)

clc
close all
fs = 5000;   % sampling frequency
ts = 1/fs;   % time sampling period
df = 1;      % frequency sampling period
N = fs/df;   % number of DFT smaples in the spectrum and the number of samples in the time domain signal
             % N should be even
fc = 1000;   % carrier frequency in Hz
t = -N/2*ts:ts:(N/2-1)*ts;  % time vector with N smaples
m = cos(200*pi*t) + rectangularPulse(-1/200,1/200, t);  % message signal

%plotting message signal
subplot(3,1,1);
plot(t,m);
grid on;
axis([-0.1 0.1 -1.5 2.5]);
title('Message or Source Signal');
xlabel('time (s)')
ylabel('Voltage (V)');

%plotting the message signal magnitude spectrum
X = (1/fs) * fft(m);          % computes the Fast Fourier Transform (FFT)  
f = 0 : df : df*(N-1);        % frequency vector for the FFT
fnew = f - fs/2;              % shift the frequency vector to go from -fs/2 to fs/2
subplot(3,1,2);
plot(fnew, fftshift(abs(X))); % plots the magnitude spectrum of the message signal
grid on;
axis([-500 500 0 0.015]);
title('Message Signal Magnitude Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

%plotting the message signal phase spectrum
subplot(3,1,3);
plot(fnew, fftshift(angle(X))); % plots the phase spectrum of the message signal
grid on;
axis([-500 500 -pi-0.2 pi+0.2]);
title('Message Signal Phase Spectrum');
xlabel('frequency (Hz)');
ylabel('Phase Spectrum');

%plotting modulated signal
s = 2*m .* cos(2*pi*fc*t);      % double-sideband (DSB) modulated signal
figure;
subplot(3,1,1);
plot(t,s);
grid on;
axis([-0.1 0.1 -4 5]);
title('Modulated Signal');
xlabel('time (s)')
ylabel('Voltage (V)');

% plotting modulated signal spectrum
S = (1/fs) * fft(s);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(3,1,2);
plot(fnew, fftshift(abs(S))); % plots the magnitude spectrum of the modulated signal
grid on;
axis([-1200 1200 0 0.015]);
xticks([-1200 -1100 -1000 -900 -800 -500 -100 0 100 500 800 900 1000 1100 1200]);
title('Modulated Signal Magnitude Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

%plotting the modulated signal phase spectrum
subplot(3,1,3);
plot(fnew, fftshift(angle(S))); % plots the phase spectrum of the message signal
grid on;
axis([-500 500 -pi-0.2 pi+0.2]);
title('Modulated Signal Phase Spectrum');
xlabel('frequency (Hz)');
ylabel('Phase Spectrum');

