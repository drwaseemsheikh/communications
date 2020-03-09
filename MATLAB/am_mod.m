% MATLAB script to show Conventional Amplitude Modulation (AM) and
% Demodulation using Carrier Recovery via Notch Filter and Envelope
% Detection

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
a = 1;                       % modulation index 0 <= a <= 1
A = 1;                       % carrier amplitude
fd = 0;                      % Doppler frequency shift in Hz
theta = pi/2;                % phase delay  

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

%plotting modulated received signal
figure;
s = A * (1 + a*x./abs(min(x))) .* cos(2*pi*(fc+fd)*t + theta);      % conventional AM signal
subplot(2,1,1);
plot(t,s);
grid on;
axis([0 2*T -2.5 2.5]);
title('Received Modulated Signal');
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
xticks([-400 -300 -(fc+fd) -200 -100 0 100 200 fc+fd 300 400]);
title('Received Modulated Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% Extracting the carrier from the received signal using a narrow-bandpass filter
% plotting the signal after LPF
B1 = 2;            % Bandwidth of the narrowband-bandpass filter. Must be even.
% DFT of an ideal narrowband-bandpass filter with bandwidth Bl Hz
Hnr = [zeros(1,(fc+fd-B1/2)/df) ones(1,B1/df+1) zeros(1,(fs/2-(fc+fd+B1/2))/df)]; 
Hnl = fliplr(Hnr);
Hn = [Hnr Hnl(2:length(Hnl)-1)];
figure;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,1);
plot(fnew, fftshift(abs(Hn))); % plots the magnitude spectrum of the narrow-bandpass filter
grid on;
axis([-400 400 0 1.2]);
xticks([-400 -300 -(fc+fd) -200 -100 0 100 200 fc+fd 300 400]);
title('Magnitude Spectrum of the Narrowband Bandpass Filter');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

%plotting the recovered carrier signal
RC = Hn.*S;
rc = fs * real(ifft(RC));
subplot(2,1,2);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(rc));
grid on;
axis([-0.2 0.2 -1.5 1.5]);
title('Recovered Carrier Signal from the Notch Filter');
xlabel('time (s)')
ylabel('Recovered Carrier');

% demodulation of the AM signal using recovered carrier signal
% plotting the input signal to the LPF
v = 2 * s .* rc;          % input signal to the the lowpass filter (LPF)
% to see the effect of phase and frequency error for coherent detection,
% comment the above line and uncomment the below line
%v = 2 * s .* cos(2*pi*fc*t);

figure;
subplot(2,1,1);
plot(t,v);
grid on;
axis([0 2*T -1 4]);
title('Received Signal after Recovered Carrier Multiplication and before LPF');
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
xticks([-800 -600 -2*(fc+fd) -400 -200 0 200 400 2*(fc+fd) 600 800]);
title('Spectrum of the Received Signal after Recovered Carrier Multiplication and before LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% plotting the signal after LPF
B2= 200;            % bandwidth of the lowpass filter
H = [ones(1,B2/df+1) zeros(1,N-2*B2/df-1) ones(1,B2/df)]; % DFT of an ideal lowpass filter with bandwidth 200 Hz
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
d = ((d-1)/a)*abs(min(x));    % rescale the message signal
figure;
subplot(2,1,1);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(d));
grid on;
axis([0 2*T -2.5 1.5]);
title('Demodulated Signal using Carrier Recovery');
xlabel('time (s)')
ylabel('Voltage(V)');

% Envelope Detection
env = abs(hilbert(s));
env = ((env-1)/a)*abs(min(x));    % rescale the message signal
subplot(2,1,2);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(env));
grid on;
axis([0 2*T -2.5 1.5]);
title('Demodulated Signal using Envelope Detection');
xlabel('time (s)')
ylabel('Voltage(V)');

% calculating power and efficiency of the AM signal
disp('power in the message signal')
px = sum(x.^2)/length(x)

disp('power in the AM modulated signal')
ps = sum(s.^2)/length(s)

disp('AM power efficiency in percent')
((1/2*a^2*px)/ps) * 100


