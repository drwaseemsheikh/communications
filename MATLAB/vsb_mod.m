% MATLAB script to show Vestigial Sideband (VSB) Modulation and Demodulation

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
m = [ones(1,length(t1)) -2*ones(1,length(t2)) zeros(1,length(t3))]; % message signal
Bm = 100;                     % bandwidth of the baseband message signal

% plotting message signal
subplot(2,1,1);
plot(t,m);
grid on;
axis([0 2*T -2.5 1.5]);
xticks([0 0.025 0.05 0.1 0.15 0.2 0.25 0.3]);
title('Message or Source Signal');
xlabel('time (s)')
ylabel('Voltage(V)');

% plotting message signal spectrum
X = (1/fs) * fft(m);          % computes the Fast Fourier Transform (FFT)  
f = 0 : df : df*(N-1);        % frequency vector for the FFT
fnew = f - fs/2;              % shift the frequency vector to go from -fs/2 to fs/2
subplot(2,1,2);
plot(fnew, fftshift(abs(X))); % plots the magnitude spectrum of the message signal
grid on;
axis([-400 400 0 0.06]);
xticks([-400 -300 -200 -100 -80 -40 0 40 80 100 200 300 400]);
title('Message Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

%plotting DSB modulated signal
figure;
sd = m.*cos(2*pi*fc*t);      % double-sideband (DSB) modulated signal
subplot(2,1,1);
plot(t,sd);
grid on;
axis([0 2*T -2.5 2.5]);
title('DSB Modulated Signal');
xlabel('time (s)')
ylabel('Voltage(V)');

% plotting DSB modulated signal spectrum
SD = (1/fs) * fft(sd);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,2);
plot(fnew, fftshift(abs(SD))); % plots the magnitude spectrum of the modulated signal
grid on;
axis([-400 400 0 0.06]);
xticks([-400 -fc-Bm -fc -fc+Bm -100 0 100 fc-Bm fc fc+Bm 400]);
title('DSB Modulated Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% VSB shaping filter on the transmistter
f1 = 0: df : (fc-0.25*Bm)/df - 1;
f2 = (fc-0.25*Bm)/df : df : (fc+0.25*Bm)/df - 1;
f3 = (fc+0.25*Bm)/df : df : (fc+Bm)-1;
f4 = (fc+Bm) : df : fs/2;

Hir = [zeros(1, length(f1)) (1/Bm)*(f2-(fc-0.25*Bm))...
    ones(1,length(f3)) zeros(1, length(f4))];
Hil = fliplr(Hir);
Hi = [Hir Hil(2:length(Hir)-1)];
figure;
subplot(3,1,1);
plot(fnew, fftshift(Hi));
axis([-400 400 0 1.2]);
xtickangle(90);
xticks([-400 -fc-Bm -300 -fc-0.25*Bm -fc -fc+0.25*Bm -200 -fc+Bm -100 0 ...
    100 fc-Bm 200 fc-0.25*Bm fc fc+0.25*Bm 300 fc+Bm 400]);
xlabel('frequency (Hz)');
ylabel('Hi(f): VSB shaping filter');
grid on;

% magnitude spectrum of the VSB modulated signal at the output of the VSB filter
SV = SD.*Hi;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(3,1,2);
plot(fnew, fftshift(abs(SV))); % plots the magnitude spectrum of the LPF
axis([-400 400 0 0.06]);
xtickangle(90);
xticks([-400 -fc-Bm -300 -fc-0.25*Bm -fc -fc+0.25*Bm -200 -fc+Bm -100 0 ...
    100 fc-Bm 200 fc-0.25*Bm fc fc+0.25*Bm 300 fc+Bm 400]);
title('Magnitude Spectrum of the VSB modulated signal at the output of VSB filter');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');
grid on;

% VSB modulated signal at the output of the VSB filter
sv = fs * real(ifft(SV));
subplot(3,1,3);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(sv));
grid on;
axis([0 2*T -2.5 1.5]);
title('VSB Modulated Signal');
xlabel('time (s)')
ylabel('Voltage(V)');


% demodulation of the VSB signal using coherent detection and equalizer
% lowpass filter
% VSB signal after multiplication with a coherent local carrier
e = 4 * sv .* cos(2*pi*fc*t);
figure;
subplot(2,1,1);
plot(t,e);
grid on;
axis([0 2*T -2.5 2.5]);
title('Received VSB signal after coherent signal multiplication and before LPF');
xlabel('time (s)')
ylabel('Voltage(V)');

% spectrum of the VSB signal after multiplication with a coherent local
% carrier
E = (1/fs) * fft(e);
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,2);
plot(fnew, fftshift(abs(E)));
grid on;
axis([-800 800 0 0.06]);
xticks([-600 -500 -400 -200 -40 0 40 200 400 500 600 800]);
title('Spectrum of the received VSB signal after coherent signal multiplication and before LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% transfer function of the equalizer lowpass filter (Ho)
His = fftshift(Hi);
Hil = [His(fc/df+1:length(His)) zeros(1,fc/df)];
Hir = [zeros(1,fc/df) His(1:length(Hi)-fc/df)];
Hod = Hil+Hir;
Hod = [zeros(1,(fs/2-Bm)/df+1) Hod((fs/2-Bm)/df+2:(Bm+fs/2)/df) zeros(1,(fs/2-Bm)/df)];
Ho = [zeros(1,(fs/2-Bm)/df+1) 1./Hod((fs/2-Bm)/df+2:(Bm+fs/2)/df) zeros(1,(fs/2-Bm)/df)];

figure;
subplot(2,1,1);
plot(fnew, Hod);
xlabel('frequency (Hz)');
ylabel('Hi(f+fc) + Hi(f-fc)), |f|<=B');
xtickangle(90);
axis([-200 200 0 1.1]);
grid on;

subplot(2,1,2);
plot(fnew, Ho);
xlabel('frequency(Hz)');
ylabel('Ho(f), Equalizer Lowpass Filter')
xtickangle(90);
axis([-200 200 0 2.1]);
grid on;

% spectrum of the VSB demodulated signal at the output of the equalizer
% lowpass filter (Ho)
D = fftshift(Ho).*E;
f = 0 : df : df*(N-1);
fnew = f - fs/2;
figure;
subplot(2,1,1);
plot(fnew, fftshift(abs(D)));
grid on;
axis([-200 200 0 0.06]);
title('Magnitude Spectrum of the VSB demodulated signal at the output of Equalizer LPF');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% demodulated VSB signal at the output of the equalizer lowpass filter
d = fs * real(ifft(D));
subplot(2,1,2);
tnew = t-(ts*length(t)/2);
plot(tnew,fftshift(d));
grid on;
axis([0 2*T -2.5 1.5]);
title('Demodulated VSB Signal');
xlabel('time (s)')
ylabel('Voltage(V)');