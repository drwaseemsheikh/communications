% MATLAB script to show FM and PM modulation and demodulation using a stair-case message
% function

clc
close all
fs = 40000;  % sampling frequency
ts = 1/fs;   % time sampling period
df = 1;      % frequency sampling period
N = fs/df;   % number of DFT smaples in the spectrum
T = 0.2;     % time interval for plotting time domain signal
T1 = 0.05;   % shorter time interval than T to show FM modulation
fc = 2000;   % carrier frequency in Hz
A = 1;       % carrier amplitdue
kf = 3000;  % frequency deviation constant ((rad/sec)/unit of message signal)
kp = 100;    % phase deviation constant (rad/unit of message signal)
N1 = 500;    % number of samples in the positive pulse
N2 = 500;    % number of samples in the negative pulse
N3 = fs-N1-N2;                 % number of samples in the 0 level of the pulse
t1 = 0:ts:(N1-1)*ts;           % time samples for the positive pulse
t2 = N1*ts:ts:(N1+N2-1)*ts;    % time samples for the negative pulse
t3 = 3*(N1+N2)*ts:ts:(N-1)*ts; % time samples for the 0 level of the pulse
t4 = [t1 t2];
t5 = t4(length(t4)) + [t1 t2];
t6 = t5(length(t5)) + [t1 t2];
t = [t4 t5 t6 t3];
m = [ones(1,length(t1)) -2*ones(1,length(t2)) ...
     0*ones(1,length(t1)) -1*ones(1,length(t2)) ...
     -0.5*ones(1,length(t1)) 0.75*ones(1,length(t2)) ...
     zeros(1,length(t3))]; % message signal

%plotting message signal
figure;
subplot(2,1,1);
plot(t,m);
grid on;
axis([0 T -2.2 1.2]);
title('Message or Source Signal');
xlabel('time (s)');
ylabel('Voltage(V)');

%plotting message signal spectrum
M = (1/fs) * fft(m);          % computes the Fast Fourier Transform (FFT)  
f = 0 : df : df*(N-1);        % frequency vector for the FFT
fnew = f - fs/2;              % shift the frequency vector to go from -fs/2 to fs/2
subplot(2,1,2);
plot(fnew, fftshift(abs(M))); % plots the magnitude spectrum of the message signal
grid on;
axis([-1000 1000 0 0.045]);
title('Message Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% computing power and 98% bandwidth of the FM signal
disp('Power in message signal using time domain (watt) = ');
pt = trapz(t,m.^2)

disp('Power in message signal spectrum using frequency domain (watt) =');
pf = trapz(f,abs(M).^2)

pfc = cumtrapz(fnew, abs(fftshift(M)).^2);
I = find(pfc >= 0.98*pf);
disp('98% bandwidth of message signal from simulation (Hz)');
BWm = fnew(I(1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FM Modulation and Demodulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting phase and frequency deviation for the FM signal
figure;
pdev1 = kf * cumtrapz(t,m);         % phase deviation for FM signal
fdev1 = (1/(2*pi)) * kf * m;        % frequency deviation for FM signal (Hz)
subplot(2,1,1);
plot(t,pdev1);
grid on;
axis([0 T -100 45]);
title('Phase Deviation of FM Signal (rad)');
xlabel('time (s)');
ylabel('phase deviation (rad)');

subplot(2,1,2);
plot(t,fdev1);
grid on;
axis([0 T -1050 550]);
title('Frequency Deviation of FM Signal (Hz)');
xlabel('time (s)');
ylabel('frequency deviation (Hz)');

%plotting FM modulated signal
figure;
s1 = A * cos(2*pi*fc*t + pdev1);         % FM modulated signal
subplot(2,1,1);
plot(t,s1);
grid on;
axis([0 T1 -1.1 1.1]);
title('FM Signal');
xlabel('time (s)');
ylabel('Voltage (V)');

% plotting FM signal spectrum
S1 = (1/fs) * fft(s1);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,2);
plot(fnew, fftshift(abs(S1))); % plots the magnitude spectrum of the modulated signal
grid on;
axis([-3500 3500 0 0.015]);
%xticks([-400 -300 -fc -200 -100 0 100 200 fc 300 400]);
title('FM Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% computing power and 98% bandwidth of the FM signal
disp('Power in the FM signal using time domain (watt) = ');
pt = trapz(t,s1.^2)

disp('Power in the FM signal spectrum using frequency domain (watt) =');
pf = trapz(f,abs(S1).^2)

disp('98% bandwidth of FM signal from simulation (Hz)');
BW_98_FM = powerbandwidth(f,fs,fc,df,S1,98)

disp('Peak frequency deviation of the FM signal (Hz)');
pfdfm = max(abs(fdev1))

disp('Deviation ratio of the FM signal');
drfm = pfdfm/BWm

disp('98% bandwidth of FM signal according to Carson''s rule (Hz)');
2 * BWm * (drfm + 1)

%plotting message signal, phase deviation, frequency deviation, and FM
%signal
figure;
subplot(4,1,1);
plot(t,m);
grid on;
axis([0 T1 -2.5 1.5]);
title('Message or Source Signal');
xlabel('time (s)');
ylabel('Voltage(V)');

subplot(4,1,2);
plot(t,pdev1);
grid on;
axis([0 T1 -100 40]);
title('Phase Deviation of FM Signal (rad)');
xlabel('time (s)');
ylabel('phase deviation (rad)');

subplot(4,1,3);
plot(t,fdev1);
grid on;
axis([0 T1 -1100 600]);
title('Frequency Deviation of FM Signal (Hz)');
xlabel('time (s)');
ylabel('frequency deviation (Hz)');

subplot(4,1,4);
plot(t,s1);
grid on;
axis([0 T1 -1.1 1.1]);
title('FM Signal');
xlabel('time (s)');
ylabel('Voltage (V)');

% plotting phase and frequency deviation for the PM signal
figure;
pdev2 = kp * m;                         % phase deviation for the PM signal (rad)
fdev2 = 1/(2*pi) * diff(pdev2)/ts;      % frequency deviation for the PM signal (Hz)
subplot(2,1,1);
plot(t,pdev2);
grid on;
axis([0 T -210 110]);
title('Phase Deviation of PM Signal (rad)');
xlabel('time (s)');
ylabel('phase deviation (rad)');

subplot(2,1,2);
plot(t(:,1:length(fdev2)),fdev2);
grid on;
axis([0 T -5000 5000]);
title('Frequency Deviation of PM Signal (Hz)')
xlabel('time (s)');
ylabel('frequency deviation (Hz)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FM Demodulation using Frequency Discriminator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1d = diff(s1)/ts;           % time difference (approx, differentation)
env = abs(hilbert(s1d));     % envelope = absolute value of positive analytic signal   
env = (env - 2*pi*fc)/kf;    % rescale envelope (wc + kf m(t))
figure;
subplot(2,1,1);
plot(t(:,1:length(env)),env);
grid on;
%axis([0 2*T -2.5 1.5]);
title('Demodulated Signal using Discriminator');
xlabel('time (s)')
ylabel('Voltage (V)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM Modulation and Demodulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plotting PM signal and its spectrum
figure;
s2 = A * cos(2*pi*fc*t + pdev2);    % PM modulated signal
subplot(2,1,1);
plot(t,s2);
grid on;
axis([0 T1 -1.1 1.1]);
title('PM Signal');
xlabel('time (s)');
ylabel('Voltage (V)');

% plotting PM signal spectrum
S2 = (1/fs) * fft(s2);          % Fast Fourier Transform of the modulated signal
f = 0 : df : df*(N-1);
fnew = f - fs/2;
subplot(2,1,2);
plot(fnew, fftshift(abs(S2))); % plots the magnitude spectrum of the modulated signal
grid on;
axis([-2300 2300 0 0.02]);
%xticks([-400 -300 -fc -200 -100 0 100 200 fc 300 400]);
title('PM Signal Spectrum');
xlabel('frequency (Hz)');
ylabel('Magnitude Spectrum');

% computing power and 98% bandwidth of the FM signal
disp('Power in the PM signal using time domain (watt) = ');
pt = trapz(t,s2.^2)

disp('Power in the PM signal spectrum using frequency domain (watt) =');
pf = trapz(f,abs(S2).^2)

disp('98% bandwidth of PM signal from simulation (Hz)');
BW_98_PM = powerbandwidth(f,fs,fc,df,S2,98)


disp('Peak frequency deviation of the PM signal (Hz)');
%pfdpm = max(abs(fdev2))   % since m(t) is not continuous, we can't apply this
%formula. Instead the frequency devaition of the PM signal is 0 since the
%information is modulated as phase shift and not as frequency change
pfdpm = 0

disp('Deviation ratio of the PM signal');
drpm = pfdpm/BWm

disp('98% bandwidth of PM signal according to Carson''s rule (Hz)');
2 * BWm * (drpm + 1)

%plotting message signal, phase deviation, frequency deviation, and PM
%signal
figure;
subplot(4,1,1);
plot(t,m);
grid on;
axis([0 T1 -2.5 1.5]);
title('Message or Source Signal');
xlabel('time (s)');
ylabel('Voltage(V)');

subplot(4,1,2);
plot(t,pdev2);
grid on;
axis([0 T1 -210 110]);
title('Phase Deviation of PM Signal (rad)');
xlabel('time (s)');
ylabel('phase deviation (rad)');

subplot(4,1,3);
plot(t(:,1:length(fdev2)),fdev2);
grid on;
axis([0 T1 -5000 5000]);
title('Frequency Deviation of PM Signal (Hz)');
xlabel('time (s)');
ylabel('frequency deviation (Hz)');

subplot(4,1,4);
plot(t,s2);
grid on;
axis([0 T1 -1.1 1.1]);
title('PM Signal');
xlabel('time (s)');
ylabel('Voltage (V)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM Demodulation using Frequency Discriminator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s2d = diff(s2)/ts;           % time difference (approx, differentation)
env = abs(hilbert(s2d));     % envelope = absolute value of positive analytic signal   
env = (env - 2*pi*fc)/kp;    % rescale envelope (wc + kp d/dtm(t))
env = env - mean(env);       % remove the DC component of the envelope
dmd = cumtrapz(t(:,1:length(env)),env(:,1:length(env)));  %integrate for PM demodulation
dmd = dmd - mean(dmd);       % remove the DC component
figure;
subplot(2,1,1);
plot(t(:,1:length(env)),env);
grid on;
%axis([0 2*T -2.5 1.5]);
title('Output of the Discriminator');
xlabel('time (s)')
ylabel('Voltage (V)');

subplot(2,1,2);
plot(t(:,1:length(dmd)),dmd);
grid on;
%axis([0 2*T -2.5 1.5]);
title('Demodulated PM Signal (Output of the Integrator)');
xlabel('time (s)')
ylabel('Voltage (V)');