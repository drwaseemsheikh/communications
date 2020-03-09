% Autocorrelation and power spectral density of polar NRZ baseband encoding
% scheme

close all;

Tb = 1/1000;                        % bit interval (sec)
Rb = 1/Tb;                          % bit rate (bits/sec)
A = 1;                              % amplitude of the wave
fs = 100 * Rb;                      % sampling frequency (samples/sec)
ts = 1/fs;                          % sampling interval (sec)
df = 1;                             % DFT frequency sampling interval
N = fs/df;                          % number of DFT samples (must be even)
t = 0 : ts : (N/2-1) * ts;          % time vector (N/2 samples because
                                    % autocorrelation vector has twice as many
                                    % samples
tp = 0 : ts : Tb-ts;                % time vector for a single pulse

p = rectangularPulse(0, Tb, tp);    % rectangular pulse
p(1) = 1;                           % change the first sample from 0.5 to 1
Np = length(t) / length(tp);        % number of pulses
x = p;
J = 30;                             % number of iterations for averaging the autocorrelation function
rx = zeros(1, 2*Np*length(x)-1);    % initialzie the autocorrelation vector
Sx = zeros(1, 2*Np*length(x)-1);    % initialize the PSD vector

for j = 1 : J
    x = p;
    for i = 2 : Np
        ai = unidrnd(2);
        
        if ai == 1
            ai = -1;
        else
            ai = 1;
        end
        
        pi = ai * p;
        x = [x pi];
    end
    
    rx = rx + xcorr(x, 'biased');
end

rx = rx/J;         % average of the autocorrelation function

figure;
subplot(3,1,1);
plot(t,x);
ylabel('signal');
xlabel('time (s)');
title('Autocorrelation and Power Spectrum of Polar Baseband Encoding');
grid on;
axis([0 0.5 -1.2 1.2])

subplot(3,1,2);
tau = -length(t)+1 : 1 : length(t)-1;
tau = ts * tau;
plot(tau, rx);
ylabel('autcorrelation funcntion');
xlabel('Shift (s)');
grid on;
axis([-0.003 0.003 0 1.1]);

subplot(3,1,3);
Sx = (1/fs) * fft(rx);                  % computes the Fast Fourier Transform (FFT)
f = 0 : df : 2 * df * (length(x) - 1);  % frequency vector for the FFT
fnew = f - fs/2;                        % shift the frequency vector to go from -fs/2 to fs/2
plot(fnew, fftshift(abs(Sx)));          % plots the magnitude spectrum of the message signal
xlabel('frequency (Hz)');
ylabel('power spectrum');
grid on;
axis([-2000 2000 0 1.5e-3])