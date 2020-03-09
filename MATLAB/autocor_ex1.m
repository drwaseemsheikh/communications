% examples of autocorrelation and power spectral density

close all;

Ns = 10;                        % number of sinuosoids in the signal
theta = 2 * pi * rand(1, Ns);   % phase delay of the sinusoids (rad)
Dfmax = 10;                     % maximum change in the sinusoid frequency (Hz)
Df = Dfmax * rand(1,Ns);        % change in the sinuosid frequency
fmin = [1 1 40 40 100 100];     % minimum frequency of a sinusoid
Pn = [0 2 0 2 0 5];             % noise power

for j = 1:length(fmin)
    fmax = fmin(j) + Dfmax;     % maximum frequency of a sinusoid
    df = 0.1;                   % frequency sampling interval
    fs = 10 * fmax;             % sampling frequency
    N = fs/df;                  % number of DFT samples
    ts = 1/fs;                  % sampling interval
    t = 0 : ts : (N-1) * ts;    % time vecotr
    Am = 1;                     % maximum amplitude of the sinusoids
    A = Am * rand(1,Ns);        % amplitude of the sinuoids
    
    x = zeros(1, length(t));    % signal
    x = x + Pn(j) * randn(1, length(x));
    
    if Pn(j) ~= 0
        ttlfig = char("Autcorrelation and Power Spectrum of a Signal with Frequency Content between " + ...
            fmin(j) + " and " + fmax + " Hz" + " and White Noise");
    else
        ttlfig = char("Autcorrelation and Power Spectrum of a Signal with Frequency Content between " + ...
            fmin(j) + " and " + fmax + " Hz");
    end
    
    for i = 1 : Ns
        x = x + A(i) * cos(2 * pi * (fmin(j) + Df(i)) * t + theta(i));
    end
    
    rx = xcorr(x, 'biased');
    tau = -length(t)+1 : 1 : length(t)-1;
    tau = ts * tau;
    
    figure;
    subplot(3,1,1);
    plot(t,x);
    ylabel('signal');
    xlabel('time (s)');
    title(ttlfig);
    grid on;
    
    subplot(3,1,2);
    plot(tau, rx);
    ylabel('autcorrelation funcntion');
    xlabel('Shift (s)');
    grid on;
    
    subplot(3,1,3);
    Sx = (1/fs) * fft(rx);                  % computes the Fast Fourier Transform (FFT)
    f = 0 : df/2 : df * (length(x) - 1);      % frequency vector for the FFT
    fnew = f - fs/2;                        % shift the frequency vector to go from -fs/2 to fs/2
    plot(fnew, fftshift(abs(Sx)));          % plots the magnitude spectrum of the message signal
    xlabel('frequency (Hz)');
    ylabel('power spectrum');
    grid on;
end