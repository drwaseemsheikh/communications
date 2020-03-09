function [bw] = powerbandwidth(f,fs,fc,df,S,x)
%PowerBandwidth Computes the x% power bandwidth of a signal 
%   f  = frequency vector
%   fs = sampling frequency
%   fc = carrier frequency
%   df = frequency width per sample
%   S  = spectrum of the signal
%   x  = percentage of total signal power contained in the bandwidth
%   bw = x% power bandwidth of the signal with spectrum S

p = trapz(f,abs(S).^2);
fp = f(1:fs/(2*df)+1);
Sp = S(1:fs/(2*df)+1);

px = 0;
i = 0;
L = length(fp);
while (px < 0.5*x/100*p) && (i <= L)
    i = i + 1;
    Strunc = Sp(fc/df+1-i : fc/df+1+i);
    px = trapz(df, abs(Strunc).^2);
end

if i > L
    bw = 0;
end

bw = (2*i + 1)*df;
end