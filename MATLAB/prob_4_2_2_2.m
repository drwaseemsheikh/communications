% Problem 4.2-2(ii) plots
close all;
clear;

fs = 4000;
Ts = 1/fs;
T = 1;
t = 0:Ts:T;
df = 1/100;
f1 = -10:df:10;
f2 = -1020:df:1020;

% message signal
figure;
m = 2*exp(-2*t);
subplot(3,1,1);
plot(t,m);
xlabel('time (s)');
ylabel('message signal');
grid on;

% message signal spectrum
M = 1./(1+j*pi*f1);
subplot(3,1,2);
plot(f1,abs(M));
xlabel('frequency (Hz)');
ylabel('magnitude spectrum (message)');
grid on;

subplot(3,1,3);
plot(f1,angle(M));
xlabel('frequency (Hz)');
ylabel('phase spectrum (message)');
grid on;

%modulated signal
figure;
s = 2*m.*cos(2000*pi*t);
subplot(3,1,1);
plot(t,s);
grid on;
axis([0 0.1 -5 5]);
xlabel('time (s)');
ylabel('modulated signal');
grid on;

%modulated signal spectrum
S = 1./(1+j*pi*(f2+1000)) + 1./(1+j*pi*(f2-1000));
subplot(3,1,2);
plot(f2,abs(S));
xlabel('frequency (Hz)');
ylabel('magnitude spectrum (modulated)');
grid on;

subplot(3,1,3);
plot(f2,angle(S));
xlabel('frequency (Hz)');
ylabel('phase spectrum (message)');
grid on;
