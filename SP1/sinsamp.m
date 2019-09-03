%
% sinsamp.m
%
% This program will sample a sinusoidal signal
% and plot the result.
%
clear all;
samptime=1/1000;
freq=10;
for ii=1:1000
    time(ii)=(ii-1)*samptime;
    x(ii)=cos(2*pi*freq*time(ii));
end
subplot(2,1,1);
plot(time,x);
title('Sampled Sinusoid with loop');
xlabel('Time (sec)');
ylabel('x(t)');
grid on;
%
% accomplished above but w/o loop
%
clear all;
samptime=1/1000;
freq=10;
time=[0:999]*samptime;
x=cos(2*pi*freq*time);
subplot(2,1,2);
plot(time,x);
title('Sampled Sinusoid w/o loop');
xlabel('Time (sec)');
ylabel('x(t)');
grid on;