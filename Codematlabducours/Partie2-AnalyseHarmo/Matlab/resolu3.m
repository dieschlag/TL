%%% S. Rossignol -- 31/07/13

clear;
close all;

lll=8192;
ll2=lll/2;
fe=8000;

timi=[-150:150]/fe;
sigi=0.08*cos(2*pi*115*timi)+1.0*cos(2*pi*185*timi);
ggg=fft(sigi,lll);
hhh=fftshift(ggg);

ttt=[0:lll-1]/lll*(length(sigi)/2-0.5);
decalage=ttt*pi*2;
fff=ggg.*exp(j*decalage);

fff=fftshift(fff);

freq=([0:lll-1]-lll/2)*fe/lll;

figure(1);
clf;
stem(timi,sigi);
title('signal','Fontsize',20);
ylim([-1.1 1.1]);
grid on;
xlabel('temps','Fontsize',20);
hold off;
print -depsc2 resolu3_1.eps

zzz=300;

figure(2);
clf;
plot(freq,abs(fff),'Linewidth',2);
title('amplitude (zoom au debut)','Fontsize',20);
grid on;
xlim([0 zzz]);
xlabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 resolu3_2.eps

