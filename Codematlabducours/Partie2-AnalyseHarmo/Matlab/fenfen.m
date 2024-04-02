%%% S. Rossignol -- 31/07/13

clear;
close all;

lll=8192;
ll2=lll/2;
fe=8000;

timi=[-150:150]/fe;


% hanning

sigi1=hanning(length(timi))';sigi1=sigi1/sum(sigi1);
ggg1=fft(sigi1,lll);
hhh1=fftshift(ggg1);

ttt=[0:lll-1]/lll*(length(sigi1)/2-0.5);
decalage=ttt*pi*2;

fff1=ggg1.*exp(j*decalage);
fff1=fftshift(fff1);

% hamming

sigi2=hamming(length(timi))';sigi2=sigi2/sum(sigi2);
ggg2=fft(sigi2,lll);
hhh2=fftshift(ggg2);

fff2=ggg2.*exp(j*decalage);
fff2=fftshift(fff2);

% blackman

sigi3=blackman(length(timi))';sigi3=sigi3/sum(sigi3);
ggg3=fft(sigi3,lll);
hhh3=fftshift(ggg3);

fff3=ggg3.*exp(j*decalage);
fff3=fftshift(fff3);

% rectangulaire

sigi4=ones(1,length(timi));sigi4=sigi4/sum(sigi4);
ggg4=fft(sigi4,lll);
hhh4=fftshift(ggg4);

fff4=ggg4.*exp(j*decalage);
fff4=fftshift(fff4);


freq=([0:lll-1]-lll/2)*fe/lll;

figure(1);
clf;
plot(timi,sigi1,'b','Linewidth',4);
hold on;
plot(timi,sigi2,'r','Linewidth',4);
plot(timi,sigi3,'g','Linewidth',4);
plot(timi,sigi4,'k','Linewidth',4);
title('signal','Fontsize',20);
ylim([-0.0003 0.0079365]);
grid on;
xlabel('temps','Fontsize',20);
legend('hanning','hamming','blackman','rectangulaire');
hold off;
print -depsc2 fenfen_1.eps

zzz=200;

figure(2);
clf;
plot(freq,20*log10(abs(fff1)),'b','Linewidth',4);
hold on;
plot(freq,20*log10(abs(fff2)),'r','Linewidth',4);
plot(freq,20*log10(abs(fff3)),'g','Linewidth',4);
plot(freq,20*log10(abs(fff4)),'k','Linewidth',4);
title('amplitude en dB (zoom au centre)','Fontsize',20);
grid on;
xlim([-zzz zzz]);
ylim([-100 1]);
xlabel('frequence (Hz)','Fontsize',20);
legend('hanning','hamming','blackman','rectangulaire');
hold off;
print -depsc2 fenfen_2.eps

