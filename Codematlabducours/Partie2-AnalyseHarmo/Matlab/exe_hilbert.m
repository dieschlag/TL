%%% S. Rossignol -- 23/08/13

clear;
close all;

fe=8000; % fréquence d'échantillonnage
f0=500;  % fréquence de la cosinusoïde

tt = [0:1/fe:1];

tfft=8192*2; % taille FFT
ff=[0:tfft-1]/tfft*fe-fe/2;


%%% Signal original

xx=cos(2*pi*f0*tt);

spec=abs(fft(xx,tfft));
spec=fftshift(spec);

figure(1);
clf;
plot(ff,20*log10(spec));
xlim([-1000 1000]);
ylim([-10 75]);
xlabel('frequence en Hz','Fontsize',25);
ylabel('module du spectre en dB','Fontsize',25);
title('spectre signal original','Fontsize',25);
grid on;
hold off;
print -depsc2 hilbertorig.eps


%%% Hilbert

hh = hilbert(xx);

spec=abs(fft(hh,tfft));
spec=fftshift(spec);

figure(2);
clf;
plot(ff,20*log10(spec));
xlim([-1000 1000]);
ylim([-10 80]);
xlabel('frequence en Hz','Fontsize',25);
ylabel('module du spectre en dB','Fontsize',25);
title('spectre signal analytique','Fontsize',25);
grid on;
hold off;
print -depsc2 hilberthilb.eps



figure(3);
clf;
plot(tt,xx,'Linewidth',3);
hold on;
plot(tt,imag(hh),'r','Linewidth',3);
xlim([0.50 0.51]);
ylim([-1.1 1.1]);
xlabel('temps en s','Fontsize',25);
title('trame x(t) bleu/sombre ; sigma(t) rouge/clair','Fontsize',25);
grid on;
hold off;
print -depsc2 hilbertsig.eps

