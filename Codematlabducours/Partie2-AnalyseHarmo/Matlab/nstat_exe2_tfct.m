% S. Rossignol -- 29/07/13

clear;
close all;

fe=16000;
mttt=0.2;
tt1=[0:1/fe:mttt/2];
tt2=[length(tt1)/fe:1/fe:mttt];
ttt=[tt1 tt2];

signal = [0.5*cos(2*pi*50*tt1) 0.5*cos(2*pi*200*tt2)];

figure(1);
clf;
plot(ttt,signal,'r','Linewidth',2)
title('signal','Fontsize',20);
ylim([-1.1 1.1]);
grid on;
xlabel('temps (s)','Fontsize',20);
hold off;
print -depsc2 nstat2_exe1_tfct.eps


tfft=8192*1; %%% si "tfft" est trop grand, "imagesc" plante
fmax=400;

%%% 30 ms
tfen=0.03;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.005));

figure(2);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('30 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 nstat2_exe2_tfct.eps


%%% 60 ms
tfen=0.06;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.005));

figure(3);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('60 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 nstat2_exe3_tfct.eps


%%% 90 ms
tfen=0.09;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.005));

figure(4);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('90 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 nstat2_exe4_tfct.eps


%%% 120 ms
tfen=0.12;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.005));

figure(5);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('120 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 nstat2_exe5_tfct.eps

