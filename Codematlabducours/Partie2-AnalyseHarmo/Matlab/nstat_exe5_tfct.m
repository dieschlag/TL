% S. Rossignol -- 29/07/13

clear;
close all;

%pkg load signal

fe=16000;

signal=[zeros(1,500) ones(1,1000) zeros(1,500)];
ttt=[0:length(signal)-1]/fe;
mttt=max(ttt);

figure(1);
clf;
plot(ttt,signal,'Linewidth',3)
title('signal','Fontsize',25);
ylim([-1.1 1.1]);
xlim([min(ttt) max(ttt)]);
%grid on;
xlabel('temps (s)','Fontsize',25);
hold off;
print -depsc2 nstat5_exe1_tfct.eps



tfft=8192*1; %%% si "tfft" est trop grand, "imagesc" plante
fmax=400;

%%% 25 ms
tfen=0.025;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.002));

figure(2);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('25 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',25);
ylabel('frequence (Hz)','Fontsize',25);
hold off;
print -depsc2 nstat5_exe2_tfct.eps


%%% 35 ms
tfen=0.035;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.002));

figure(3);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('35 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',25);
ylabel('frequence (Hz)','Fontsize',25);
hold off;
print -depsc2 nstat5_exe3_tfct.eps


%%% 45 ms
tfen=0.045;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.002));

figure(4);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('45 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',25);
ylabel('frequence (Hz)','Fontsize',25);
hold off;
print -depsc2 nstat5_exe4_tfct.eps


%%% 55 ms
tfen=0.055;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.002));

figure(5);
clf;
imagesc(tspec+tfen/2, fspec, (log10(abs(Sspec))));
title('55 ms ; |FFT| ; freq.>0 ; zoom au debut','Fontsize',25);
ylim([0 fmax]);
%xlim([min(tspec+tfen/2) max(tspec+tfen/2)]);
xlim([0 mttt]);
grid on;
xlabel('temps (s)','Fontsize',25);
ylabel('frequence (Hz)','Fontsize',25);
hold off;
print -depsc2 nstat5_exe5_tfct.eps

