%%% S. Rossignol -- 24/08/13

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


%%% Hilbert

sigcomp=hilbert(signal);

for ii=1:length(sigcomp)-1
  freqinst(ii) = (angle(sigcomp(ii+1))-angle(sigcomp(ii)))/2/pi*fe;

  %%% une légère correction pour cas aberrants
  if (freqinst(ii)<0.)
    freqinst(ii)=freqinst(ii-1);
  end;
end;

figure(2);
clf;
plot(ttt(1:end-1),freqinst,'r','Linewidth',3)
title('frequence instantanee','Fontsize',20);
ylim([40 220]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 freqinst.eps


%%% spectrogrammes

tfft=8192*16;
fmax=400;

%%% 30 ms
tfen=0.03;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.001));

for ii=1:size(Sspec,2)
  [vv,pp]=max(abs(Sspec(:,ii)));
  f0spec1(ii)=fspec(pp);
end;

figure(3);
clf;
plot(tspec+tfen/2, f0spec1,'Linewidth',3);
title('frequence spectrogramme (30ms)','Fontsize',20);
xlim([0 mttt]);
ylim([40 220]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 freqspec.eps


%%% 60 ms
tfen=0.06;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.001));

for ii=1:size(Sspec,2)
  [vv,pp]=max(abs(Sspec(:,ii)));
  f0spec2(ii)=fspec(pp);
end;

figure(4);
clf;
plot(tspec+tfen/2, f0spec2,'Linewidth',3);
title('frequence spectrogramme (60ms)','Fontsize',20);
xlim([0 mttt]);
ylim([40 220]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;


%%% 90 ms
tfen=0.09;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.001));

for ii=1:size(Sspec,2)
  [vv,pp]=max(abs(Sspec(:,ii)));
  f0spec3(ii)=fspec(pp);
end;

figure(5);
clf;
plot(tspec+tfen/2, f0spec3,'Linewidth',3);
title('frequence spectrogramme (90ms)','Fontsize',20);
xlim([0 mttt]);
ylim([40 220]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;


%%% 120 ms
tfen=0.12;lfen=round(fe*tfen);
[Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.001));

for ii=1:size(Sspec,2)
  [vv,pp]=max(abs(Sspec(:,ii)));
  f0spec4(ii)=fspec(pp);
end;

figure(6);
clf;
plot(tspec+tfen/2, f0spec4,'Linewidth',3);
title('frequence spectrogramme (120ms)','Fontsize',20);
xlim([0 mttt]);
ylim([40 220]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;

