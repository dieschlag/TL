%%% S. Rossignol -- 24/08/13

clear;
close all;

fe=16000;
mttt=0.2;
tt1=[0:1/fe:mttt/2];
tt2=[length(tt1)/fe:1/fe:mttt];
ttt=[tt1 tt2];

signal = [0.5*cos(2*pi*50*tt1) 0.5*cos(2*pi*200*tt2)];
signal = signal + 0.004*randn(1,length(ttt));

figure(1);
clf;
plot(ttt,signal,'Linewidth',2,'r')
title('signal','Fontsize',20);
ylim([-1.1 1.1]);
grid on;
xlabel('temps (s)','Fontsize',20);
hold off;
print -depsc2 sigbruit.eps


%%% Hilbert

sigcomp=hilbert(signal);

figure(2);
clf;
plot(ttt,signal,'Linewidth',2,'r');
hold on;
plot(ttt,imag(sigcomp),'Linewidth',2,'r');
title('signal','Fontsize',20);
ylim([-1.1 1.1]);
grid on;
xlabel('temps (s)','Fontsize',20);
hold off;

for ii=1:length(sigcomp)-1
  freqinst(ii) = (arg(sigcomp(ii+1))-arg(sigcomp(ii)))/2/pi*fe;

  %%% une légère correction pour cas aberrants
  if (freqinst(ii)<0.)
    freqinst(ii)=freqinst(ii-1);
  end;
end;

figure(3);
clf;
plot(ttt(1:end-1),freqinst,'Linewidth',3,'r')
title('frequence instantanee','Fontsize',20);
ylim([0 260]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 freqinstbruit.eps


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

figure(4);
clf;
plot(tspec+tfen/2, f0spec1,'Linewidth',3);
title('frequence spectrogramme (30ms)','Fontsize',20);
xlim([0 mttt]);
ylim([0 260]);
grid on;
xlabel('temps (s)','Fontsize',20);
ylabel('frequence (Hz)','Fontsize',20);
hold off;
print -depsc2 freqspecbruit.eps
