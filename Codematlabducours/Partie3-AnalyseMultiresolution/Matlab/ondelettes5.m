%%% Ondelettes continues
%%%
%%% S. Rossignol -- 2020/2021

clear all;
close all;


%%% signal
%pkg load signal

fe=16000;   %%% fréquence d'écahntillonnage des sons
flim=400;   %%% fréquence max limite observée 
            %%%   => modifiée parfois plus loin
taille=100; %%% 'taille' de l'ondelette
            %%%   => modifiée parfois plus loin
name1='o1'; %%% ondelette de Morlet

dosig=7;
if (dosig==0)
  %%% un sinus tout seul
  mttt=0.1;
  ttt=[0:1/fe:mttt];
  signal = 0.5*cos(2*pi*100*ttt);
  name2='s0';
  sound(signal,fe);
elseif (dosig==1)
  %%% nstat_exe1_tfct.m
  mttt=0.2;
  ttt=[0:1/fe:mttt];
  signal = 0.5*cos(2*pi*50*ttt) + 0.5*cos(2*pi*200*ttt); %%% 2 sinus
  name2='s1';
  flim=1000;

  %%% pour avoir un son suffisamment long pour être audible
  mt2=10*mttt;
  tt2=[0:1/fe:mt2];
  sig2 = 0.5*cos(2*pi*50*tt2) + 0.5*cos(2*pi*200*tt2);
  sound(sig2*0.5,fe);
elseif (dosig==2)
  %%% nstat_exe2_tfct.m
  mttt=0.2;
  tt1=[0:1/fe:mttt/2];
  tt2=[length(tt1)/fe:1/fe:mttt];
  ttt=[tt1 tt2];
  signal = [0.5*cos(2*pi*50*tt1) 0.5*cos(2*pi*200*tt2)];
  name2='s2';
  taille=12;

  %%% pour avoir un son suffisamment long pour être audible
  mt2=10*mttt;
  ttt1=[0:1/fe:mt2/2];
  ttt2=[length(ttt1)/fe:1/fe:mt2];
  sig2 = [0.5*cos(2*pi*50*ttt1) 0.5*cos(2*pi*200*ttt2)];
  sound(sig2*0.5,fe);
elseif (dosig==3)
  %%% nstat_exe3_tfct.m
  mttt=0.2;
  ttt=[0:1/fe:mttt];
  f0=50+(200-50)/0.2/2.*ttt; %%% pour le /2., voir l'exercice à la fin du cours de Sig1
  signal = 0.5*cos(2*pi*f0.*ttt);
  name2='s3';
  taille=5;

  %%% pour avoir un son suffisamment long pour être audible
  mt2=10*mttt;
  tt2=[0:1/fe:mt2];
  f0=50+(200-50)/0.2/2.*tt2;
  sig2= 0.5*cos(2*pi*f0.*tt2);
  sound(sig2*0.5,fe);
elseif (dosig==4)
  %%% nstat_exe4_tfct.m => pas très intéressant ici

  name2='s4';
elseif (dosig==5)
  %%% nstat_exe5_tfct.m
  signal=[zeros(1,500) ones(1,1000) zeros(1,500)];
  mttt=length(signal)/fe;
  name2='s5';

  %%% son pas très audible
elseif (dosig==6)
  %%% signal sonore~: un sinus avec un vibrato
  f0c=440;  %%% fréquence centrale du signal (en Hz)
  Avib=20;  %%% amplitude de la modulation (en Hz)
  fvib=5;   %%% fréquence de la modulation (en Hz)
  phivib=0; %%% phase initiale de la modulation

  mttt=0.4;
  ttt=[0:1/fe:mttt]'; %%% on va faire un signal long de 0.4 s

  phi = 2*pi*(f0c*ttt + Avib/2/pi/fvib*sin(2*pi*fvib*ttt+phivib));
  signal=0.5*cos(phi);

  flim=800;

  name2='s6';
  sound(signal,fe);
elseif (dosig==7)
  %%% signal sonore~: quatre sinus avec un vibrato
  f0c=440;  %%% fréquence centrale du signal (en Hz)
  Avib=20;  %%% amplitude de la modulation (en Hz)
  fvib=5;   %%% fréquence de la modulation (en Hz)
  phivib=0; %%% phase initiale de la modulation

  mttt=2.2; %%% 0.2 : trop court pour être bien perçu
  ttt=[0:1/fe:mttt]'; %%% on va faire un signal long de 0.4 s

  phi = 2*pi*(f0c*ttt + Avib/2/pi/fvib*sin(2*pi*fvib*ttt+phivib));
  signal=0.5*cos(phi);
  nii=4;
  for ii=2:nii
    phi = 2*pi*(ii*f0c*ttt + ii*Avib/2/pi/fvib*sin(2*pi*fvib*ttt+phivib));
    signal = signal + 0.5*cos(phi);
  end;
  signal = signal/max(abs(signal));

  flim=2000;

  name2='s7';
  sound(signal*0.5,fe);
end;


y=signal;
nn=length(signal);

fmin=0.001;
fmax=0.499;

tfft1=2^17;
ffs=[0:tfft1-1]/tfft1*fe;

tscalo=8192; %%% pour se mettre à égalité avec le spectrogramme (ci-dessous)
fff = linspace(fmin,fmax,tscalo);
aaa = max(fff)./fff;

W=zeros(tscalo,nn+1);

for ii=tscalo:-1:1
  aa=aaa(ii);

  at = taille*aa;
  xx = -round(at):round(at);
  g  = exp(-(2*log(10)/at^2)*xx.^2).*exp(i*2*pi*fff(ii)*xx);

  % calcul des coefficients d'ondelettes à l’échelle a
  wa=conv(y,g);
  lwa=length(wa);
  W(ii,1:nn+1)=abs(wa(lwa/2-nn/2:lwa/2+nn/2))/sqrt(aa);

  % essai d'équivalence canal -> fréquence
  ssp=fft(g,tfft1);
  ssp=ssp(1:tfft1/2);
  [mm,posmm]=max(ssp);
  freq(ii)=ffs(posmm);
end;

%%%figures
colormap jet;

figure(1);
clf;
imagesc([0 nn]/fe,fff*fe,W(1:tscalo,1:nn));
ylim([fmin*fe flim]);
title('Ondelette de Morlet','Fontsize',15);
xlabel('temps (s)','Fontsize',15);
ylabel('frequence (Hz)','Fontsize',15);
hold off;

figure(2);
clf;
imagesc([0 nn]/fe,[1:tscalo],W(1:tscalo,1:nn));
title('coefficients ondelettes');
xlabel('temps (en s)');
ylabel('numero ordre ondelette');
hold off;


%%% spectrogrammes supplémentaires
if dosig>=6
  tfft=8192; %%% si "tfft" est trop grand, "imagesc" plante
  tfen=0.04;lfen=round(fe*tfen);
  [Sspec, fspec, tspec] = specgram(signal, tfft, fe, hanning(lfen), lfen-round(fe*0.002));
  figure(3);
  clf;
  imagesc(tspec+tfen/2, fspec, ((abs(Sspec))));
  title('spectrogramme','Fontsize',15);
  ylim([0 flim]);
  xlim([0 mttt]);
  xlabel('temps (s)','Fontsize',15);
  ylabel('frequence (Hz)','Fontsize',15);
  hold off;
  print(['specc' name2 '.eps'],'-depsc2');
end;

doit=0; %%% ça prend du temps de sauver le 'surface' en postscript : pourquoi ?
        %%% => c'est mieux quand c'est un 'imagesc'
if doit==1
  figure(1)
  print(['ondelettenew' name1 name2 '.eps'], '-depsc2');
end;


fprintf(1,'\nErreur max en frequence : %f pour cent\n', max(abs(freq/fe-fff)./fff)*100.);


%%% on regarde un spectre et un scale
if dosig==7
  for jj=1:5:size(Sspec,2)
    spepe=abs(Sspec(:,jj));
    figure(5);
    clf;
    grid on;
    hold on;
    for ii=1:nii
      plot([ii*f0c ii*f0c],[min(spepe) max(spepe)],'r');
      plot([ii*(f0c+Avib) ii*(f0c+Avib)],[min(spepe) max(spepe)],'k');
      plot([ii*(f0c-Avib) ii*(f0c-Avib)],[min(spepe) max(spepe)],'k');
    end;
    plot(fspec,spepe);
    xlim([0 2000]);
    hold off;
    %print -depsc2 spespevib1.eps
    pause;
  end;

  for jj=1:100:size(W,2)
    scalo=W(:,jj);
    figure(4);
    clf;
    grid on;
    hold on;
    for ii=1:nii
      plot([ii*f0c ii*f0c],[min(scalo) max(scalo)],'r');
      plot([ii*(f0c+Avib) ii*(f0c+Avib)],[min(scalo) max(scalo)],'k');
      plot([ii*(f0c-Avib) ii*(f0c-Avib)],[min(scalo) max(scalo)],'k');
    end;
    plot(freq,scalo);
    xlim([0 2000]);
    ylim([0 max(scalo)]);
    hold off;
    print -depsc2 scascavib1.eps
    pause;
  end;
end;

