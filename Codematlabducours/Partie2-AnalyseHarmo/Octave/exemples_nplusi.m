%%% S. Rossignol -- 11/2008


clear all;
close all;


%%% influence :
%%% 1- de la position fr\'equentielle
%%% 2- du 0-padding
%%% 3- du fen\^etrage
%%% 4- du fen\^etrage avec plusieurs sinuso\"{\i}des
do1 = 1;
do2 = 1; %%% needs 1
do3 = 1;
do4 = 1;


%%%%%% pour 1, 2 et 3
fe = 512;
nsig=1024;
ttt = [1:nsig]/fe;
nfi=10000;
ffi = [0:nfi-1]/nfi*fe;


%%%%%%
f0a = 50.00;
f0b = 50.25;
f0c = 51.75;
f0d = 53.25;
signala = cos(2*pi*f0a*ttt);
signalb = cos(2*pi*f0b*ttt);
signalc = cos(2*pi*f0c*ttt);
signald = cos(2*pi*f0d*ttt);



%%%
if (do1==1)


%%% 1 %%%


spectre = abs(fft(signala))*2/length(signala);
fff = [0:length(spectre)-1]/length(spectre)*fe;
rspectre = spectre;
rfff = fff;

speth = abs(sinc((nsig/fe)*(ffi-f0a)));

figure(1);
clf;
plot(ttt(1:100),signala(1:100), 'LineWidth', 3);
xlabel('temps en secondes','FontSize',20);
ylabel('amplitude','FontSize',20);
title('portion du signal : sinus de frequence 50 Hz','FontSize',20);
hold off;
print -depsc2 sigech1.eps

figure(2);
clf;
plot(ffi,speth,'r-.', 'LineWidth', 3);
hold on;
plot(fff,spectre,'-o', 'markersize', 5, 'LineWidth', 3);
%xlim([0 fe/2]);
xlim([45 55]);
xlabel('frequence (Hz)','FontSize',20);
ylabel('amplitude','FontSize',20);
title('fe=512 Hz, f0=50 Hz, Nsig=1024, pas de 0-padding','FontSize',20);
legend('TF','TFD');
grid on;
hold off;
print -depsc2 spec1.eps


%%%%%%
spectre = abs(fft(signalb))*2/length(signalb);
fff = [0:length(spectre)-1]/length(spectre)*fe;

speth = abs(sinc((nsig/fe)*(ffi-f0b)));

figure(3);
clf;
plot(ttt(1:100),signalb(1:100), 'LineWidth', 3);
xlabel('temps en secondes','FontSize',20);
ylabel('amplitude','FontSize',20);
title('portion du signal : sinus de frequence 50.25 Hz','FontSize',20);
hold off;
print -depsc2 sigech2.eps

figure(4);
clf;
plot(ffi,speth,'r-.', 'LineWidth', 3);
hold on;
plot(fff,spectre,'-o','markersize', 5, 'LineWidth', 3);
%xlim([0 fe/2]);
xlim([45 55]);
xlabel('frequence (Hz)','FontSize',20);
ylabel('amplitude','FontSize',20);
title('fe=512 Hz, f0=50.25 Hz, Nsig=1024, pas de 0-padding','FontSize',20);
legend('TF','TFD');
grid on;
hold off;
print -depsc2 spec2.eps

spectresans0 = spectre;
fffsans0     = fff;

figure(5);
clf;
%plot(ffi,speth,'r-.', 'LineWidth', 3);
plot(fff,spectre,'-or','markersize', 5, 'LineWidth', 3);
hold on;
plot(rfff,rspectre,'-og','markersize', 5, 'LineWidth', 3);
xlim([45 55]);
xlabel('frequence (Hz)','FontSize',20);
ylabel('amplitude','FontSize',20);
title('fe=512 Hz, Nsig=1024, pas de 0-padding','FontSize',20);
legend('f=50.25 Hz','f=50 Hz');
grid on;
hold off;
print -depsc2 spec3.eps

end; %%% of 1



%%%
if do2==1


%%% 2 %%%


nfft=4096;

%%%%%%
spectre = abs(fft(signala,nfft))*2/length(signala);
fff = [0:length(spectre)-1]/length(spectre)*fe;

speth = abs(sinc((nsig/fe)*(ffi-f0a)));

figure(6);
clf;
plot(ffi,speth,'r-.', 'LineWidth', 3);
hold on;
plot(fff,spectre,'-o', 'LineWidth', 3);
%xlim([0 fe/2]);
xlim([45 55]);
xlabel('frequence (Hz)','FontSize',20);
ylabel('amplitude','FontSize',20);
title('fe=512 Hz, f0=50 Hz, Nsig=1024, Nfft=4096','FontSize',20);
legend('TF','TFD');
grid on;
hold off;


%%%%%%
spectre = abs(fft(signalb,nfft))*2/length(signalb);
fff = [0:length(spectre)-1]/length(spectre)*fe;

speth = abs(sinc((nsig/fe)*(ffi-f0b)));

figure(7);
clf;
plot(ffi,speth,'r-.', 'LineWidth', 3);
hold on;
plot(fff,spectre,'-o','markersize', 5, 'LineWidth', 3);
xlim([45 55]);
xlabel('frequence (Hz)','FontSize',20);
ylabel('amplitude','FontSize',20);
title('fe=512 Hz, f0=50.25 Hz, Nsig=1024, Nfft=4096','FontSize',20);
legend('TF','TFD');
grid on;
hold off;
print -depsc2 spec4.eps

figure(14);
clf;
plot(fffsans0,spectresans0,'-or','markersize', 5, 'LineWidth', 3);
hold on;
plot(fff,spectre,'-og','markersize', 5, 'LineWidth', 3);
xlim([45 55]);
xlabel('frequence (Hz)','FontSize',20);
ylabel('amplitude','FontSize',20);
title('fe=512 Hz, Nsig=1024, f0=50.25 Hz','FontSize',20);
legend('pas de 0-padding','avec 0 padding');
grid on;
hold off;
print -depsc2 spec5.eps

figure(15);
clf;
plot(ttt,signalb,'-or','markersize', 3, 'LineWidth', 3);
hold on;
xlabel('temps (s)','FontSize',20);
title('fe=512 Hz, Nsig=1024, f0=50.25 Hz','FontSize',20);
legend('pas de 0-padding');
grid on;
hold off;
print -depsc2 signal_nplusi1.eps

figure(16);
clf;
plot([1:4096]/fe,[signalb zeros(1,3072)],'-or','markersize', 3, 'LineWidth', 3);
hold on;
xlabel('temps (s)','FontSize',20);
title('fe=512 Hz, Nsig=1024, Nfft=4096, f0=50.25 Hz','FontSize',20);
legend('0-padding');
grid on;
hold off;
print -depsc2 signal_nplusi2.eps

end; %%% of 2



%%%
if do3==1


%%% 3 %%%


  %%%%%%
  tfft=length(signala);

  spectre = abs(fft(signala.*hanning(length(signala))',tfft))/length(signala);
  fff = [0:length(spectre)-1]/length(spectre)*fe;

  sperec = abs(fft(signala,tfft))/length(signala);

  figure(8);
  clf;
  plot(fff,sperec,'r-.*', 'markersize', 0, 'LineWidth', 3);
  hold on;
  plot(fff,spectre,'-o', 'markersize', 0, 'LineWidth', 3);
  xlim([45 55]);
  xlabel('frequence (Hz)','FontSize',20);
  ylabel('amplitude','FontSize',20);
  title('fe=512 Hz, f0=50 Hz, Nsig=1024, pas de 0-padding','FontSize',20);
  legend('rectangulaire','Hanning');
  grid on;
  hold off;
  print -depsc2 ex3_spec1.eps

  figure(9);
  clf;
  plot(fff,20*log10(sperec),'r-.*', 'markersize', 0, 'LineWidth', 3);
  hold on;
  plot(fff,20*log10(spectre),'-o', 'markersize', 0, 'LineWidth', 3);
  xlim([45 55]);
  ylim([-110 0]);
  xlabel('frequence (Hz)','FontSize',20);
  ylabel('amplitude (dB)','FontSize',20);
  title('fe=512 Hz, f0=50 Hz, Nsig=1024, pas de 0-padding','FontSize',20);
  legend('rectangulaire','Hanning');
  grid on;
  hold off;
  print -depsc2 ex3_spec2.eps

  %%%%%%
  tfft=length(signala)*4;

  spectre2 = abs(fft(signala.*hanning(length(signala))',tfft))/length(signala);
  fff2 = [0:length(spectre2)-1]/length(spectre2)*fe;

  sperec2 = abs(fft(signala,tfft))/length(signala);

  figure(10);
  clf;
  plot(fff2,sperec2,'r-.*', 'markersize', 0, 'LineWidth', 3);
  hold on;
  plot(fff2,spectre2,'-o', 'markersize', 0, 'LineWidth', 3);
  xlim([45 55]);
  xlabel('frequence (Hz)','FontSize',20);
  ylabel('amplitude','FontSize',20);
  title('fe=512 Hz, f0=50 Hz, Nsig=1024, NFFT=4096','FontSize',20);
  legend('rectangulaire','Hanning');
  grid on;
  hold off;
  print -depsc2 ex3_spec3.eps

  figure(11);
  clf;
  plot(fff2,20*log10(sperec2),'r-.*', 'markersize', 0, 'LineWidth', 3);
  hold on;
  plot(fff2,20*log10(spectre2),'-o', 'markersize', 0, 'LineWidth', 3);
  xlim([45 55]);
  ylim([-110 0]);
  xlabel('frequence (Hz)','FontSize',20);
  ylabel('amplitude (dB)','FontSize',20);
  title('fe=512 Hz, f0=50 Hz, Nsig=1024, NFFT=4096','FontSize',20);
  legend('rectangulaire','Hanning');
  grid on;
  hold off;
  print -depsc2 ex3_spec4.eps

end; %%% of 3



if do4==1


  %%% 4 %%%
  ttt = ([1:nsig]-(1+nsig)/2)/fe;
  signala = cos(2*pi*f0a*ttt);
  signalc = cos(2*pi*f0c*ttt);
  signald = cos(2*pi*f0d*ttt);

  nfft=4096*2;


  %%%%%%
  aa = 1.00;
  ac = 0.06;
  ad = 0.04;
  signals = aa*signala+ac*signalc+ad*signald;


  %%%%%%
  win = ones(1,length(signals));


  %%%%%%
  spectre = abs(fft(signals.*win,nfft))/length(signals);
  fff = [0:length(spectre)-1]/length(spectre)*fe;

  figure(12);
  clf;grid on;hold on;
  xlim([45 59]);
  ylim([-80 5]);
%  xlim([48 52]);   %%% pour r\'egler
%  ylim([-27 -18]); %%% pour r\'egler
  xlabel('frequence (Hz)','FontSize',20);
  ylabel('amplitude (dB)','FontSize',20);
  title('fe=512 Hz, 3 sinus, Nsig=1024, NFFT=8192, Rectangulaire','FontSize',20);
  nfig=6;
  mksize= 25;
  for ii=1:nfig
    if ii==1
      plot(fff,20*log10(spectre),'-', 'LineWidth', 3);
      legend('Spectre (rectangulaire)','location','northeast');
    elseif ii==2
      plot(f0a, 20*log10(aa/2), '-.or', 'markersize', mksize, 'LineWidth', 3);
      legend('Spectre (rectangulaire)','position composante 1','location','northeast');
    elseif ii==3
      plot(f0d, 20*log10(ad/2),   'or', 'markersize', mksize, 'LineWidth', 3);
      legend('Spectre (rectangulaire)','position composante 1','position composante 2','location','northeast');
    elseif ii==4
      plot([fff(1) 50 fff(length(fff))],[-19 -19 -19], '-k', 'LineWidth', 4);   %%% seuil
      legend('Spectre (rectangulaire)','position composante 1','position composante 2','Seuil 1 - trop haut : la composante 2 manque','location','northeast');
    elseif ii==5
      plot([fff(1) 50 fff(length(fff))],[-26.5 -26.5 -26.5], '-k', 'LineWidth', 4); %%% seuil
      plot(51.25, -23.9, '*g', 'markersize', mksize, 'LineWidth', 3);
      plot(48.75, -24.0, '*g', 'markersize', mksize, 'LineWidth', 3);
      plot(50.72, -19.5, '*g', 'markersize', mksize, 'LineWidth', 3);
      plot(49.28, -19.5, '*g', 'markersize', mksize, 'LineWidth', 3);
      legend('Spectre (rectangulaire)','position composante 1','position composante 2','Seuil 1 - trop haut : la composante 2 manque','Seuil 2 - trop bas : il y a 4 fausses alarmes','fausse alarme 1','fausse alarme 2','fausse alarme 3','fausse alarme 4','location','northeast');
    elseif ii==6
      plot(f0c, 20*log10(ac/2), '--or', 'markersize', mksize, 'LineWidth', 3);
      legend('Spectre (rectangulaire)','position composante 1','position composante 2','Seuil 1 - trop haut : la composante 2 manque','Seuil 2 - trop bas : il y a 4 fausses alarmes','fausse alarme 1','fausse alarme 2','fausse alarme 3','fausse alarme 4','position composante 3','location','northeast');
    end;
    %plot(fff,20*log10(spectre),'-', 'LineWidth', 3); %%% redessin

    command = ['print ' 'ex4_spec1-' num2str(ii) '.eps -depsc2'];
    eval(command);
    %pause;
  end;
  hold off;
  print -depsc2 ex4_spec1.eps


  %%%%%%
  win = hanning(length(signals))';


  %%%%%%
  spectre = abs(fft(signals.*win,nfft))/length(signals);
  fff = [0:length(spectre)-1]/length(spectre)*fe;

  figure(13);
  clf;grid on;hold on;
  plot(fff,20*log10(spectre),'-', 'LineWidth', 3);
  plot([fff(1) 50 fff(length(fff))],[-42 -42 -42], '-k', 'LineWidth', 4); %%% seuil
  plot(f0a, 20*log10(aa/4), '-.or', 'markersize', mksize, 'LineWidth', 3);
  plot(f0c, 20*log10(ac/4),   'or', 'markersize', mksize, 'LineWidth', 3);
  plot(f0d, 20*log10(ad/4), '--or', 'markersize', mksize, 'LineWidth', 3);
  plot(fff,20*log10(spectre),'-', 'LineWidth', 3); %%% redessin
  xlim([45 59]);
  ylim([-80 5]);
  xlabel('frequence (Hz)','FontSize',20);
  ylabel('amplitude (dB)','FontSize',20);
  title('fe=512 Hz, 3 sinus, Nsig=1024, NFFT=8192, Hanning','FontSize',20);
  legend('Spectre (Hanning)','Seuil - correct','position composante 1','position composante 2','position composante 3');
  hold off;
  print -depsc2 ex4_spec2.eps

end;

