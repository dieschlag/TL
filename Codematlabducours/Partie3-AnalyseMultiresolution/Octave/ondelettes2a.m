%%% Quelques exemples de transformation en ondelettes (ondelettes continues)
%%%
%%% S. Rossignol -- 04/08/2013, puis 2021

clear all;
close all;


feo=1000; %%% fréquence d'échantillonnage de l'ondelette


%%% 1 -- génératrice ondelette

sigi=[exp(-[-49.5:1/feo:49.5].^2/400)];


%%% 2 -- ondelette 1 : dérivée première de la génératrice

sigid = ([sigi 0] - [0 sigi])*feo; %%% dérivée numérique
sigid = sigid(2:length(sigi));
sigid=sigid/max(sigid); %%% normalisation
timid=([0:length(sigid)-1]-length(sigid)/2+0.5)/feo; %%% instants d'échantillonnage de l'ondelette


%%% 3 -- ondelette 2 : dérivée seconde de la génératrice

sigidd = ([sigid 0] - [0 sigid])*feo; %%% dérivée numérique de la dérivée
sigidd = sigidd(2:length(sigid));
sigidd=-sigidd/max(sigidd); %%% normalisation
timidd=([0:length(sigidd)-1]-length(sigidd)/2+0.5)/feo; %%% instants d'échantillonnage de l'ondelette


%%% on choisit l'une des deux ondelettes (une seule ligne ci-dessous doit être décommentée)
wave=sigid;timi=timid;name1='o1'; %%% 'name1' sert pour la formation du nom du fichier où est sauvée la figure ci-dessous
%wave=sigidd;timi=timidd;name1='o2';


%%%
figure(2);
clf;
plot(timi,wave,'Linewidth',3)
title('ondelette','Fontsize',25);
ylim([min(wave) max(wave)]);
hold on;
grid on;
xlabel('temps','Fontsize',25);
hold off;


%%% signal

fe=16000; %%% fréquence d'échantillonnage du signal traité

dosig=5;   %%% 5 signaux différents sont analysables : choisir lequel ici (le 4 est enlevé)
if (dosig==1)
  %%% voir nstat_exe1_tfct.m pour l'équivalent tfct (Transformée de Fourier à Court Terme)
  mttt=0.2;          %%% longueur du signal en seconde
  ttt=[0:1/fe:mttt]; %%% instants d'échantillonnage
  signal = 0.5*cos(2*pi*50*ttt) + 0.5*cos(2*pi*200*ttt); %%% signal 1. : somme de deux cosinus
  name2='s1'; %%% 'name2' sert pour la formation du nom du fichier où est sauvée la figure ci-dessous
elseif (dosig==2)
  %%% nstat_exe2_tfct.m
  mttt=0.2;
  tt1=[0:1/fe:mttt/2];
  tt2=[length(tt1)/fe:1/fe:mttt];
  ttt=[tt1 tt2];
  signal = [0.5*cos(2*pi*50*tt1) 0.5*cos(2*pi*200*tt2)]; %%% signal 2. : succession de deux sinus
  name2='s2';
elseif (dosig==3)
  %%% nstat_exe3_tfct.m
  mttt=0.2;
  ttt=[0:1/fe:mttt];
  f0=50+(200-50)/0.2/2.*ttt; %%% pour le /2., voir l'exercice à la fin du cours
  signal = 0.5*cos(2*pi*f0.*ttt);                        %%% signal 3. : chirp de fréquence
  name2='s3';
elseif (dosig==4)
  %%% nstat_exe4_tfct.m => pas très intéressant ici

  name2='s4';
elseif (dosig==5)
  %%% nstat_exe5_tfct.m
  signal=[zeros(1,500) ones(1,1000) zeros(1,500)];      %%% signal 5. : une porte
  mttt=length(signal)/fe;
  name2='s5';
end;
lsig=length(signal); %%% longueur du signal en nombre d'échantillons


%%%

speci=[]; %%% initialisation de la matrice qui comprendra les spectres

nanal=400;          %%% nombre de niveau d'analyse
sdeb=30;            %%% niveau minimal (la correspondance en fréquence n'est pas faite)
sfin=40000;         %%% niveau maximal (la correspondance en fréquence n'est pas faite)
ratio=(sfin/sdeb)^(1/nanal);
saut=sdeb/ratio;
for ii=1:nanal
  saut=round(saut*ratio);

  ondelette=wave(1:saut:end);
  ondelette=ondelette*(ii^2)/(nanal^2); %%% ondelette dilatée

  %figure(2);
  %clf;
  %plot(ondelette);
  %hold on;
  %plot(signal(1:length(ondelette)));
  %hold off;
  %drawnow;
  %pause;

  spec=conv(signal,ondelette); %%% convolution entre le signal et l'ondelette (souvenez-vous des SLI)
  lspec=length(spec);          %%% longueur du spectre obtenu
  diff=round((lspec-lsig)/2);  %%% demi-différence de longueur entre le signal et le spectre
                               %%% note : la différence vient de la convolution

  spec=spec(diff:diff+lsig-1); %%% il faut bricoler pour que tous les spectres aient la même longueur,
                               %%% afin de pouvoir les concaténer (ci-dessous)
  speci=[speci spec'];
end;

figure(3);
clf;
title('module analyse en ondelettes','Fontsize',25);
hold on;
imagesc([0 mttt],[0 nanal-1],abs(speci)');
xlim([0 mttt]);
ylim([0 nanal-1]);
xlabel('temps (s)','Fontsize',25);
ylabel('frequence','Fontsize',25);
hold off;
print(['ondelettes2_1' name1 name2 '.eps'], '-depsc2');  %%% pour sauver automatiquement les figures
                                                         %%% dans un fichier qu'on peut inclure après
                                                         %%% dans un rapport ou une présentation, il faut
                                                         %%% utiliser la commande "print" : voir l'aide

