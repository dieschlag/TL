%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to study the behaviour of the scalogram
%
% 16/04/2007 version 1 -- 2021 version 2
% St\'ephane Rossignol
%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;


%%% initializations
%path(LOADPATH,'~/tftb-0.1/mfiles'); %%% old way
addpath('tftb-0.2/mfiles');

pkg load signal


%%% initialisations
mapmap = lagis_colormap_matlab1();


whichsig = 6;  %%% de -1 à 9
               %%% attention : il y a du mélange dans l'ordre de présentation de exemples dans les slides
taillespectre=1024;


if whichsig==-1
  %%% example coming from "tutorial.pdf", pages 52...
  tsig = 128;
  sig1 = real(anapulse(tsig)); %%% une impulsion (le signal a une partie imaginaire, qu'on vire)

  figure(1);
  clf;
  plot(sig1);
  title('signal');
  xlim([1 tsig]);
  hold off;
  print -depsc2 tftbs-1_1.ps

  fullscalogram = lagis_mytfrscalo(sig1, 1:tsig, 6, 0.05, 0.45, taillespectre, 1);
  scalogram = fullscalogram; %%% pas d'inversion
  figure(2);
  clf;
  imagesc(scalogram);
  title('impulsion - scalogram de Morlet','FontSize',20);
  colormap(mapmap);
  hold off;
  print -depsc2 tftbs-1_2.ps

  fullscalogram2 = tfrscalo(sig1, 1:tsig, 6, 0.05, 0.45, taillespectre, 1);
  scalogram2 = fullscalogram2; %%% pas d'inversion
  figure(3);
  clf;
  imagesc(scalogram2);
  title('impulsion - scalo Morlet original','FontSize',20);
  colormap(mapmap);
  hold off;
  print -depsc2 tftbs-1_2a.ps

  h1 = blackman(51);
  fullspectrogram = tfrsp(sig1, 1:tsig, 2*taillespectre, h1); %% 2* : comme ça le spectro et la scalo ont la même taille
  spectrogram = fullspectrogram(taillespectre:-1:1,:);
  figure(4);
  clf;
  imagesc([1:tsig],([taillespectre:-1:1]-1)/taillespectre/2,spectrogram);
  colormap(mapmap);
  title('impulsion - spectrogramme','FontSize',20);
  hold off;
  print -depsc2 tftbs-1_3.ps
end;


if whichsig==0
  %%% example coming from "tutorial.pdf", pages 52...
  tsig = 128;
  sig2 = real(fmconst(tsig, 0.1)) + real(fmconst(tsig, 0.2));   %%% 2 sinus

  fullscalogram = lagis_mytfrscalo(sig2, 1:tsig, 60, 0.05, 0.45, taillespectre, 1);
  scalogram = fullscalogram; %%% pas d'inversion
  figure(1);
  clf;
  imagesc(scalogram);
  title('2 sinus - scalogram de Morlet','FontSize',20);
  hold off;
  print -depsc2 tftbs0_1.ps

  fullscalogram2 = tfrscalo(sig2, 1:tsig, 60, 0.05, 0.45, taillespectre, 1);
  scalogram2 = fullscalogram2; %%% pas d'inversion
  figure(2);
  clf;
  imagesc(scalogram2);
  title('2 sinus - scalo Morlet original','FontSize',20);
  hold off;
  print -depsc2 tftbs0_1a.ps

  h1 = blackman(51);
  fullspectrogram = tfrsp(sig2, 1:tsig, 2*taillespectre, h1); %% 2* : comme ça le spectro et la scalo ont la même taille
  spectrogram = fullspectrogram(taillespectre:-1:1,:);
  figure(3);
  clf;
  imagesc([1:tsig],([taillespectre:-1:1]-1)/taillespectre/2,spectrogram);
  title('2 sinus - spectrogramme','FontSize',20);
  hold off;
  print -depsc2 tftbs0_2.ps
end;


%%% example coming from "tutorial.pdf", pages 52...
%%% and my own examples


%%% 02/05/07 : pour le "spectrogramme" seulement, afin d'am\'eliorer le contraste des images : 
%%%            "imagesc" donne des images beaucoup plus belles que celles sauvées par "saveimage"
%%% => mais pourquoi avoir besoin de faire ça ???
contrast = 1.0; % 1.0 : par d\'efaut, il n'y a pas de contraste !!!
%%% C'EST MODIFIE PLUS LOIN


sizesig    = 1280;                 % instead of 128
sizewin    =  121;                 % instead of  23
sizemorlet =   60; %round(sqrt(sizesig)); % instead of 6
if whichsig==6
  sizemorlet = 150;
end;


%%% le signal est construit
switch whichsig
  case 1
    %%% 1 sinus

    sig3 = real(fmconst(sizesig, 0.1));
    nomdusignal='un sinus';
    nomfig1='tftbs1_1.ps';
    nomfig2='tftbs1_2.ps';
    nomfig2a='tftbs1_2a.ps';

  case 2
    %%% chirp en fréquence

    sig3 = real(fmlin(sizesig, 0.05, 0.45));
    nomdusignal='un chirp';
    nomfig1='tftbs2_1.ps';
    nomfig2='tftbs2_2.ps';
    nomfig2a='tftbs2_2a.ps';

  case 3
    %%% 2 chirps

    sig3 = real(fmlin(sizesig, 0.01, 0.41) + fmlin(sizesig, 0.05, 0.45));
    nomdusignal='deux chirps';
    nomfig1='tftbs3_1.ps';
    nomfig2='tftbs3_2.ps';
    nomfig2a='tftbs3_2a.ps';

  case 4
    %%% signal avec vibrato

    midfreq = 0.0294; % fr\'equence centrale du premier harmonique
    amplmod = 0.0012; % amplitude de la modulation de fr\'equence
    nperiods = 5;     % nombre de p\'eriodes du vibrato dans la fen\^etre

    % param\`etres : taillesignal, pluspetitefreq, plusgrandefreq, p\'eriodemodulation, centremodulation, freqat0, sensdevariation
    sig3 = fmsin(sizesig, 1*(midfreq-amplmod), 1*(midfreq+amplmod), sizesig/nperiods, sizesig/2, 1*midfreq, 1);
    for ii=2:16
      sig3 = sig3 + fmsin(sizesig, ii*(midfreq-amplmod), ii*(midfreq+amplmod), sizesig/nperiods, sizesig/2, ii*midfreq, 1);
    end;
    sig3 = real(sig3);
    nomdusignal='vibrato';
    nomfig1='tftbs4_1.ps';
    nomfig2='tftbs4_2.ps';
    nomfig2a='tftbs4_2a.ps';

    contrast = 3.0;

  case 5
    %%% signal avec vibrato

    midfreq = 0.13000; % fr\'equence centrale du premier harmonique
    amplmod = 0.02;   % amplitude de la modulation de fr\'equence
    nperiods = 5;     % nombre de p\'eriodes du vibrato dans la fen\^etre

    % param\`etres : taillesignal, pluspetitefreq, plusgrandefreq, p\'eriodemodulation, centremodulation, freqat0, sensdevariation
    sig3 = fmsin(sizesig, 1*(midfreq-amplmod), 1*(midfreq+amplmod), sizesig/nperiods, sizesig/2, 1*midfreq, 1);
    for ii=2:3
      sig3 = sig3 + fmsin(sizesig, ii*(midfreq-amplmod), ii*(midfreq+amplmod), sizesig/nperiods, sizesig/2, ii*midfreq, 1);
    end;
    sig3 = real(sig3);
    nomdusignal='vibrato';
    nomfig1='tftbs5_1.ps';
    nomfig2='tftbs5_2.ps';
    nomfig2a='tftbs5_2a.ps';

    contrast = 3.0;

  case 6
    %%% ce signal est compos\'e de la somme :
    %%% - de deux sinuso\"{\i}des en basses fr\'equences et de fr\'equences tr\`es proches :
    %%%   pour les s\'eparer avec le spectrogramme, il faut prendre une GRANDE fen\^etre d'analyse
    %%% - d'une sinuso\"{\i}de en haute fr\'equence modul\'ee \`a une cadence tr\`es rapide et avec une grande amplitude :
    %%%   pour voir la modulation avec le spectrogramme, il faut prendre une PETITE fen\^etre d'analyse

    midfreq = 0.40; % fr\'equence centrale du premier harmonique
    amplmod = 0.04; % amplitude de la modulation de fr\'equence

    % param\`etres : taillesignal, pluspetitefreq, plusgrandefreq, p\'eriodemodulation, centremodulation, freqat0, sensdevariation
    sig3 = real(fmsin(sizesig, 1*(midfreq-amplmod), 1*(midfreq+amplmod), sizesig/5.0, sizesig/2, 1*midfreq, 1) + ...
           fmconst(sizesig, 0.05) + fmconst(sizesig, 0.06));
    nomdusignal='2 sin proches+1 sin module';
    nomfig1='tftbs6_1.ps';
    nomfig2='tftbs6_2.ps';
    nomfig2a='tftbs6_2a.ps';

    contrast = 1.5;

  case 7
    %%% 2 chirps

    sig3 = real(fmlin(sizesig, 0.01, 0.41) + fmlin(sizesig, 0.41, 0.01));
    nomdusignal='deux chirps';
    nomfig1='tftbs7_1.ps';
    nomfig2='tftbs7_2.ps';
    nomfig2a='tftbs7_2a.ps';

  case 8
    %%% signal avec vibrato

    midfreq = 0.08; % fr\'equence centrale du premier harmonique
    amplmod = 0.03; % amplitude de la modulation de fr\'equence

    % param\`etres : taillesignal, pluspetitefreq, plusgrandefreq, p\'eriodemodulation, centremodulation, freqat0, sensdevariation
    sigc = fmsin(sizesig, 1*(midfreq-amplmod), 1*(midfreq+amplmod), sizesig/5.0, sizesig/2, 1*midfreq, 1) + ...
           fmsin(sizesig, 2*(midfreq-amplmod), 2*(midfreq+amplmod), sizesig/5.0, sizesig/2, 2*midfreq, 1) + ...
           fmsin(sizesig, 3*(midfreq-amplmod), 3*(midfreq+amplmod), sizesig/5.0, sizesig/2, 3*midfreq, 1) + ...
           fmsin(sizesig, 4*(midfreq-amplmod), 4*(midfreq+amplmod), sizesig/5.0, sizesig/2, 4*midfreq, 1);
    sig3 = real(sigc);
    nomdusignal='vibrato';
    nomfig1='tftbs8_1.ps';
    nomfig2='tftbs8_2.ps';
    nomfig2a='tftbs8_2a.ps';

  case 9
    %%% un sinus avec une perturbation

    sig3 = 0.04*real(fmconst(sizesig, 0.1));

    sig3(302) = 0.125;
    sig3(303) = 0.25;
    sig3(304) = 0.5;
    sig3(305) = 1;
    sig3(306) = 0.5;
    sig3(307) = 0.25;
    sig3(308) = 0.125;

    nomdusignal='sinus avec perturbation';
    nomfig1='tftbs9_1.ps';
    nomfig2='tftbs9_2.ps';
    nomfig2a='tftbs9_2a.ps';

    %sig3(605) = 1;
    %sig3(905) = 1;

end;


if whichsig>0
  %%% spectrogramme
  h1 = blackman(sizewin);
  h1 = h1/sqrt(sum(h1.*h1));
  %figure(1);
  %plot(h1);
  %
  fullspectrogram = tfrsp(sig3, 1:sizesig, 2*taillespectre, h1); %% 2* : comme ça le spectro et la scalo ont la même taille
  %
  spectrogram = fullspectrogram(taillespectre:-1:1,:);
  spectrogram = spectrogram/max(max(spectrogram))*255;
  %
  spectrogramcontrasted = fullspectrogram(taillespectre/2:-1:1,:).^contrast;
  spectrogramcontrasted = spectrogramcontrasted/max(max(spectrogramcontrasted))*255;
  %
  figure(4);
  clf;
  imagesc([1:sizesig],([taillespectre:-1:1]-1)/taillespectre/2,spectrogram);
  title([nomdusignal '- spectrogramme'],'FontSize',20);
  colormap(mapmap);
  hold off;
  print(nomfig1,'-depsc2');

  if whichsig==6
    % on prend une fen\^etre deux fois plus grande pour essayer de s\'eparer les deux sinuso\"{\i}des
    % en basses fr\'equences
    h2 = blackman(3*sizewin);
    h2 = h2/sqrt(sum(h2.*h2));
    fullspectrogram = tfrsp(sig3, 1:sizesig, 2*taillespectre, h2); %% 2* : comme ça le spectro et la scalo ont la même taille
    %
    spectrogram = fullspectrogram(taillespectre:-1:1,:);
    spectrogram = spectrogram/max(max(spectrogram))*255;
    %
    spectrogramcontrasted = fullspectrogram(sizesig:-1:1,:).^contrast;
    spectrogramcontrasted = spectrogramcontrasted/max(max(spectrogramcontrasted))*255;
    %
    figure(5);
    clf;
    imagesc([1:sizesig],([taillespectre:-1:1]-1)/taillespectre/2,spectrogram);
    title([nomdusignal '- spectrogramme ; trames plus grandes'],'FontSize',15);
    colormap(mapmap);
    hold off;
    print('tftbs6_3.ps','-depsc2');
  end;


  %%% Morlet
  fullscalogrammorlet = lagis_mytfrscalo(sig3, 1:sizesig, sizemorlet, 0.0025, 0.4975, taillespectre, 1);
  %fullscalogrammorlet = lagis_mytfrscalo(sig3, 1:sizesig, sizemorlet, 0.05, 0.45, taillespectre, 1);
  scalogrammorlet = fullscalogrammorlet; %%% pas d'inversion
  scalogrammorlet = scalogrammorlet/max(max(scalogrammorlet))*255;
  figure(6);
  clf;
  imagesc(scalogrammorlet);
  title([nomdusignal '- scalogram de Morlet'],'FontSize',20);
  colormap(mapmap);
  hold off;
  print(nomfig2,'-depsc2');

  %%% scalogramme original
  fullscalogrammorlet2 = tfrscalo(sig3, 1:sizesig, sizemorlet, 0.0025, 0.4975, taillespectre, 1);
  scalogrammorlet2 = fullscalogrammorlet2; %%% pas d'inversion
  scalogrammorlet2 = scalogrammorlet2/max(max(scalogrammorlet2))*255;
  figure(7);
  clf;
  imagesc(scalogrammorlet2);
  title([nomdusignal '- scalogram de Morlet'],'FontSize',20);
  colormap(mapmap);
  hold off;
  print(nomfig2a,'-depsc2');
end;

