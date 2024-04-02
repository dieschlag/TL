%%% St\'ephane Rossignol -- 2021

clear all;
close all;

addpath('tftb-0.2/mfiles');

%pkg load signal


%%% initialisations
mapmap = lagis_colormap_matlab1();
taillespectre=1024;

[sig1,fe] = audioread('parenthese.wav');
tsig=length(sig1);

taillewav=100;

fullscalogram = lagis_mytfrscalo(sig1, 1:tsig, taillewav, 0.01, 0.45, taillespectre, 1);
scalogram = fullscalogram; %%% pas d'inversion
figure(2);
clf;
imagesc(20*log10(scalogram));
title('parenthese - scalogram de Morlet','FontSize',20);
colormap(mapmap);
ylim([0 200]);
hold off;
print -depsc2 wavparent1.ps

fullscalogram2 = tfrscalo(sig1, 1:tsig, taillewav, 0.01, 0.45, taillespectre, 1);
scalogram2 = fullscalogram2; %%% pas d'inversion
figure(3);
clf;
imagesc(20*log10(scalogram2));
title('parenthese - scalo Morlet original','FontSize',20);
colormap(mapmap);
ylim([0 400]);
hold off;
print -depsc2 wavparent2.ps

h1 = blackman(641);
fullspectrogram = tfrsp(sig1, 1:tsig, 2*taillespectre, h1);
spectrogram = fullspectrogram(taillespectre:-1:1,:);
figure(4);
clf;
imagesc([1:tsig],([1:-1:taillespectre]-1)/taillespectre/2,20*log10(spectrogram(1024:-1:1,:)));
colormap(mapmap);
title('parenthese - spectrogramme','FontSize',20);
ylim([0 200]);
hold off;
print -depsc2 wavparent3.ps

