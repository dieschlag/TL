%%% Basis pursuit avec ``ltfat'' (sous octave/matlab, donc)
%%%
%%% St\'ephane Rossignol -- 2021

clear all;
close all;

pkg load ltfat


%%% décommenter une des lignes pour traiter tel ou tel signal
%[f,fs] = greasy;ns=1;
[f,fs]=audioread('parenthese.wav');ns=2;


%%% 'frame' (=dictionnaire) choisie

%% 'frame' de Gabor, avec une redondance de 8
%F = frame('dgtreal','gauss',64,512);
%F = frame('dgtreal','gauss',512,1024);

%% 'frame' de Gabor => ce qu'on a utilisé pour le 'matching pursuit'
%F = frame('dgt','hann',256,1024);
F = frame('dgtreal','hann',256,1024);


%%%
% Solveur du problème de basis pursuit
[c,~,~,frec,cdd] = franabp(F,f);


%%% Plotte les coefficients extraits
figure(1);
clf;
plotframe(F,c,'dynrange',50);
hold on;
title('Coefficients extraits');
hold off;


### tri des coefficients parcimonieux
[dd,ooo]=sort(abs(c),'descend');


%%% on ne garde que les 'maxit' plus grands coefficients
%%% => pour le 'matching pursuit', ça se faisait directement dans 'franamp'
%%% décommenter une de ces lignes
%maxit = 10;  %%% pas intelligible
%maxit = 50;  %%% pas très intelligible
maxit = 100; %%% c'est déjà intelligible !!!
%maxit = 500; %%% encore un peu 'bizarre'
%maxit = 1000; %%% encore un peu 'bizarre'
%maxit = 2000; %%% quasi nickel
%maxit = 4000; %%% nickel
%maxit = 6000; %%% nickel
%maxit = 7000; %%% nickel
%maxit = 7279; %%% nickel
%maxit = 7500; %%% nickel
%maxit = 8000; %%% nickel


cc=c;
cc(ooo(maxit+1:end))=0; %%% tous les plus petits coefficients sont mis à zéro
fprintf(1,'nombre de coefs > 0 : %d\n', sum(abs(cc)>0)); %%% vérification


### resynthése
frec = frsyn(F,cc);
frec=frec(1:length(f));
sound(frec,fs);


###
figure(3);
clf;
plot(f);
hold on;
plot(frec,'r');
title('signal original (bleu) et signal reconstruit (rouge)');
hold off;


%%% Erreur de reconstruction (devrait être proche de 0).
%%% 'frec' est obtenu en appliquant l'operateur de resynthèse de la
%%% 'frame' 'F' sur les coefficients extraits 'c'.
fprintf(1,'erreur de reconstruction : %f\n',norm(f-frec));


%%% Décroissance des coefficients triés en amplitudes
%%% (valeur absolue) (compressibilité des coefficients)
figure(4);
clf;
semilogx([sort(abs(c),'descend')/max(abs(c))]);
title('coefficients plottes en ordre decroissant');
hold on;
hold off;

