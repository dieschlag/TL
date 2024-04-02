%%% Algorithme de Hérault-Jutten (exemple du cours)
%%%
%%% Stéphane Rossignol -- 24/05/2011 ; 2021


%%%%%% on peut utiliser ce code en script (STANDALONE) ou en fonction (FUNCTION)
%%%%%% il faut choisir ci-dessous en décommentant/commentant les lignes qu'il faut

%%% STANDALONE
clear all;
close all;
nmaxiter=20000; %%% nombre max d'itérations

%%% FUNCTION
%function [succes,a12,a21] = jutten1(nmaxiter)



ecoute=1;       %%% on écoute les divers sons ou pas
nminiter=10000; %%% nombre minimum d'itérations


%%% les signaux sont construits

nkk=5000;     % nombre d'échantillons dans les signaux
kk=[0:nkk];

fe=8000;     % fréquence d'échantillonnage
Te=1/fe;     % période d'échantillonnage

u1=zeros(length(kk),1);
u2=zeros(length(kk),1);

f1=50 + 950*rand(1,1);   % fréquence fondamentale du premier signal harmonique
f2=50 + 950*rand(1,1);   % fréquence fondamentale du deuxième signal harmonique
npmax=min([ floor(fe/2./f1) floor(fe/2./f2) ]); % nombre maximum de partiels possible
npp=round(rand(1,1)*npmax);                     % nombre de partiels par signal
if (npp==0)
  npp=1;
end;


ampg=2; %%% pourquoi les performances dépendent-elles de ça ?
a1(1)=(0.5+rand(1,1)/2)*ampg;
a2(1)=(0.5+rand(1,1)/2)*ampg;
u1= a1(1)*sin(2.*pi*(f1)*kk*Te+2.*pi*rand(1,1));
u2= a2(1)*sin(2.*pi*(f2)*kk*Te+2.*pi*rand(1,1));
for np=2:npp
  a1(np)=rand(1,1)/np*ampg; % les amplitudes des partiels aussi sont tirées aléatoirement
  u1= a1(np)*sin(2.*pi*(np*f1)*kk*Te+2.*pi*rand(1,1));
  a2(np)=rand(1,1)/np*ampg; % les amplitudes des partiels aussi sont tirées aléatoirement
  u2= a2(np)*sin(2.*pi*(np*f2)*kk*Te+2.*pi*rand(1,1));
end;
u1=u1-mean(u1);
u2=u2-mean(u2);


%%% la matrice de mélange est tirée aléatoirement

a12=0.5*rand(1,1);
a21=0.5*rand(1,1);
y1 = 1.0*u1 + a12*u2;
y2 = a21*u1 + 1.0*u2;

pourson=max(abs([y1 y2]));

if ecoute==1
  fprintf(1,'on ecoute les signaux avant le melange\n');
  sound(u1/pourson,fe);pause(2);
  sound(u2/pourson,fe);pause(2);
end;

if ecoute==1
  fprintf(1,'on ecoute les signaux apres le melange\n');
  sound(y1/pourson,fe);pause(2);
  sound(y2/pourson,fe);pause(2);
end;


%%% conditions 1 de convergence : en pratique, non disponibles
%mean(u1.^4)*mean(u2.^4)
%9* mean( (u1.^2) .* (u2.^2) )^2


%%% main code
w12=zeros(nmaxiter,1); % initialisation
w21=zeros(nmaxiter,1); % initialisation
w12(1)=0.1;%(rand(1,1)-0.5)*0; %%% condition 2 de convergence : c'est plutôt w12*w21<1 ???
w21(1)=0.5;%(rand(1,1)-0.5)*0;

jj=1;
ll=0;
finish=0;
while finish==0
  jj=jj+1;

  cu1=y1;
  cu2=y2-w21(jj-1)*cu1;
  dcu1=1.;
  dcu2=1.;
  ii=0;
  while ( (dcu1>1e-5 || dcu2>1e-5) && ii<1000 )
    ii=ii+1;
    ocu1=cu1;cu1=y1-w12(jj-1)*cu2;dcu1=mean(abs(cu1-ocu1));
    ocu2=cu2;cu2=y2-w21(jj-1)*cu1;dcu2=mean(abs(cu2-ocu2));
  end;

  tau = 1.8e-2;
  dw12 = tau*mean( (cu1.^3).*cu2      );
  dw21 = tau*mean( cu1     .*(cu2.^3) );

  w12(jj)=w12(jj-1)+dw12;
  w21(jj)=w21(jj-1)+dw21;


  %%% cette condition d'arrêt n'est pas raisonnable, puisque 'u1' and 'u2'
  %%% ne sont pas connus
  aa1=sum((cu1-u1).^2)/sum(u1.^2) + sum((cu2-u2).^2)/sum(u2.^2); %%% cu1 correspond à u1 et cu2 à u2
  aa2=sum((cu1-u2).^2)/sum(u2.^2) + sum((cu2-u1).^2)/sum(u1.^2); %%% cu1 correspond à u2 et cu2 à u1
%  if (aa1<1.e-3 || aa2<1.e-3)
%    finish=1;
%  end;

  %%% c'est une condition d'arrêt plus raisonnable :
  %%% il faut que 'w12' et 'w21' ne varient plus beaucoup
  scond = abs(dw12)+abs(dw21);
  if ( scond < 1e-7 && jj>nminiter )
    finish=1;
  end;
  %%% ... ou qu'on ait ateint 'nmaxiter'
  if (jj==nmaxiter)
    finish=1;
  end;


  ll=ll+1;
  if (ll==500 || finish==1) %%% toutes les 500 itérations on écrit et plotte où on en est
    ll=0;
    fprintf(1, '(%d)\n', jj);
    fprintf(1, 'amplitudes:\n');
    fprintf(1, '%f  ', a1);
    fprintf(1, '\n');
    fprintf(1, '%f  ', a2);
    fprintf(1, '\n');
    fprintf(1, 'A trouver:  a12=%.10f  a21=%.10f     (f1=%f f2=%f)\n', a12, a21, f1, f2);
    fprintf(1, 'Trouves :   w12=%.10f  w21=%.10f\n', w12(jj), w21(jj));
    fprintf(1, 'condition d arret : %e\n', scond);

    figure(1);
    clf;
    plot(w12(1:jj));
    hold on;
    plot(w21(1:jj));
    plot([0 jj],[a12 a12],'r');
    plot([0 jj],[a21 a21],'r');
    title('evolution des estimations des parametres de melange');
    hold off;
    drawnow;

    fprintf(1, 'erreurs (1/RSB lineaire) entre signaux a trouver et signaux estimes : \ncu1=u1: %f ou cu1=u2: %f\ncu2=u1: %f ou cu2=u2: %f\n\n', 
       sum((cu1-u1).^2)/sum(u1.^2), sum((cu1-u2).^2)/sum(u2.^2), sum((cu2-u1).^2)/sum(u1.^2), sum((cu2-u2).^2)/sum(u2.^2));
    fflush(1);
  end;

end;


%%% on teste si ça a marché (est-ce que "ui" et "cui" sont presque pareils)
%%% note : 'aa1' et 'aa2' sont calculés plus haut

succes=0;
if (aa1<1.e-2) %%% on autorise 1% d'erreur (soit un RSB de 20 dB)
  succes=1;
end;
if (aa2<1.e-2) %%% on autorise 1% d'erreur (soit un RSB de 20 dB)
  succes=1;
end;

if succes==1
  tailleframe=0.02;
  mini=min([u1 u2 y1 y2 cu1 cu2]);
  maxi=max([u1 u2 y1 y2 cu1 cu2]);

  pourslides=0; %%% pour mes slides et le cours

  figure(1);
  clf;
  plot(u1(1:fe*tailleframe));
  title('premier signal avant melange','Fontsize',20);
  xlim([1 fe*tailleframe]);
  ylim([mini maxi]);
  hold off;
  if pourslides==1
    audiowrite('juttenson1.wav',u1,fe);
    print juttenfig1.eps -depsc2
  end;

  figure(2);
  clf;
  plot(u2(1:fe*tailleframe));
  title('second signal avant melange','Fontsize',20);
  xlim([1 fe*tailleframe]);
  ylim([mini maxi]);
  hold off;
  if pourslides==1
    audiowrite('juttenson2.wav',u2,fe);
    print juttenfig2.eps -depsc2
  end;

  figure(3);
  clf;
  plot(y1(1:fe*tailleframe));
  title('premier melange','Fontsize',20);
  xlim([1 fe*tailleframe]);
  ylim([mini maxi]);
  hold off;
  if pourslides==1
    audiowrite('juttenson3.wav',y1,fe);
    print juttenfig3.eps -depsc2
  end;

  figure(4);
  clf;
  plot(y2(1:fe*tailleframe));
  title('second melange','Fontsize',20);
  xlim([1 fe*tailleframe]);
  ylim([mini maxi]);
  hold off;
  if pourslides==1
    audiowrite('juttenson4.wav',y2,fe);
    print juttenfig4.eps -depsc2
  end;

  figure(5);
  clf;
  plot(cu1(1:fe*tailleframe));
  title('premier signal apres demelange','Fontsize',20);
  xlim([1 fe*tailleframe]);
  ylim([mini maxi]);
  hold off;
  if pourslides==1
    audiowrite('juttenson5.wav',cu1,fe);
    print juttenfig5.eps -depsc2
  end;

  figure(6);
  clf;
  plot(cu2(1:fe*tailleframe));
  title('second signal apres demelange','Fontsize',20);
  xlim([1 fe*tailleframe]);
  ylim([mini maxi]);
  hold off;
  if pourslides==1
    audiowrite('juttenson6.wav',cu2,fe);
    print juttenfig6.eps -depsc2
  end;

  if ecoute==1
    fprintf(1,'on ecoute les signaux apres le demelange\n');
    sound(cu1/pourson,fe);pause(2);
    sound(cu2/pourson,fe);pause(2);
  end;

  %pause;
end;

