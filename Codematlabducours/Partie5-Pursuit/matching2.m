%%% Matching pursuit avec ``ltfat'' (sous octave/matlab, donc)
%%%
%%% St\'ephane Rossignol -- 2021

clear all;
close all;

pkg load ltfat


%%% décommenter une des lignes pour traiter tel ou tel signal
[f,fs] = greasy;ns=1;   %%% son de ltfat
%[f,fs]=audioread('parenthese.wav');ns=2;   %%% mon son


F = frame('dgt','hann',256,1024); %%% dictionnaire utilisé


%%% décommenter une de ces lignes
%maxit = 10;  %%% pas intelligible
%maxit = 50;  %%% presque intelligible
%maxit = 100; %%% c'est déjà intelligible !!!
maxit = 500; %%% encore un peu 'bizarre'
%maxit = 2000; %%% quasi nickel
%maxit = 4000; %%% nickel
%maxit = 6000; %%% nickel
%maxit = 7000; %%% nickel
%maxit = 7279; %%% nickel
%maxit = 7500; %%% nickel
%maxit = 8000; %%% nickel


%%% Matching Pursuit (MP)

[c2,frec2,info2] = franamp(F,f,'mp','maxit',maxit);

figure(1);
clf;
plot(20*log10([info2.relres]));
title('erreur residuelle en dB');
hold off;
print -depsc2 match1a.eps

figure(2);
clf;
plot([0:length(f)-1]/fs,f);
title('signal original');
hold off;
print -depsc2 match2a.eps

figure(3);
clf;
plot([0:length(f)-1]/fs,real(frec2));
title('signal reconstruit');
hold off;
print -depsc2 match3a.eps

figure(4);
clf;
plot([0:length(f)-1]/fs,f-real(frec2));
title('difference entre le signal et le signal reconstruit');
hold off;
print -depsc2 match4a.eps


if ns==1 && maxit==10
  figure(5);
  clf;
  plot(abs(c2));
  title('coefficients c2');
  hold off;
  print -depsc2 match5a.eps
end;


%%% écoute : son original et son reconstruit

pause(1);
sound(f,fs);
pause(1)
sound(frec2,fs);


if ns==2
  fprintf(1,'%d %d\n',maxit,sum(c2>0));
end;

