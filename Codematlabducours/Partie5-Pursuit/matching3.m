%%% Matching pursuit avec ``ltfat'' (sous octave/matlab, donc)
%%%
%%% Stéphane Rossignol -- 2021

clear all;
close all;

pkg load ltfat


%%% décommenter une des lignes pour traiter tel ou tel signal
%[f,fs] = greasy;ns=1;
[f,fs]=audioread('parenthese.wav');ns=2;

F = frame('dgt','hann',256,1024);


%%% décommenter une de ces lignes
%maxit = 10;  %%% pas intelligible
%maxit = 50;  %%% presque intelligible
%maxit = 100; %%% c'est déjà intelligible !!!
%maxit = 500; %%% encore un peu 'bizarre'
%maxit = 2000; %%% quasi nickel
maxit = 4000; %%% nickel
%maxit = 6000; %%% nickel
%maxit = 7000; %%% nickel
%maxit = 7279; %%% nickel
%maxit = 7500; %%% nickel
%maxit = 8000; %%% nickel


%%% Matching Pursuit (MP)

[c2,frec2,info2] = franamp(F,f,'mp','maxit',maxit);


%%% Orthogonal Matching Pursuit (OMP)

[c1,frec1,info1] = franamp(F,f,'omp','cg','maxit',maxit);


figure(1);
clf;
plot(20*log10([info1.relres]));
title('erreur residuelle en dB');
hold off;
print -depsc2 match1omp.eps

figure(2);
clf;
plot([0:length(f)-1]/fs,f);
title('signal original');
hold off;
print -depsc2 match2omp.eps

figure(3);
clf;
plot([0:length(f)-1]/fs,real(frec1));
title('signal reconstruit');
hold off;
print -depsc2 match3omp.eps

figure(4);
clf;
plot([0:length(f)-1]/fs,f-real(frec1));
title('difference entre le signal et le signal reconstruit');
hold off;
if ns==2 && maxit==4000
  print -depsc2 match4omp.eps
end;


figure(5);
clf;
plot(20*log10(info1.relres),'b');
hold on;
plot(20*log10(info2.relres),'r');
legend({'OMP','MP'});
hold off;
if ns==2 && maxit==4000
  print -depsc2 match5omp.eps
end;


%%% écoute : son original et son reconstruit

pause(1);
sound(f,fs);
pause(1)
sound(frec2,fs);


if ns==2
  fprintf(1,'%d %d\n',maxit,sum(c1>0));
end;

