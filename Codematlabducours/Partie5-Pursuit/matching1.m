%%% Matching pursuit avec ``ltfat'' (sous octave/matlab, donc)
%%%
%%% Stéphane Rossignol -- 2021

clear all;
close all;

pkg load ltfat


[f,fs] = greasy;    %%% son de ltfat


meth=1;
if (meth==1)
  %%% 1. on ne prend que la TFCT avec Hanning
  F = frame('dgt','hann',256,1024); %%% dictionnaire utilisé par le matching pursuit
else
  %%% 2. on prend la TFCT plus des ondelettes
  F1 = frame('dgt','hann',256,1024);
  F2 = frame('fwt', 'sym20', 2, 'zero');  %%% attention : ne pas prendre J trop grand
  F = frame('fusion',1,F1,F2);
end;

maxit = 4000;       %%% nombre maximum d'itérations (~N dans les slides)


%%% Matching Pursuit (MP)

[c2,frec2,info2] = franamp(F,f,'mp','maxit',maxit);

figure(1);
clf;
plot(20*log10([info2.relres]));
title('erreur residuelle en dB');
hold off;
print -depsc2 match1.eps

figure(2);
clf;
plot([0:length(f)-1]/fs,f);
title('signal original');
hold off;
print -depsc2 match2.eps

figure(3);
clf;
plot([0:length(f)-1]/fs,real(frec2));
title('signal reconstruit');
hold off;
print -depsc2 match3.eps

figure(4);
clf;
plot([0:length(f)-1]/fs,f-real(frec2));
title('difference entre le signal et le signal reconstruit');
hold off;
print -depsc2 match4.eps

