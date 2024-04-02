%%% Entropie de l'exmple du cours (Notions d’entropie - Exemple)
%%%
%%% S. Rossignol -- 2021

clear all;
close all;


%%% probabilités des 'sons'
pri = [0.50  0.29  0.06  0.05  0.04  0.03  0.02  0.01];
cumpri = cumsum(pri);
fprintf(1,'sumpri=%f\n',sum(pri));

entrop=0;
for ii=1:8
  entrop = entrop - pri(ii)*log2(pri(ii));
end;

fprintf(1,'entropie : %f\n',entrop);


%%% longueurs des symboles
longueur=[1 2 3 4 5 6 7 7];

fprintf(1,'longueur moyenne symbole : %f\n',sum(pri.*longueur));


%%% un exemple
code1='0';
code2='10';
code3='110';
code4='1110';
code5='11110';
code6='111110';
code7='1111110';
code8='1111111';

phrase='';
phrasecodee='';

nsons=100000;  %%% nombre de 'sons' tirés au hasard (en respectant les probas)
mots=rand(1,nsons);
for ii=1:nsons
  if mots(ii)<=cumpri(1)
    phrase=[phrase 'a'];
    phrasecodee=[phrasecodee code1];
  elseif mots(ii)<=cumpri(2)
    phrase=[phrase 'b'];
    phrasecodee=[phrasecodee code2];
  elseif mots(ii)<=cumpri(3)
    phrase=[phrase 'c'];
    phrasecodee=[phrasecodee code3];
  elseif mots(ii)<=cumpri(4)
    phrase=[phrase 'd'];
    phrasecodee=[phrasecodee code4];
  elseif mots(ii)<=cumpri(5)
    phrase=[phrase 'e'];
    phrasecodee=[phrasecodee code5];
  elseif mots(ii)<=cumpri(6)
    phrase=[phrase 'f'];
    phrasecodee=[phrasecodee code6];
  elseif mots(ii)<=cumpri(7)
    phrase=[phrase 'g'];
    phrasecodee=[phrasecodee code7];
  else
    phrase=[phrase 'h'];
    phrasecodee=[phrasecodee code8];
  end;
end;

fprintf(1,'\nlongueur de la phrase : %d sons ; longueur du code : %d bits ; ratio : %f\n', nsons, length(phrasecodee), length(phrasecodee)/nsons);

