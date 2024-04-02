%%% SVM 2 classes -- 3 classes ; la stratégie 1 contre tous est utilisée (voir slide 222)
%%%
%%% Stéphane Rossignol -- 2021-2022


clear all; %%% toujours mettre cette ligne au début d'un script
close all; %%% toujours mettre cette ligne au début d'un script


%%%%% les données qui suivent (jusqu'au label 'fin des fausses données') sont des fausses données
%%%%% => pour le vrai problème des chants d'oiseaux, il faut modifier
%%%%%    les lignes suivantes (notamment, remplacer les '...stdd*randn(1,nn)'
%%%%%    par les features données dans les fichiers 'featuresoiseau1.txt' etc.), puis
%%%%%    par les features extraits par vous
%%%%% => attention : ne pas inverser vecteurs lignes et vecteurs colonnes !!!!!!!!!!
%%%%%    il faut vérifier que les vecteurs sont dans le bon sens (avec la fonction 'size',
%%%%%    par exemple)

nn=100;   %%% nombre de points par classe
          %%% => pour le problème des chants d'oiseaux, il faut prendre 'nn=1000'
stdd=0.1; %%% éparpillement des points de chaque classe

%%% classe 1
cl1_d1=stdd*randn(1,nn);
cl1_d2=stdd*randn(1,nn);

%%% classe 2
cl2_d1=1+stdd*randn(1,nn);
cl2_d2=0+stdd*randn(1,nn);

%%% classe 3
cl3_d1=sqrt(2)/2+stdd*randn(1,nn);
cl3_d2=sqrt(2)/2+stdd*randn(1,nn);

%%%%%% fin des fausses données


%%% normalisation des 2 features en moyenne et en écart-type
%%%%%% pour les chants d'oiseaux, il faut normaliser avec les 1000 points de 'SONS'
%%%%%% et il faut normaliser les features obtenus pour 'SONS-VC' avec les moy1...
%%%%%% obtenus ici (avec la base d'apprentissage)
moy1=mean([cl1_d1 cl2_d1 cl3_d1]);
moy2=mean([cl1_d2 cl2_d2 cl3_d2]);
std1=std([cl1_d1 cl2_d1 cl3_d1]);
std2=std([cl1_d2 cl2_d2 cl3_d2]);
cl1_d1=(cl1_d1-moy1)/std1;
cl1_d2=(cl1_d2-moy2)/std2;
cl2_d1=(cl2_d1-moy1)/std1;
cl2_d2=(cl2_d2-moy2)/std2;
cl3_d1=(cl3_d1-moy1)/std1;
cl3_d2=(cl3_d2-moy2)/std2;


%%%%%% attention : ici, pour le problème des chants d'oiseaux, il faut découper la base d'apprentissage
%%%%%% (obtenue à partir de 'SONS') en deux parties : 80 % pour faire l'apprentissage, et 20 % pour tester


%%% les 2 paramètres SVM à régler
sig2=0.1; %% 0.3  semble pas mal pour les 'fausses données'
nu=0.01;  %% 0.01 semble pas mal pour les 'fausses données'


figure(1);
clf;
grid on;
hold on;
plot(cl1_d1,cl1_d2,'or');
plot(cl2_d1,cl2_d2,'*b');
plot(cl3_d1,cl3_d2,'+g');
hold off;
drawnow;


%%%%%% normalement, vous n'avez rien à modifier ci-dessous

%%% calcul de la solution 1 : optimisation sous contraintes, méthode LOQO

%% 1er SVM : cl1 versus cl2+cl3

fprintf(1,'-----> classe1 versus classe2 + classe3\n');
nn1=nn;
nn2=nn+nn;
donnees1 = [[cl1_d1 cl2_d1 cl3_d1]' [cl1_d2 cl2_d2 cl3_d2]']'; %%% données
classes1 = [ones(1,nn1) -ones(1,nn2)];             %%% labels des 2 classes

[alphaloqo1,yloqo1,bbbopt1] = uncontretousoutil(donnees1,nn1,nn2,classes1,sig2,nu);
tmpv2_1 = alphaloqo1.*classes1';
fprintf(1,'pressez enter pour continuer\n');
pause;

%% 2ème SVM : cl2 versus cl1+cl3

fprintf(1,'-----> classe2 versus classe1 + classe3\n');
nn1=nn;
nn2=nn+nn;
donnees2 = [[cl2_d1 cl1_d1 cl3_d1]' [cl2_d2 cl1_d2 cl3_d2]']'; %%% données
classes2 = [ones(1,nn1) -ones(1,nn2)];             %%% labels des 2 classes

[alphaloqo2,yloqo2,bbbopt2] = uncontretousoutil(donnees2,nn1,nn2,classes2,sig2,nu);
tmpv2_2 = alphaloqo2.*classes2';
fprintf(1,'pressez enter pour continuer\n');
pause;

%% 3ème SVM : cl3 versus cl1+cl2

fprintf(1,'-----> classe3 versus classe1 + classe2\n');
nn1=nn;
nn2=nn+nn;
donnees3 = [[cl3_d1 cl1_d1 cl2_d1]' [cl3_d2 cl1_d2 cl2_d2]']'; %%% données
classes3 = [ones(1,nn1) -ones(1,nn2)];             %%% labels des 2 classes

[alphaloqo3,yloqo3,bbbopt3] = uncontretousoutil(donnees3,nn1,nn2,classes3,sig2,nu);
tmpv2_3 = alphaloqo3.*classes3';
fprintf(1,'pressez enter pour continuer\n');
pause;

%%%%%% à partir de là, vous allez avoir des choses à modifier pour le problème des
%%%%%% chants d'oiseaux


%%%%%% cette partie du code doit être modifiée pour le problème avec les chants d'oiseaux
%%%%%% => il n'y a pas de pavage à faire

%%% pavage du plan, pour montrer les frontières entre les 3 classes
%%% => je prends 'nzz' par 'nzz' points, couvrant tout l'espace à deux dimensions,
%%%    et je regarde à quelle classe chacun de ces points appartient
nzz=500;

mind1c=min(donnees1(1,:))-0.0;
maxd1c=max(donnees1(1,:))+0.0;
pasd1=(maxd1c-mind1c)/nzz;

mind2c=min(donnees1(2,:))-0.0;
maxd2c=max(donnees1(2,:))+0.0;
pasd2=(maxd2c-mind2c)/nzz;

ppd1res = [];
ppd2res = [];
ppd1=mind1c;

%%%%%% pour les chants d'oiseaux, ces deux boucles doivent être enlevées
%%%%%% à la place de 'ppd1' et 'ppd2', il faut mettre un vrai point à classer
%%%%%% c'est-à-dire les points pour obtenir : 
%%%%%%   Erreur d’apprentissage (ou risque empirique)         : les points d'apprentissage (80% de SONS)
%%%%%%   Erreur de généralisation (ou risque espéré, ou réel) : les points de test (20% de SONS)
%%%%%%   Erreur de validation croisée                         : les points de validation croisée (SONS-VC)
for ii=1:nzz
  ppd2=mind2c;
  for jj=1:nzz
     clnew1 = tmpv2_1'*lagis_rbf_gaussien([ppd1 ppd2]', donnees1, sig2)' + bbbopt1;
     clnew2 = tmpv2_2'*lagis_rbf_gaussien([ppd1 ppd2]', donnees2, sig2)' + bbbopt2;
     clnew3 = tmpv2_3'*lagis_rbf_gaussien([ppd1 ppd2]', donnees3, sig2)' + bbbopt3;

     %%% on détermine le maximum, pour trouver à quelle classe appartient le point [ppd1 ppd2]
     %%% voir slide 222 (Stratégie Un contre tous)
     if clnew1>clnew2 && clnew1>clnew3
       clnew(ii,jj)=1;   %%% classe 1 : pour les chants d'oiseaux, il faut déterminer si ça correspond
                         %%%    à l'oiseau 1, à l'oiseau 2 ou à l'oiseau 3
     elseif clnew2>clnew1 && clnew2>clnew3
       clnew(ii,jj)=0;   %%% classe 2 : pour les chants d'oiseaux, il faut déterminer si ça correspond
                         %%%    à l'oiseau 1, à l'oiseau 2 ou à l'oiseau 3
     elseif clnew3>clnew1 && clnew3>clnew2
       clnew(ii,jj)=-1;  %%% classe 3 : pour les chants d'oiseaux, il faut déterminer si ça correspond
                         %%%    à l'oiseau 1, à l'oiseau 2 ou à l'oiseau 3
     else
       fprintf(1,'probleme\n');
     end;

     if ii==1
       ppd2res=[ppd2res ppd2];
     end;
    ppd2=ppd2+pasd2;
  end;
  ppd1res=[ppd1res ppd1];
  ppd1=ppd1+pasd1;
end;

%%% plot pour la route
figure(7);
clf;
title('frontiere entre les 3 classes');
hold on;
imagesc(ppd1res,ppd2res,-clnew');
grid on;
hold on;
plot(cl1_d1,cl1_d2,'*r');
plot(cl2_d1,cl2_d2,'ob');
plot(cl3_d1,cl3_d2,'+g');
xlabel('dimension 1');
ylabel('dimension 2');
xlim([min(ppd1res) max(ppd1res)]);
ylim([min(ppd2res) max(ppd2res)]);
hold off;
colormap('cool');
drawnow;
