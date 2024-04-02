%%%%%%%%%%%%%%%%%%%%
% Test l'optimisation LOQO seule, cas "SVM 1 classe"
% Code Lagis et CentraleSupélec
%
% Stéphane Rossignol - 14/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%


%%%%%% On peut utiliser ce code sous deux formes : un script, ou une fonction

%%% utilisation en tant que 'script' : décommenter les 2 lignes ci-dessous
clear all;
close all;
%%% utilisation en tant que fonction (utilisé par "script_incre") : décommenter la ligne ci-dessous
%function [] = test_optimisation_loqo_svm1classe ()



%%% initialisations
mmm =   300;      % nombre de points dans la classe
ndimensions = 2;    % nombre de dimensions
%%%sig2 = 2;        % 'sig2' est estimé ou forcé un peu plus bas
nu   = 0.1;
use_octave = 1;     % utilisé pour LOQO (garder la valeur à 1)
calcul_loqo = 1;    % calcul aussi avec l'optimiseur LOQO (garder la valeur à 1)
fractionseuil=20.0; % tous les "alpha>max(alpha)/fractionseuil" sont gardés


%%% noyaux disponibles
RBFGAUSSIEN          = 0;  %%% le seul disponible pour CentraleSupélec pour le moment
QUADRATIQUERATIONNEL = 1;
POLYNOMIAL           = 2;
MULTIQUADRATIQUE     = 3;


%%% formation des données
numexemple=2;
switch numexemple
  case 1
    % données random
    randn("seed",1);
    rand("seed",1);
    data = randn(ndimensions,mmm);

  case 2
    % problème des oiseaux
    datatout1 = load('featuresoiseau1.txt');
    datatout2 = load('featuresoiseau2.txt');
    datatout3 = load('featuresoiseau3.txt');

    % normalisation des données
    moy1=mean([datatout1(:,1)' datatout2(:,1)' datatout3(:,1)']);
    moy2=mean([datatout1(:,2)' datatout2(:,2)' datatout3(:,2)']);
    std1=std([datatout1(:,1)' datatout2(:,1)' datatout3(:,1)']);
    std2=std([datatout1(:,2)' datatout2(:,2)' datatout3(:,2)']);

    datatout=datatout1; %%% on classifie les données du premier oiseau
    mmm=1000;           %%% on prend mmm données sur 1000
    data2=datatout(1:mmm,:)';
    data=data2;
    data(1,:)=data2(1,:)-moy1;
    data(1,:)=data(1,:)/std1;
    data(2,:)=data2(2,:)-moy2;
    data(2,:)=data(2,:)/std2;
end;


%%% on essaie d'estimer le 'sig2' à utiliser : on utilise ici le critère de Smith
sig2th=1;
for ii=1:ndimensions
   sig2th=sig2th*(max(data(ii,:))-min(data(ii,:)));
end;
sig2th=sig2th^(1/ndimensions);
sig2th=2*(sig2th^2);
sig2=sig2th;

%%% on peut forcer le 'sig2' à une certaine valeur, en décommentant cette ligne
%%% => par exemple, partir de 0.1 et augmenter jusqu'à 'sig2th' qui ici vaut environ 40
sig2=10;


%%% optimisation avec l'algorithme LOQO
if (calcul_loqo==1)
  %%% Attention : il faut vérifier que 'data' a la bonne forme
  %%% => c'est-à-dire, ne pas intervertir 'dimensions' et 'nombre de données d'apprentissage'

  %%%
  fprintf(1,'\nOn demarre le calcul des SVM avec LOQO\n');

  [wwxLOQO, bbbLOQO, alphaLOQO] = lagis_svm_1classe_noyau_nd (data, sig2, nu, use_octave, RBFGAUSSIEN);
end;


calcon=1; %%% calcul et plot de la fonction de décision obtenue
if calcon == 1
  %%%
  nzz=100;
  for ii=1:ndimensions
    minxx(ii) =min(data(ii,:));
    maxxx(ii) =max(data(ii,:));
    minxxz(ii)=minxx(ii)-1.0;
    maxxxz(ii)=maxxx(ii)+1.0;
  end;

  for ii=1:nzz
    xxt(1,1) = minxxz(1) + (maxxxz(1)-minxxz(1))/(nzz-1)*(ii-1);
    for jj=1:nzz
      xxt(2,1) = minxxz(2) + (maxxxz(2)-minxxz(2))/(nzz-1)*(jj-1);

      xxd1(ii,jj) = xxt(1,1);
      xxd2(ii,jj) = xxt(2,1);

      if (calcul_loqo==1)
        zzdLOQO(ii,jj) = alphaLOQO*lagis_rbf_gaussien(xxt, data, sig2)' + bbbLOQO;
      end;
    end;
  end;

  if (calcul_loqo==1)
    figure(2);
    clf;
    grid on;
    hold on;
    plot(data(1,:),data(2,:),'ob');
    contour(xxd1,xxd2,zzdLOQO,[0 0]);
    title('Frontiere obtenue avec optimisation LOQO');
    hold off;
  end;

  if (calcul_loqo==1)
    figure(4);
    clf;
    grid on;
    hold on;
    plot(alphaLOQO);
    title('alphas obtenus avec optimisation LOQO');
    hold off;
  end;
end;


%%% attention : cette boucle peut prendre du temps,
%%% quand le nombre mmm est grand (>40000) !!!!
nvecteursupport = sum(alphaLOQO>max(alphaLOQO)/fractionseuil);
fprintf(1,"\n\nsig2th=%f\n\nnombre de vecteurs supports : %d\n\n", sig2th,nvecteursupport);

test_asma=1;  %%% on regarde la répartition des échantillons, entre hors-classe, en-classe et vecteurs supports
if test_asma==1
   for ii=1:mmm
      fdecision(ii) = alphaLOQO*lagis_rbf_gaussien(data(:,ii), data, sig2)'  + bbbLOQO;
      if fdecision(ii)>0
         decision(ii) =  1;
      else
         decision(ii) = -1;
      end;
   end;

   nvecteursupport = sum(alphaLOQO>max(alphaLOQO)/fractionseuil);
   appartient    = sum(decision==1);
   appartientpas = sum(decision==-1);
   %fprintf(1,"\n\nsig2th=%f\n\nnombre de vecteurs supports : %d\n\n", sig2th, nvecteursupport);
   fprintf(1,"      %d echantillons appartiennent a la classe\n", appartient);
   fprintf(1,"      %d echantillons n appartiennent a la classe\n\n", appartientpas);

   figure(5);
   clf;
   subplot(211);
   grid on;
   hold on;
   plot(fdecision);
   title('LOQO - fonction de decision pour chaque echantillon');
   hold off;
   subplot(212);
   grid on;
   hold on;
   plot(decision);
   title('LOQO - decision pour chaque echantillon');
   hold off;
end;

