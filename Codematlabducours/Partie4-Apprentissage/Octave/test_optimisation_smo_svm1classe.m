%%%%%%%%%%%%%%%%%%%%
% Test l'optimisation SMO, cas "SVM 1 classe" -- test_optimisation_smo
% On peut comparer aussi à ce que donne l'optimiseur LOQO
% Code Lagis et CentraleSupélec
%
% TODO :
% ------
%
% - 14/04/06
%   déséquilibre entre SMO et LOQO pour l'appel à la fonction d'optimisation
%   => voir : "lagis_optimisation_smo" versus "lagis_svm_1classe_noyau_nd"
%      et de plus leurs "retours" sont différents
%
% - 27/04/06 - 31/08/06
%   Il y a un étrange "bbbSMO = -bbbSMOf;" à faire. Ca vient de quoi ?
%   Sans aucun doute des différences de notations
%
% Stéphane Rossignol - 14/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%


%%%%%% On peut utiliser ce code sous deux formes : un script, ou une fonction

%%% utilisation en tant que 'script' : décommenter les 2 lignes ci-dessous
clear all;
close all;
%%% utilisation en tant que fonction (utilisé par "script_incre") : décommenter la ligne ci-dessous
%function [] = test_optimisation_smo_svm1classe ()


%%% 2021 : obsolète
%path(DEFAULT_LOADPATH,'/home/rossigno/LAGIS-CODE/SMO-1-CLASSE');


%%% initialisations
mmm =   300;      % nombre de points dans la classe ; peut être chnagé plus loin
ndimensions = 2;    % nombre de dimensions
%%%sig2 = 2;        % 'sig2' est estimé ou forcé un peu plus bas
nu   = 0.1;
utilise_code_c = 0; % SMO en C (1) ou en octave (0) ; note : le code en C n'est pas fourni pour la ST7
use_octave = 1;     % utilisé pour LOQO (garder la valeur à 1)
calcul_loqo = 1;    % calcul aussi avec l'optimiseur LOQO (mettre à 1 pour ça)
fractionseuil=20.0; % tous les "alpha>max(alpha)/fractionseuil" sont gardés


%%% noyaux disponibles
RBFGAUSSIEN          = 0;  %%% le seul disponible pour CentraleSupélec pour le moment
QUADRATIQUERATIONNEL = 1;
POLYNOMIAL           = 2;
MULTIQUADRATIQUE     = 3;


%%% formation des données
numexemple=1;
switch numexemple
  case 1
    % données random
    randn("seed",1);
    rand("seed",1);
    data = randn(ndimensions,mmm);

  case 2
    % problème des oiseaux
    datatout = load('featuresoiseau1.txt');
    mmm=1000;
    data2=datatout(1:mmm,:)';
    data=data2;
    data(1,:)=data2(1,:)-mean(data2(1,:));
    data(1,:)=data(1,:)/std(data(1,:));
    data(2,:)=data2(2,:)-mean(data2(2,:));
    data(2,:)=data(2,:)/std(data(2,:));
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
sig2=100;


%%% optimisation avec l'algorithme LOQO
if (calcul_loqo==1)
  %%% Attention : il faut vérifier que 'data' a la bonne forme
  %%% => c'est-à-dire, ne pas intervertir 'dimensions' et 'nombre de données d'apprentissage'

  %%%
  fprintf(1,'\nOn demarre le calcul des SVM avec LOQO\n');

  [wwxLOQO, bbbLOQO, alphaLOQO] = lagis_svm_1classe_noyau_nd (data, sig2, nu, use_octave, RBFGAUSSIEN);
end;


%%% optimisation avec l'algorithme SMO
if utilise_code_c==0
  %%% SMO - code octave (lent)

  %%% Attention : il faut vérifier que 'data' a la bonne forme
  %%% => c'est-à-dire, ne pas intervertir 'dimensions' et 'nombre de données d'apprentissage'

  %%%
  fprintf(1,'\nOn demarre le calcul des SVM avec SMO\n');

  topt=clock;
  [alphaSMO, bbbSMOf] = lagis_smo_1classe(data, nu, sig2, topt,'');
  eopt=etime(clock,topt);
else
  %%% SMO - code C (plus rapide que le code octave)
  % fichier de données formé

  %versionpg = "";
  %versionpg = "-v2";
  %versionpg = "-v3";
  %versionpg = "-v4";

  %versionpg = "-v5";  % version stable à utiliser comme référence
  %% => -v6 est dédiée au calcul des modèles pour le projet du Lagis
  versionpg = "-v7";   % version stable la plus rapide (et de développement)
                       % ATTENTION : "lagis_smo_1classe-v7-incre.c" compile aussi sur "lagis_smo_1classe-v7" !!!!
                       % => voir "script_incre.m"
  %versionpg = "-v8";

  %%% séquentiel
  %versionpg = "-v9a";   % version "SMO séquentiel", version "a"
  %versionpg = "-v9c";   % version "SMO séquentiel", version "c"

  fprintf(1,"La version du programme C utilisee est : %s\n",versionpg);

  if strcmp(versionpg,"-v9c")~=1
    ff=fopen(["SMO-1-CLASSE/donnees_pour_smo-c" versionpg ".txt"],"wt");
    fprintf(ff,"%d ",ndimensions);
    fprintf(ff,"%d ",mmm);
    fprintf(ff,"%f ",nu);
    fprintf(ff,"%f ",sig2);
    for ii=1:ndimensions
      fprintf(ff,"%f ",data(ii,:));
    end;
    fclose(ff);
  else
    %%% ff = ensemble complet
    %%% fg = première moitié
    %%  fh = seconde moitié
    fg=fopen(["SMO-1-CLASSE/donnees_pour_smo-c" versionpg "-1.txt"],"wt");
    fprintf(fg,"%d ",ndimensions);
    fprintf(fg,"%d ",mmm/2);
    fprintf(fg,"%f ",nu);
    fprintf(fg,"%f ",sig2);
    for ii=1:ndimensions
      fprintf(fg,"%f ",data(ii,1:mmm/2));
    end;
    fclose(fg);

    fh=fopen(["SMO-1-CLASSE/donnees_pour_smo-c" versionpg "-2.txt"],"wt");
    fprintf(fh,"%d ",ndimensions);
    fprintf(fh,"%d ",mmm/2);
    fprintf(fh,"%f ",nu);
    fprintf(fh,"%f ",sig2);
    for ii=1:ndimensions
      fprintf(fh,"%f ",data(ii,mmm/2+1:mmm));
    end;
    fclose(fh);

    ff=fopen(["SMO-1-CLASSE/donnees_pour_smo-c" versionpg "-0.txt"],"wt");
    fprintf(ff,"%d ",ndimensions);
    fprintf(ff,"%d ",mmm);
    fprintf(ff,"%f ",nu);
    fprintf(ff,"%f ",sig2);
    for ii=1:ndimensions
      fprintf(ff,"%f ",data(ii,:));
    end;
    fclose(ff);
  end;

  % calcul
  if strcmp(versionpg,"-v9c")~=1
    commande=["SMO-1-CLASSE/lagis_smo_1classe" versionpg];
    topt=clock;  
    system(commande);
    eopt=etime(clock,topt);
  else
    commande1=["SMO-1-CLASSE/lagis_smo_1classe" versionpg " 1"];
    commande2=["SMO-1-CLASSE/lagis_smo_1classe" versionpg " 2"];
    commande0=["SMO-1-CLASSE/lagis_smo_1classe" versionpg " 0"];

    topt=clock;  
    system(commande1);
    system(commande2);
    system(commande0);
    eopt=etime(clock,topt);
  end;

  % fichier de résultat lu
  if strcmp(versionpg,"")==1
    load SMO-1-CLASSE/resultats_de_smo-c.txt;
    alphaSMO = resultats_de_smo_c(1:mmm)';
    bbbSMOf = resultats_de_smo_c(mmm+1);
  elseif strcmp(versionpg,"-v2")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v2.txt;
    alphaSMO = resultats_de_smo_c_v2(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v2(mmm+1);
  elseif strcmp(versionpg,"-v3")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v3.txt;
    alphaSMO = resultats_de_smo_c_v3(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v3(mmm+1);
  elseif strcmp(versionpg,"-v4")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v4.txt;
    alphaSMO = resultats_de_smo_c_v4(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v4(mmm+1);
  elseif strcmp(versionpg,"-v5")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v5.txt;
    alphaSMO = resultats_de_smo_c_v5(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v5(mmm+1);
  elseif strcmp(versionpg,"-v7")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v7.txt;
    alphaSMO = resultats_de_smo_c_v7(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v7(mmm+1);
  elseif strcmp(versionpg,"-v8")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v8.txt;
    alphaSMO = resultats_de_smo_c_v8(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v8(mmm+1);
  elseif strcmp(versionpg,"-v9a")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v9a.txt;
    alphaSMO = resultats_de_smo_c_v9a(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v9a(mmm+1);
  elseif strcmp(versionpg,"-v9c")==1
    load SMO-1-CLASSE/resultats_de_smo-c-v9c.txt;
    alphaSMO = resultats_de_smo_c_v9c(1:mmm)';
    bbbSMOf = resultats_de_smo_c_v9c(mmm+1);
  end;
end;
bbbSMO = -bbbSMOf;     %%%%% Inversion à cause des notations différentes

fprintf(1,'\ntemps requis pour optimiser les SVM sur %d points avec SMO : %f\n',mmm,eopt);


%%% pour "script_incre.m"
%ftmp=fopen("SMO-1-CLASSE/incre_tmp.txt","w");
%fprintf(ftmp,'%f',eopt);
%fclose(ftmp);

if (calcul_loqo==1)
   fprintf(1,'\nOn compare ce qu on obtient avec LOQO et ce qu on obtient avec SMO\n');
   fprintf(1,'bbbSMOf=%f bbbSMO=%f - bbbLOQO=%f\n',bbbSMOf,bbbSMO,bbbLOQO);
   fprintf(1,'sum(abs(alphaLOQO-alphaSMO))/mmm=%f   max(abs(alphaLOQO-alphaSMO))=%f\n',sum(abs(alphaLOQO-alphaSMO))/mmm,max(abs(alphaLOQO-alphaSMO)));
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

      zzdSMO(ii,jj)  = alphaSMO*lagis_rbf_gaussien(xxt, data, sig2)'  + bbbSMO;
      if (calcul_loqo==1)
        zzdLOQO(ii,jj) = alphaLOQO*lagis_rbf_gaussien(xxt, data, sig2)' + bbbLOQO;
      end;
    end;
  end;


  %%%
  figure(1);
  clf;
  grid on;
  hold on;
  plot(data(1,:),data(2,:),'ob');
  contour(xxd1,xxd2,zzdSMO,[0 0]);
  title('Frontiere obtenue avec optimisation SMO');
  hold off;

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
    figure(3);
    clf;
    subplot(211);
    grid on;
    hold on;
    plot(alphaSMO,'r');
    title('alphas obtenus avec optimisation SMO');
    hold off;
    subplot(212);
    grid on;
    hold on;
    plot(alphaLOQO,'k');
    title('alphas obtenus avec optimisation LOQO');
    hold off;

    figure(4);
    clf;
    grid on;
    hold on;
    plot(alphaLOQO);
    title('alphas obtenus avec optimisation LOQO');
    hold off;
  else
     figure(3);
     clf;
     grid on;
     hold on;
     plot(alphaSMO,'r');
     title('alphas obtenus avec optimisation SMO');
     hold off;
  end;
end;


%%% attention : cette boucle peut prendre du temps,
%%% quand le nombre mmm est grand (>40000) !!!!
nvecteursupport = sum(alphaSMO>max(alphaSMO)/fractionseuil);
fprintf(1,"\n\nsig2th=%f\n\nnombre de vecteurs supports : %d\n\n", sig2th,nvecteursupport);

test_asma=1;  %%% on regarde la répartition des échantillons, entre hors-classe, en-classe et vecteurs supports
if test_asma==1
   for ii=1:mmm
      fdecision(ii) = alphaSMO*lagis_rbf_gaussien(data(:,ii), data, sig2)'  + bbbSMO;
      if fdecision(ii)>0
         decision(ii) =  1;
      else
         decision(ii) = -1;
      end;
   end;

   nvecteursupport = sum(alphaSMO>max(alphaSMO)/fractionseuil);
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
   title('SMO - fonction de decision pour chaque echantillon');
   hold off;
   subplot(212);
   grid on;
   hold on;
   plot(decision);
   title('SMO - decision pour chaque echantillon');
   hold off;
end;

