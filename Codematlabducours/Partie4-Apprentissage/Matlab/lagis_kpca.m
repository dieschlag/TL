%%%%%%%%%%%%%%%%%%%%
% K-PCA
% Code Lagis
%
% Le noyau utilisé ici est :
%    - le RBF gaussien de sigma sqrt(sig2/2)
%    - le polynomial de degré sig2
% => le choix du noyau peut-être décisif : voir la différence entre
%    test_pcakpca_ecole.m (RBFgaussien meilleur) et 
%    test_mfcc_pcakpca.m  (polynomial meilleur)
%
% Arguments :
% -----------
% 
% - "data_center" doit compter 
%     * 'n_dimensions'        lignes et 
%     * 'taille_echantillons' colonnes
%     => par exemple, size(data_center)=> (2, 300)
%   et chaque ligne doit être centrée
% - "sig2" est le paramètre "sig" pour "RBFgaussien" 
%           et le paramètre "p"   pour "polynomial"
%    il est calculé séparément (lagis_sig2 pour RBFgaussien)
%    mais ce n'est pas encore très efficace, alors tâtonnez !!!
% - "nomnoyau" est égal à "RBFgaussien" ou "polynomial"
% - on est sous "octave" ou pas
%
% liste TODO
% ----------
%
% 1) 24/01/2006
%    Des extensions à d'autres noyaux sont à implémenter.
%
% 2) 25/01/2006
%    La normalisation des vecteurs propres manque =>
%    mais les résultats obtenus sont les mêmes qu'avec LS_SVM, donc ???
%
% 3) 25/01/2006
%    Essayer de ne pas utiliser les matrices KKmean et m1m,
%    afin que de pouvoir utiliser des plus grands vecteurs de données
%
%
% Stéphane Rossignol - 24/01/2006 ; et 2021
%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
% Temps de calcul pour les "large-scale K-PCA"
% Exemple d'école "test_pcakpca_ecole.m"
%
% - il y a 3 classes
% - colonne 1 : complete et sans GM ; colonnes suivantes : neig=3
% - machine   : 193.51.52.10
%
% nbre échantillons par classe : temps en secondes (pas exactement CPU)
% nnn =  100 :     1.18 /     0.51 /    0.51(GM) /   0.35(GM+) / --(K)
% nnn =  200 :    10.52 /     2.82 /    2.82(GM) /   0.93(GM+) / --(K)
% nnn =  400 :    79.77 /    20.11 /   19.84(GM) /   4.00(GM+) / --(K)
% nnn =  800 :   614.90 /   154.97 /  152.67(GM) /  25.45(GM+) / --(K)
% nnn = 1600 : pas fait /  1213.61 / 1197.79(GM) / 183.18(GM+) / --(K)
% nnn = 2000 : pas fait / pas fait /    pas fait / 343.03(GM+) / --(K)
% nnn = 2200 : pas fait / pas fait /    pas fait / 456.42(GM+) / --(K)
% nnn = 2300 : pas fait / pas fait /    pas fait / 519.29(GM+) / --(K)
% nnn = 2400 : pas fait / pas fait /    pas fait / 588.70(GM+) / --(K)
% nnn = 2500 : pas fait / pas fait /    pas fait / 658.79(GM+) / --(K)
% nnn = 2550 : pas fait / pas fait /    pas fait /    voir (1) /  9092.89(K)
%              8976.92(K+) [21/02/06]
% nnn = 2600 : pas fait / pas fait /    pas fait /    pas fait /  9583.28(K) [15/02/06]
% nnn = 2700 : pas fait / pas fait /    pas fait /    pas fait / 10579.02(K) [16/02/06]
% nnn = 2800 : pas fait / pas fait /    pas fait /    pas fait / 11770.87(K) [16/02/06]
% nnn = 2900 : pas fait / pas fait /    pas fait /    pas fait / 13135.16(K) [16/02/06]
% nnn = 3000 : pas fait / pas fait /    pas fait /    pas fait / 14682.68(K) [20/02/06]
%              14371.16(K+) [21/02/06]
% nnn = 3100 : pas fait / pas fait /    pas fait /    pas fait /    pas fait
% nnn = 3200 : pas fait / voir (1) /    voir (1) /    pas fait /    pas fait
% nnn = 3400 : pas fait / pas fait /    pas fait /    pas fait /    pas fait
%              20998.43(K+) [21/02/06]
% nnn = 3500 : pas fait / pas fait /    pas fait /    pas fait /    pas fait
%              22845.46(K+) [23/02/06]
% nnn = 3600 : pas fait / pas fait /    pas fait /    pas fait /    pas fait
%              24811.90(K+) [24/02/06]
% nnn = 3700 : pas fait / pas fait /    pas fait /    pas fait /    pas fait
%              voir (1)
% nnn = 4000 : pas fait / pas fait /    pas fait /    pas fait /    pas fait
%              voir (1)
% nnn = 4600 : pas fait / pas fait /    pas fait /    pas fait /    pas fait
% nnn = 5000 : pas fait / pas fait /    pas fait /    pas fait / voir (2)
% nnn = 9000 : pas fait / pas fait /    pas fait /    pas fait / voir (2)
% nnn =10000 : pas fait / pas fait /    pas fait /    pas fait / voir (2)
%
% (1) out of memory => du coup, les aménagements :
%     - (GM),  modifications "grosses matrices" effectuées dans ce fichier
%     - (GM+), en étendant les 'clear' à "lagis_pca_ecole.m" aussi
%     - (K),   pour les grosses matrices KKmean est calculée non
%              matriciellement
%
% (2) Maximum variable size allowed by the program is exceeded.
%     KK=zeros(taille,taille);
%
% Remarque 07/02/2006 :
% ---------------------
%
% - pour nnn=2000, j'ai obtenu :
%   * temps requis pour calculer la KKmean : 320.80 (91.6 %)
%   * temps requis pour calculer la K-PCA :  350.26
% - pour nnn=2500, j'ai obtenu :
%   * temps requis pour calculer la KKmean : 613.22
%   * temps requis pour calculer la K-PCA :  658.79 (93.1 %)
% => donc, il faut travailler ce calcul de la KKmean !!!
%
%%%%%%%%%%%%%%%%%%%%

function [eig_vec_kpca,abs_val_eig] = lagis_kpca (data_center, sig2, nomnoyau, use_octave)


%%% K-PCA
% initialisations
fprintf(1,'Formation de quelques matrices\n');
taille      = size(data_center,2);
ndimensions = size(data_center,1);


%%% K-PCA
% calcul de la matrice de Gram
KK=zeros(taille,taille);
if strcmp(nomnoyau,'RBFgaussien')==1
   % noyau/kernel : exp(-||xi - xj||^2/sig2)
   for ii=1:taille
     for jj=1:taille
       phiphi = (data_center(1,ii)-data_center(1,jj))^2;
       for kk=2:ndimensions
         phiphi = phiphi + (data_center(kk,ii)-data_center(kk,jj))^2;
       end;
       KK(ii,jj) =  exp(-phiphi/sig2);
     end;
   end;
elseif strcmp(nomnoyau,'polynomial')==1
   % noyau/kernel : (xi.xj+1)^2
   for ii=1:taille
     for jj=1:taille
       phiphi = 1.0 + data_center(1,ii)*data_center(1,jj);
       for kk=2:ndimensions
         phiphi = phiphi + data_center(kk,ii)*data_center(kk,jj);
       end;
       KK(ii,jj) =  phiphi^sig2;
     end;
   end;
else
   fprintf(1,'\nNoyau indefini : doit etre RBF_gaussien ou polynomial\n');
end;


%%% K-PCA
% calcul KKmean
fprintf(1,'Calcul de KKmean\n');

limitetaille=7500; % c'est la valeur de 'limitetaille' qui convient
%limitetaille=200;

%%%tt1=cputime;
if taille<=limitetaille     %--test--%
  %%% K-PCA
  % avec formation de la matrice 1m
  m1m = ones(taille,taille)/(taille);

  KKmean = KK - m1m*KK - KK*m1m + m1m*KK*m1m;     %--test--%
%  KKmean1 = KK - m1m*KK - KK*m1m + m1m*KK*m1m;    %--test--%

  %%% K-PCA - (GM) grosses matrices
  % on n'a plus besoin de cette matrice
  clear m1m;
else    %--test--%
  %%% K-PCA
  % sans la formation de la matrice m1m

  %%% K-PCA
  %
  % A ENLEVER QUAND (K+) SERA OPERATIONNEL
  %
  % teste si KKmean peut-être formée
  % c'est de la perte de temps, mais ça permet d'être sûr
  % => si ce n'est pas possible, Matlab/Octave s'échappe 
%  KKmean=KK;
%  clear KKmean;

  %%% K-PCA
  % calculs intermédiaires
  ccc = 0;
  for kk=1:taille
    ccc = ccc+sum(KK(kk,:))/taille;
  end;
  ccc = ccc/taille;

  %%% K-PCA
  % boucles
  fprintf(1,'\n');
  tt2=clock;
  for ii=1:taille
    aaa = sum(KK(ii,:))/taille - ccc;

    % 0)
    % Le but de ces modifications est avant tout de permettre une "taille" aussi
    % grande que possible, avant que d'essayer de faire gagner du temps

    % 1)
    % On forme KKmean au fur et à mesure : KKmean = KK-aaa fait à
    % cet endroit du programme et pas dans la boucle jj ne ferait peut-être
    % pas gagner de temps car alors KKmean serait formée d'un coup
    % => tel qu'il est fait, le programme ralentit très progressivement

    % 2)
    % De la même façon, sum(KK(:,jj))/taille pourrait n'être calculé
    % qu'une seule fois, mais alors il faudrait un vecteur supplémentaire

    % 3)
    % Malheureusement, à cause de sum(KK(:,jj))/taille, KK ne peut pas
    % être effacée au fur et à mesure, suivant ii

    for jj=1:taille
      %KKmean(ii,jj) = KK(ii,jj) - sum(KK(:,jj))/taille - aaa;

      %%% (K+) : chaque colonne de KKmean est sauvée dans un fichier
      %%% ce qui permet de former KKmean après avoir désalloué KK
      KKtmp(jj) = KK(ii,jj) - sum(KK(:,jj))/taille - aaa;
    end;

    %%% Chaque colonne est sauvée dans un fichier
    nomcolonne=['LOGFILES/colonne' num2str(ii)];
    fid=fopen(nomcolonne,'w');
    fwrite(fid,KKtmp,'double');
    fclose(fid);

    %%%
    % tout ça, ça fait perdre du temps...
    % mais c'est pour vérifier où l'on en est dans le calcul...
    % celui qui attend doit prendre patience, et doit être aidé à
    % prendre patience
    ee2=etime(clock,tt2);
    te=ee2/ii*taille;
    fprintf(1,'temps restant %f sur %f (%f %%)\n',te-ee2,te,(te-ee2)/te*100);
    tsave=clock;
    fid = fopen('LOGFILES/exp.txt','wt');
    fprintf(fid,'temps restant %f sur %f (%f %%) [%d.%d]\n',te-ee2,te,(te-ee2)/te*100,tsave(4),tsave(5));
    fclose(fid);
  end;
end;    %--test--%
%%%ee1=cputime-tt1;
%%%fprintf(1,'temps requis pour calculer la KKmean : %f\n',ee1);


%%% K-PCA - grosses matrices, toujours
% Pour tester
% Les lignes %--test--% en commentaire sont à décommentariser
% et inversement
%KKmean
%pause;
%KKmean1
%diffmat = abs(KKmean1-KKmean)./abs(KKmean1);    %--test--%
%fprintf(1,'\n moyenne:%e max:%e min:%e\n',mean(mean(diffmat)),max(mean(diffmat)),min(mean(diffmat)));
    %--test--%


%%% K-PCA - (GM) grosses matrices
% on n'a plus besoin de cette matrice
clear KK;


%%% K-PCA - (K+) les colonnes de KKmean sont relues
%%% et KKmean est formée : ceci après que KK soit désallouée,
%%% ce qui permet d'atteindre des "taille"s supérieures
%%% Note : KKmean est symétrique => KKmean(ii,jj)=KKmean(jj,ii)
if taille>limitetaille
  fprintf(1,'\n\nFormation de KKmean\n\n');
  KKmean=[];
  tt2=clock;
  for ii=1:taille
    nomcolonne=['LOGFILES/colonne' num2str(ii)];
    fid=fopen(nomcolonne,'r');
    KKtmp=fread(fid,'double');
    KKmean=[KKmean KKtmp];
    fclose(fid);

    ee2=etime(clock,tt2);
    te=ee2/ii*taille;
    fprintf(1,'temps restant %f sur %f (%f %%)\n',te-ee2,te,(te-ee2)/te*100);
    tsave=clock;
    fid = fopen('LOGFILES/exp.txt','wt');
    fprintf(fid,'temps restant %f sur %f (%f %%) [%d.%d]\n',te-ee2,te,(te-ee2)/te*100,tsave(4),tsave(5));
    fclose(fid);
  end;
  clear KKtmp;
end;


%%% K-PCA - (GM) grosses matrices
% pour essayer avec des grosses matrices...
% pack %-> regarder comment ça marche et ce que ça fait exactement
% A FAIRE


%%% K-PCA
% valeurs propres et vecteurs propres
if taille>limitetaille
  fprintf(1,'KKmean est obtenue ; on entame le calcul des valeurs et vecteurs propres\n');
end;


complete=0;
if use_octave==1
  complete=1; % "eigs" n'existe pas sous "octave" (en 2006, en tout cas, c'était le cas)
end;
if complete==1
  neig=taille;
  [eig_vec_kpca,eig_val_kpca] = eig(KKmean);
else
  neig=24;
  [eig_vec_kpca,eig_val_kpca] = eigs(KKmean,neig);
end;


%%% K-PCA - quelques libérations de mémoire
prr=0;
if prr==0
  clear KKmean;
end;


%%% K-PCA
% réarrangement des valeurs propres
for ii=1:neig
  abs_val_eig(ii) = abs(eig_val_kpca(ii,ii));
end;
[abs_val_eig_sort,pos] = sort(abs_val_eig,'descend');

%%%
faitfig=0;
if faitfig==1
  figure(1)
  if use_octave==0
    clf;
  else
    %clearplot; %%% en 2006, clf sous octave n'existait pas
    clf;
  end;
  grid on;
  hold on;
  plot(abs_val_eig);
  title('24 plus grandes valeurs propres -- K-PCA');
  %zoom on; % "octave" se plaint
  hold off;
  print -dpsc2 figures/donnees_lanni_valeurspropreskpca.ps
end;

%fprintf(1,"Valeurs propres classiques :\n");
%fprintf(1,"\t%f\n",abs_val_eig_sort);


%%% Ajout du 25/09/06
% PRR (méthode de Padé-Rayleigh-Ritz)
% => au sujet de cette méthode, voir plutôt "lagis_kpca_prr.m"
if prr==1
  if taille>5
    taille_proj=5;
  else
    taille_proj=taille-2;
  end;

  minerreur=1e10;
  for ii=1:10
    [eig_vec_prr, eig_val_prr] = lagis_propres_prr(KKmean,taille_proj);
  
    %%% ce n'est pas la méthode standard pour valider Padé-Rayleigh-Ritz
    for ii=1:taille_proj
      eig_val_prr2(ii) = abs(eig_val_prr(ii,ii));
    end;
    [eig_val_prrsort,pos] = sort(eig_val_prr2,'descend');
    erreur = sum(abs(abs_val_eig_sort(1:taille_proj)-eig_val_prrsort))/sum(abs_val_eig_sort(1:taille_proj));

    if erreur<minerreur
      minerreur = erreur;
      fprintf(1,"\nValeurs propres Padé-Rayleigh-Ritz (erreur=%f):\n",minerreur);
      fprintf(1,"%.3f ",eig_val_prrsort);

      fprintf(1,"\nAttention : on considère les résultats de Padé-Rayleigh-Ritz\n");
      eig_vec_kpca = eig_vec_prr;
      abs_val_eig = eig_val_prr2;
    end;
  end;

  fprintf(1,"\nValeurs propres classiques :\n");
  fprintf(1,"%.3f ",abs_val_eig_sort(1:taille_proj));
end;

