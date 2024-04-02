%%%%%%%%%%%%%%%%%%%%
% Qualité des mesures de dissimilarité
% Code Lagis et CentraleSupélec
%
% Une partie du code vient de "test_optimisation_smo.m"
% Une autre vient de "lagis_evaluation.m" (code Lagis indisponible pour CentraleSupélec)
%
%
% NOTES :
% -------
%
%   - voir "lagis_qualite_mesures.m" pour le cas à une classe
%     original ; bien évidemment, le cas à une classe est
%     reproductible en mettant "nclasses = 1" (et tous les autres
%     paramètres comme dans "lagis_qualite_mesures.m")
%
%
% TODO :
% ------
%
% - 23/10/2006
%   Il y a des TODO ci-dessous
%
% Stéphane Rossignol - 07/09/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%


clear all;
close all;


%%% chemin
%path(DEFAULT_LOADPATH,'/home/rossigno/LAGIS-CODE/SMO-1-CLASSE');
%
% sur le portable (au Lagis)
%addpath '/home/rossigno/LAGIS-CODE/SMO-1-CLASSE';


%%% initialisations
nclasses = 2;      % nombre de classes
                   % pour le moment :
                   % - s'il y a 1 classe,  on reproduit les exemples de "lagis_qualite_mesures.m"
                   % - s'il y a 2 classes, la première varie selon "nchoix" et la seconde est fixe
                   % - s'il y a plus de 2 classes : implémentations à faire, donc a priori
                   %   pas la peine d'essayer pour le moment
mmmmax   = 2000;   % nombre de points par classe maximum
                   % (le nombre "mmm" est possiblement modifié plus loin)
ndimensions = 27;  % nombre de dimensions : si l'on veut visualiser le contour,
                   % il semble difficile de prendre autre chose que 2 ;
                   % sinon, pas de problème, sauf peut-être quand je considérerai autre
                   % chose que la "config. 1" ci-dessous
nchoix = 1;        % pour chacun des tests, le nombre de possibilité prises en compte
                   % attention : utiliser plus que "nchoix=4" implique des modifications
                   % plus ou moins importantes du code
preambule = "figures/"; % où on sauve les figures


%%% caractéristiques des différentes classes
%% config. 1 : les nuages sont égrenées le long de l'axe des x (dimension 1)
%%             et tout le reste est identique
moyenne = zeros(ndimensions, nclasses);
for cc=1:nclasses
  moyenne(1,cc) = 3*(cc-1);
end;


%%% autres initialisations
utilise_code_c = 0;   % code C (non fourni pour la ST7 de CentraleSupélec)
use_octave = 1;       % pour LOQO
fractionseuil=20.0;   % tous les "alpha>max(alpha)/fractionseuil" sont gardés
matid = zeros(ndimensions, ndimensions);   % pour la prise en compte des moyennes
for ii=1:ndimensions
  matid(ii,ii)=1.0;
end;


%%% noyau
RBFGAUSSIEN          = 0;  % le seul disponible à CentraleSupélec
QUADRATIQUERATIONNEL = 1;
POLYNOMIAL           = 2;
MULTIQUADRATIQUE     = 3;


%%% données
randn("seed",1);
rand("seed",1);
datatout=[];
for cc=1:nclasses
  matmoy = matid;
  for ii=1:ndimensions
    matmoy(ii,ii)=moyenne(ii,cc);
  end;
  datatout = [datatout randn(ndimensions,mmmmax)+(ones(mmmmax,ndimensions)*matmoy)'];
end;


%%% sig2 estimé -- pas utilisé => à refaire dans le cas >1 classe
%%sig2th=1;
%%for ii=1:ndimensions
%%   sig2th=sig2th*(max(datatout(ii,:))-min(datatout(ii,:)));
%%end;
%%sig2th=sig2th^(1/ndimensions);
%%sig2th=2*(sig2th^2);
%%sig2=sig2th;

%%% sig2 forcé
sig2=10.0;



%%%%%% Pour choisir un test (une 'simu'), décommentariser une des
%%%%%% "boucle"s ci-dessous ET commentariser les autres
%%%%%% => malheureusement, avec un "switch", ça ne 
%%%%%%    marche pas et il n'y a pas de "#ifdef" sous octave :-(

%%% 27/10/06 : quand "ndimensions=27", alors si "sig2=80"
%%% pour une classe et "sig2=1" pour l'autre, les mesures ne
%%% marchent plus
%%% -> phénomène moins facile à montrer quand on est
%%%    en 2 dimensions (prendre 0.01 et 2.0) que quand 
%%%    on est en 27 dimensions


%%% boucle sur "sig2"
simu=1;
faitcontour=1;

%% dim=2
%sig2res=[1 4 8 25];   % pour le cas où "nchoix=4"
%sig2res = 0.01;       % pour le cas où "nchoix=1", ce qui permet d'avoir une 
                       % figure 6 correcte, et les autres lisibles
%% dim=27
%sig2res =  1.0;
sig2res = 80.0;

mmm(1) = 100; %%% nombre de points par classe
mmm(2) = 100; %%% nombre de points par classe
nu = 0.1;
for nn=1:nchoix  %%% la fin de cette boucle intervient très bas
  %%% les paramètres de la seconde classe sont fixes : c'est la première qui varie
  sig2(1)=sig2res(nn);
  %sig2(1)=2;
  %sig2(1)
  if nclasses>1
    %% dim=2
    %sig2(2)=0.001; % avec une aussi petite valeur, ça explose même en dimension 2
    %sig2(2) = 2.0;

    %% dim=27
    sig2(2)=80.0;
  end;


%%% boucle sur "mmm"
%simu=2;
%faitcontour=1;
%mmmres=[2 10 100 1000];
%alpha = zeros(nclasses,max(mmmres));
%nu   = 0.1;
%sig2(1) = 8;
%if sig2(1)~=8
%  fprintf(1,'attention : pour reproduire le cas 1 classe, il faut sig2=8 ici\n');
%end;
%for nn=1:nchoix  %%% la fin de cette boucle intervient très bas
%  %%% les paramètres de la seconde classe sont fixes : c'est la première qui varie
%  mmm(1)=mmmres(nn); %mmm(1)=50;
%  mmm(1)
%  if nclasses>1
%    mmm(2)=100;
%  end;


%%% boucle sur "nu"
%simu=3;
%faitcontour=0;
%nures=[0.001 0.01 0.1 0.3];
%mmm(1) = 2000;
%sig2(1)=8;
%for nn=1:nchoix  %%% la fin de cette boucle intervient très bas
%  nu=nures(nn);
%  nu


%%%%%%
%%%%%%


  %%% note : pour bien initialiser "minxx" et "maxxx", il faut cette boucle 
  %%%        "cc" supplémentaire
  %%%
  %%% on fait les points à classer
  %%% => ATTENTION : ici, ça ne devrait plus changer par rapport à la
  %%%    version "1 classe" : puisque tous les points ne sont plus pris en compte
  %%%    pour déterminer par où passent les points à classer, comme je
  %%%    faisais dans la première version de ce programme
  %%fprintf(1,"attention : petite diff'erence par rapport au cas `a 1 classe\n");
  npointsaclasser=100;
  datatmp=[];
  for cc=1:nclasses
    datatmp = [datatmp datatout(:,(cc-1)*mmmmax+1:(cc-1)*mmmmax+mmm(cc))];
  end;
  for ii=1:ndimensions
    minxxres(ii)=min(datatmp(ii,:));
    maxxxres(ii)=max(datatmp(ii,:));
  end;
  minxx = minxxres(1) - 2.0;
  maxxx = maxxxres(1) + 2.0;
  for ii=1:npointsaclasser
    point(1,ii) = (ii-1)/(npointsaclasser-1)*(maxxx-minxx) + minxx;
  end;
  for jj=2:ndimensions
    meanyy = mean(datatmp(jj,:));
    for ii=1:npointsaclasser
      point(jj,ii) = meanyy;
    end;
  end;
  clear datatmp;


  %%% boucle sur le nombre de classes
  for cc=1:nclasses
    %%%
    data1=datatout(:,(cc-1)*mmmmax+1:(cc-1)*mmmmax+mmm(cc));


    %%% estimation de "sig2" avec la méthode de Manuel/Jaakkola
    sumnorme=0;
    kk=1;
    for ii=1:mmm(cc)-1
      for jj=ii+1:mmm(cc)
        ddd=data1(:,ii)-data1(:,jj);
        norme=sqrt(sum(ddd.*ddd));
        sumnorme=sumnorme+norme;
        kk=kk+1;
     end;
   end;
   sigestim2=sumnorme/(kk-1)/2;
%   sigestim2=mean(norme)/2;
%   sigestim2=median(norme)/2;
   fprintf(1,"classe:%d   sig2 estim'e=%f\n",cc,2*(sigestim2^2));
   %sig2(1)=26;
   %sig2(2)=26;


    %%% optimisation SMO
    if utilise_code_c==0
      %% code "octave" (lent)
      topt=clock;
      [alphaSMO(cc,:), bbbSMOf(cc)] = lagis_smo_1classe(data1, nu, sig2(cc), topt, '');
      eopt=etime(clock,topt);
    else
      %% code C
      % fichier de données formé
      versionpg = "-v3";

      ff=fopen(["SMO-1-CLASSE/donnees_pour_smo-c" versionpg ".txt"],"wt");
      fprintf(ff,"%d ",ndimensions);
      fprintf(ff,"%d ",mmm(cc));
      fprintf(ff,"%f ",nu);
      fprintf(ff,"%f ",sig2(cc));
      for ii=1:ndimensions
        fprintf(ff,"%f ",data1(ii,:));
      end;
      fclose(ff);
      % calcul
      commande=["SMO-1-CLASSE/lagis_smo_1classe" versionpg];
      topt=clock;  
      system(commande);
      eopt=etime(clock,topt);

      load SMO-1-CLASSE/resultats_de_smo-c-v3.txt;
      alphaSMO(cc,1:mmm(cc)) = resultats_de_smo_c_v3(1:mmm(cc))';
      bbbSMOf(cc)    = resultats_de_smo_c_v3(mmm(cc)+1);
    end;
    bbbSMO(cc) = -bbbSMOf(cc);   %%%%% Inversion à cause des notations différentes
    fprintf(1,'temps requis pour optimiser les SVM sur %d points avec SMO : %f\n', mmm(cc),eopt);


    %%%
    test_asma=1;
    if test_asma==1
      for ii=1:mmm(cc)
        fdecision(ii) = alphaSMO(cc,:)*lagis_rbf_gaussien(data1(:,ii), data1, sig2(cc))' + bbbSMO(cc);
        if fdecision(ii)>0
          decision(ii) =  1;
        else
          decision(ii) = -1;
        end;
      end;

      nvecteursupport = sum(alphaSMO(cc,:)>max(alphaSMO(cc,:))/fractionseuil);
      appartient    = sum(decision==1);
      appartientpas = sum(decision==-1);
      fprintf(1,"\n\nsig2=%f\n\nnombre de vecteur supports : %d\n\n %d echantillons appartiennent a la classe\n %d echantillons n appartiennent pas a la classe\n\n", sig2(cc),nvecteursupport,appartient,appartientpas);
    end;


    %%% détermination de "a1a1" pour la deuxième mesure
    a1a1(cc)=0.0;
    for kk=1:mmm(cc)
      for ll=1:mmm(cc)
        if alphaSMO(cc,kk)>max(alphaSMO(cc,:))/fractionseuil && alphaSMO(cc,ll)>max(alphaSMO(cc,:))/fractionseuil
          a1a1(cc) = a1a1(cc) + alphaSMO(cc,kk)*alphaSMO(cc,ll)*lagis_rbf_gaussien (data1(:,kk), data1(:,ll), sig2(cc));
        end;
      end;
    end;
    a1a1(cc)=sqrt(a1a1(cc));
    c1p1avacos(cc) = -bbbSMO(cc)/a1a1(cc);
    if c1p1avacos(cc)==1
       %%c1p1(cc)=1/(pi/2);     % Il faudrait faire ça
       c1p1(cc) = pi/2;      % certains font plutôt ça (Lionel)
    else
       c1p1(cc)=acos(c1p1avacos(cc));
    end;


    %%% on classe les points à classer
    for ii=1:npointsaclasser
      %%%%%% méthode 1
      %%% score méthode de Lionel numéro 1 du 01/09/06 : 
      %%% fonction de décision f(x) utilisée
      sortiem1(cc+(nn-1)*nchoix,ii) = -50*log(sum(alphaSMO(cc,:).*lagis_rbf_gaussien (point(:,ii), data1, sig2(cc))) + bbbSMO(cc) + 1.0);


      %%%%%% méthode 2
      %%% sorties méthode 2 : "I=acos(c1c2)/acos(c1p1)*normalisation" utilisée
      a1a2=sum(alphaSMO(cc,:).*lagis_rbf_gaussien (point(:,ii), data1, sig2(cc)));
      c1c2tmp = a1a2/a1a1(cc);

      %% on tient compte de "minc1p1" plus tard
      %sortiem2(cc+(nn-1)*nchoix,ii) = acos(c1c2tmp)/acos(c1p1tmp)*2/pi*minc1p1;
      sortiem2(cc+(nn-1)*nchoix,ii) = acos(c1c2tmp)/c1p1(cc)*2/pi;


      %%%%%% méthode 3
      %%% sorties méthode de Lionel numéro 3 du 01/09/06

      %% on tient compte de "bbbMAX" plus tard
      %sortiem3(cc+(nn-1)*nchoix,ii) = 100*( -log(sum(alphaSMO(cc,:).*lagis_rbf_gaussien (point(:,ii), data1, sig2(cc)))) + log(-bbbSMO) - log(-bbbMAX));
      sortiem3(cc+(nn-1)*nchoix,ii) = 100*( -log(sum(alphaSMO(cc,:).*lagis_rbf_gaussien (point(:,ii), data1, sig2(cc)))) + log(-bbbSMO(cc)) );
    end;

  end;   %%% boucle "cc"


  %%% plot du contour
  if faitcontour==1 && ndimensions==2
    %%%
    nzz=100;
    for ii=1:ndimensions
      minxxz(ii)=minxxres(ii)-1.0;
      maxxxz(ii)=maxxxres(ii)+1.0;
    end;

    for ii=1:nzz
      xxt(1,1) = minxxz(1) + (maxxxz(1)-minxxz(1))/(nzz-1)*(ii-1);
      for jj=1:nzz
        xxt(2,1) = minxxz(2) + (maxxxz(2)-minxxz(2))/(nzz-1)*(jj-1);

        xxd1(ii,jj) = xxt(1,1);
        xxd2(ii,jj) = xxt(2,1);

        for dd=1:nclasses
          data1=datatout(:,(dd-1)*mmmmax+1:(dd-1)*mmmmax+mmm(dd));
          zzdSMO(ii+(dd-1)*nzz,jj) = alphaSMO(dd,1:mmm(dd))*lagis_rbf_gaussien(xxt, data1, sig2(dd))' + bbbSMO(dd);
        end;
      end;
    end;


    %%%
    figure(5);
    clg;
    grid on;
    hold on;
    for dd=1:nclasses
      contour(xxd1, xxd2, zzdSMO((dd-1)*nzz+1:dd*nzz,:), [0 0]);
    end;
    %title('Fonctions de décision obtenue avec optimisation SMO');
    hold off;
    switch simu
      case 1
        chaine=[preambule "contour.nc" sprintf("%d",nclasses) ".sig2." sprintf("%f",sig2(1)) ".ps"];
        print(chaine, "-dpsc2");

      case 2
        chaine=[preambule "contour.nc" sprintf("%d",nclasses) ".mmm." sprintf("%f",mmm(1)) ".ps"];
        print(chaine, "-dpsc2");

      case 3
        chaine=[preambule "contour.nc" sprintf("%d",nclasses) ".nu." sprintf("%f",nu) ".ps"];
        print(chaine, "-dpsc2");

      otherwise
        fprintf(1,"Il manque quelque chose - 1\n");
    end;
  end; %%% condition "faitcontour"


  %%% c'est à cet endroit qu'il faut, si nécessaire (si "nclasses>1"),
  %%% calculer "bbbMAX" et "minc1p1", et en tenir compte
  
  if nclasses>1
    % calcul de "minc1p1"
    minc1p1=min(c1p1);
  
    % calcul de "bbbMAX" (voir "lagis_evaluation.m")
    bbbtmp=bbbSMO==0;
    bbbtmp2=bbbSMO-1000*bbbtmp;
    bbbMAX = max(bbbtmp2);

    for cc=1:nclasses
      sortiem2(cc+(nn-1)*nchoix,:) = sortiem2(cc+(nn-1)*nchoix,:)*minc1p1;
      sortiem3(cc+(nn-1)*nchoix,:) = sortiem3(cc+(nn-1)*nchoix,:) - 100*log(-bbbMAX);
    end;
  end;

end;   %%% boucle "nn"


%%% où sont classés les points à classer
if nclasses==2
  for nn=1:nchoix
    for ii=1:npointsaclasser
      % pour la mesure m1
      if sortiem1(1+(nn-1)*nchoix,ii) < sortiem1(2+(nn-1)*nchoix,ii)
        classem1((nn-1)*nchoix+1,ii) = 1;
      else
        classem1((nn-1)*nchoix+1,ii) = 2;
      end;

      % pour la mesure m2
      if sortiem2(1+(nn-1)*nchoix,ii) < sortiem2(2+(nn-1)*nchoix,ii)
        classem2((nn-1)*nchoix+1,ii) = 1;
      else
        classem2((nn-1)*nchoix+1,ii) = 2;
      end;

      % pour la mesure 3
      if sortiem3(1+(nn-1)*nchoix,ii) < sortiem3(2+(nn-1)*nchoix,ii)
        classem3((nn-1)*nchoix+1,ii) = 1;
      else
        classem3((nn-1)*nchoix+1,ii) = 2;
      end;
    end;
  end;

  % notes :
  % 1) ne marche correctement qui si "nchoix=1"
  % 2) "subplot" ne marche pas correctement sous "octave", d'où l'utilisation des petits 
  %    shifts
  % 3) adapté pour le cas où on fait varier "sig2"
  figure(6);
  clf;
  grid on;
  hold on;
  plot(point(1,:), classem1(1,:),    'ro--');
  plot(point(1,:), classem2(1,:)+0.05,'b*--');
  plot(point(1,:), classem3(1,:)+0.10,'g+--');
  ylabel("appartenance du point à la classe 1 ou à la classe 2");
  xlabel("m1, puis m2 déplacé de 0.05 vers le haut, et m3 de 0.1");
  hold off;
  chichi=strcat(preambule, 'classes.nc', sprintf("%d",nclasses), '.sig2.', sprintf("%f",sig2(1)), '.ps');
  print(chichi, "-dpsc2");
end;


%%% autres figures imprimées
chaine1='ob';
chaine2='+g';
figure(1);
clf;
grid on;
hold on;
for cc=1:nclasses
  switch cc
    case 1
      chainefig1=chaine1;

    case 2
      chainefig1=chaine2;

    otherwise
      fprintf(1,"ajouter des styles de ligne, puisqu il y a plus de classes que prevu\n");
      chainefig1=chaine1;
  end;

  data1=datatout(:,(cc-1)*mmmmax+1:(cc-1)*mmmmax+mmm(cc));
  plot(data1(1,:), data1(2,:), chainefig1);
  hold on;
end;
plot(point(1,:), point(2,:), '*r');
hold off;
chichi=strcat(preambule, 'donnees.nc', sprintf("%d",nclasses), '.ps');
print(chichi, "-dpsc2");


%%% pour les styles de lignes
%%% 1. : classe 1
%%% 2. : classe 2
%%% => s'il y en a plus de 2, ajouter des "chaine"s
chainea(1,:) = 'r+;  ';
chaineb(1,:) = 'g*;  ';
chainec(1,:) = 'ko;  ';
chained(1,:) = 'bx;  ';
chainea(2,:) = 'r+--;';
chaineb(2,:) = 'g*--;';
chainec(2,:) = 'ko--;';
chained(2,:) = 'bx--;';


%%% plots
switch simu
  %%% "sig2"
  case 1
    figure(2);
    clf;
    grid on;
    hold on;
    xlabel("mesure 1");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem1(cc,:),chainea(cc,1:2));%strcat(chainea(cc,:)));%, num2str(sig2res(1)), ';'));
          case 2
            plot(point(1,:),sortiem1(cc+nchoix,:),chaineb(cc,1:2));%strcat(chaineb(cc,:)));%, num2str(sig2res(2)), ';'));
          case 3
            plot(point(1,:),sortiem1(cc+2*nchoix,:),chainec(cc,1:2));%strcat(chainec(cc,:)));%, num2str(sig2res(3)), ';'));
          case 4
            plot(point(1,:),sortiem1(cc+3*nchoix,:),chained(cc,1:2));%strcat(chained(cc,:)));%, num2str(sig2res(4)), ';'));
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm1.nc', sprintf("%d",nclasses), '.sig2.', sprintf("%f",sig2(1)), '.ps');
    print(chichi, "-dpsc2");

    figure(3);
    clf;
    grid on;
    hold on;
    xlabel("mesure 2");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem2(cc,:),chained(cc,1:2));%[chained(cc,:) num2str(sig2res(1)) ';']);
          case 2
            plot(point(1,:),sortiem2(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(sig2res(2)) ';']);
          case 3
            plot(point(1,:),sortiem2(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(sig2res(3)) ';']);
          case 4
            plot(point(1,:),sortiem2(cc+3*nchoix,:),chained(cc,1:2));%,[chained(cc,:) num2str(sig2res(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm2.nc', sprintf("%d",nclasses), '.sig2.', sprintf("%f",sig2(1)), '.ps');
    print(chichi,"-dpsc2");

    figure(4);
    clf;
    grid on;
    hold on;
    xlabel("mesure 3");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem3(cc,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(sig2res(1)) ';']);
          case 2
            plot(point(1,:),sortiem3(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(sig2res(2)) ';']);
          case 3
            plot(point(1,:),sortiem3(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(sig2res(3)) ';']);
          case 4
            plot(point(1,:),sortiem3(cc+3*nchoix,:),chained(cc,1:2));%[chained(cc,:) num2str(sig2res(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm3.nc', sprintf("%d",nclasses), '.sig2.', sprintf("%f",sig2(1)), '.ps');
    print(chichi,"-dpsc2");

  %%% "mmm"
  case 2
    figure(2);
    clf;
    grid on;
    hold on;
    xlabel("mesure 1");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem1(cc,:),chainea(cc,1:2));%,[chainea(cc,:) num2str(mmmres(1)) ';']);
          case 2
            plot(point(1,:),sortiem1(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(mmmres(2)) ';']);
          case 3
            plot(point(1,:),sortiem1(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(mmmres(3)) ';']);
          case 4
            plot(point(1,:),sortiem1(cc+3*nchoix,:),chained(cc,1:2));%,[chained(cc,:) num2str(mmmres(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm1.nc', sprintf("%d",nclasses), '.mmm.', sprintf("%f",mmm), '.ps');
    print(chichi, "-dpsc2");

    figure(3);
    clf;
    grid on;
    hold on;
    xlabel("mesure 2");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem2(cc,:),chainea(cc,1:2));%,[chainea(cc,:) num2str(mmmres(1)) ';']);
          case 2
            plot(point(1,:),sortiem2(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(mmmres(2)) ';']);
          case 3
            plot(point(1,:),sortiem2(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(mmmres(3)) ';']);
          case 4
            plot(point(1,:),sortiem2(cc+3*nchoix,:),chained(cc,1:2));%,[chained(cc,:) num2str(mmmres(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm2.nc', sprintf("%d",nclasses), '.mmm.', sprintf("%f",mmm), '.ps');
    print(chichi,"-dpsc2");

    figure(4);
    clf;
    grid on;
    hold on;
    xlabel("mesure 3");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem3(cc,:),chainea(cc,1:2));%,[chainea(cc,:) num2str(mmmres(1)) ';']);
          case 2
            plot(point(1,:),sortiem3(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(mmmres(2)) ';']);
          case 3
            plot(point(1,:),sortiem3(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(mmmres(3)) ';']);
          case 4
            plot(point(1,:),sortiem3(cc+3*nchoix,:),chained(cc,1:2));%,[chained(cc,:) num2str(mmmres(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm3.nc', sprintf("%d",nclasses), '.mmm.', sprintf("%f",mmm), '.ps');
    print(chichi,"-dpsc2");

  %%% "nu"
  case 3
    figure(2);
    clf;
    grid on;
    hold on;
    xlabel("mesure 1");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem1(cc,:),chainea(cc,1:2));%,[chainea(cc,:) num2str(nures(1)) ';']);
          case 2
            plot(point(1,:),sortiem1(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(nures(2)) ';']);
          case 3
            plot(point(1,:),sortiem1(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(nures(3)) ';']);
          case 4
            plot(point(1,:),sortiem1(cc+3*nchoix,:),chained(cc,1:2));%,[chained(cc,:) num2str(nures(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm1.nc', sprintf("%d",nclasses), '.nu.', sprintf("%f",nu), '.ps');
    print(chichi, "-dpsc2");

    figure(3);
    clf;
    grid on;
    hold on;
    xlabel("mesure 2");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem2(cc,:),chainea(cc,1:2));%,[chainea(cc,:) num2str(nures(1)) ';']);
          case 2
            plot(point(1,:),sortiem2(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(nures(2)) ';']);
          case 3
            plot(point(1,:),sortiem2(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(nures(3)) ';']);
          case 4
            plot(point(1,:),sortiem2(cc+3*nchoix,:),chained(cc,1:2));%,[chained(cc,:) num2str(nures(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm2.nc', sprintf("%d",nclasses), '.nu.', sprintf("%f",nu), '.ps');
    print(chichi,"-dpsc2");

    figure(4);
    clf;
    grid on;
    hold on;
    xlabel("mesure 3");
    for cc=1:nclasses
      for dd=1:nchoix
        switch dd
          case 1
            plot(point(1,:),sortiem3(cc,:),chainea(cc,1:2));%,[chainea(cc,:) num2str(nures(1)) ';']);
          case 2
            plot(point(1,:),sortiem3(cc+nchoix,:),chaineb(cc,1:2));%,[chaineb(cc,:) num2str(nures(2)) ';']);
          case 3
            plot(point(1,:),sortiem3(cc+2*nchoix,:),chainec(cc,1:2));%,[chainec(cc,:) num2str(nures(3)) ';']);
          case 4
            plot(point(1,:),sortiem3(cc+3*nchoix,:),chained(cc,1:2));%,[chained(cc,:) num2str(nures(4)) ';']);
        end;
      end;
    end;
    hold off;
    chichi=strcat(preambule, 'm3.nc', sprintf("%d",nclasses), '.nu.', sprintf("%f",nu), '.ps');
    print(chichi,"-dpsc2");

  otherwise
    fprintf(1,"Il manque quelque chose - 2\n");
end; %%%



%%% fin