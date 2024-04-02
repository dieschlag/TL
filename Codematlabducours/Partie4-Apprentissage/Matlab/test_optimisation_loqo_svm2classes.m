%%%%%%%%%%%%%%%%%%%%
% Test l'optimisation LOQO seule, cas "SVM 2 classes"
% SVM 2 classes - avec noyau autre que le produit scalaire
% Code Lagis et CentraleSupélec
%
%
% liste TODO
% ----------
%
% Stéphane Rossignol - (voir "test_svm_loqo_2classes_gaussien.m") 24/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

clear all;
close all;


%%% ne pas changer ces paramètres
use_octave  = 1; % utilisé pour LOQO
calcul_loqo = 1; % calcul avec l'optimiseur LOQO


%%% dans le but de tester, je m'arrange pour toujours avoir les mêmes données
%%% et pour toujours passer par le même chemin d'optimisation
%%% => commenter ces lignes si on ne veut pas de ça
if (use_octave==0)
   rand('state',1);
else
   rand("seed",1);
   randn("seed",1);
end;


%%% formation des données
nn = 100;     % nombre de points par classe ; modifié possiblement plus loin
numexemple=3; % deux exemples sont possibles (1 et 2) : tester les deux
switch numexemple
  case 1
    %%% classe 1
    cl1_d1 = rand(1,nn);      % première dimension de la première classe
    cl1_d2 = rand(1,nn);      % deuxième dimension de la première classe

    a=-1.0;
    b= 2.0;               % ~contrôle la distance entre les deux classes

    for ii=1:nn
      while (cl1_d2(ii)>a*cl1_d1(ii)+b)
        cl1_d2(ii)=cl1_d2(ii)/2.0;
      end;
    end;

    %%% classe 2
    cl2_d1 = rand(1,nn);      % première dimension de la deuxième classe
    cl2_d2 = rand(1,nn);      % deuxième dimension de la deuxième classe

    dist = 0.1;

    for ii=1:nn
      while (cl2_d2(ii)<a*cl2_d1(ii)+b+dist)
        cl2_d2(ii)=cl2_d2(ii)*2.0;
      end;
    end;

    sig2 = 0.01; %%% paramètre du noyau RBF gaussien : jouez avec ce paramètre, pour voir son influence

  case 2
    %%% classe 1
    cl1_d1 = 0.0 + 0.1*randn(1,nn);
    cl1_d2 = 0.0 + 0.2*randn(1,nn);

    %%% classe 2
    cl2_d1 = 0.0 + 0.1*randn(1,nn);
    cl2_d2 = 0.8 + 0.2*randn(1,nn);

    sig2 = 0.1; %%% paramètre du noyau RBF gaussien : jouez avec ce paramètre, pour voir son influence

  %%% pour le problème des oiseaux
  case 3
    don1=load('featuresoiseau1.txt');
    don2=load('featuresoiseau2.txt');
    don3=load('featuresoiseau3.txt');

    %nn=length(don1);
    nn=1000;

    %%% classe 1 : oiseau 1
    cl1_d1 = don1(1:nn,1)';
    cl1_d2 = don1(1:nn,2)';

    %%% classe 2 : oiseau 2
    cl2_d1 = don2(1:nn,1)';
    cl2_d2 = don2(1:nn,2)';

    %%% classe 3 : oiseau 3 (pas utilisé ici, à part pour normaliser)
    cl3_d1 = don3(1:nn,1)';
    cl3_d2 = don3(1:nn,2)';

    moy1=mean([cl1_d1 cl2_d1 cl3_d1]);
    std1=std([cl1_d1 cl2_d1 cl3_d1]);

    moy2=mean([cl1_d2 cl2_d2 cl3_d2]);
    std2=std([cl1_d2 cl2_d2 cl3_d2]);

    %%% classe 1
    cl1_d1 = (don1(1:nn,1)'-moy1)/std1;
    cl1_d2 = (don1(1:nn,2)'-moy2)/std2;

    %%% classe 2
    cl2_d1 = (don3(1:nn,1)'-moy1)/std1;
    cl2_d2 = (don3(1:nn,2)'-moy2)/std2;

    sig2 = 2.;
end;

mind1 = min([cl1_d1 cl2_d1]);
maxd1 = max([cl1_d1 cl2_d1]);
mind2 = min([cl1_d2 cl2_d2]);
maxd2 = max([cl1_d2 cl2_d2]);


%%% calcul de la solution : optimisation sous contraintes, méthode LOQO
donnees = [[cl1_d1 cl2_d1]' [cl1_d2 cl2_d2]']'; %%% données
classes = [ones(1,nn) -ones(1,nn)];             %%% labels des 2 classes

if (calcul_loqo==1)
  ttloqo=clock;
  for ii=1:2*nn
    if (ii<=nn)
      vaii1=cl1_d1(ii);
      vaii2=cl1_d2(ii);
    else
      vaii1=cl2_d1(ii-nn);
      vaii2=cl2_d2(ii-nn);
    end;
    for jj=1:2*nn
      if (jj<=nn)
        vajj1=cl1_d1(jj);
        vajj2=cl1_d2(jj);
      else
        vajj1=cl2_d1(jj-nn);
        vajj2=cl2_d2(jj-nn);
      end;

      % noyau RBF gaussien
      matsca(ii,jj) = lagis_rbf_gaussien([vaii1 vaii2]', [vajj1 vajj2]', sig2);
    end;
  end;
  ppx=-ones(1,nn);
  ppy= ones(1,nn);
  matpp=[ppx ppy]'*[ppx ppy];

  c = -ones(2*nn,1);
  H = matpp.*matsca;
  A = [ppx ppy];
  b = 0;
  l = zeros(2*nn,1);
  u = 16.0*ones(2*nn,1);
  [alphaloqo,yloqo] = pr_loqo3(c, H, A, b, l, u, use_octave);
  eeloqo=etime(clock,ttloqo);
else
  eeloqo=-1;
end;

fprintf(1,'taille des donnees : %d\n',2*nn);
fprintf(1,'temps requis pour LOQO : %f\n',eeloqo);

if (calcul_loqo==1)
  alphaxloqo = alphaloqo(1:nn)';        % alphas pour la première classe
  alphayloqo = alphaloqo(nn+1:2*nn)';   % alphas pour la deuxième classe

  figure(2);
  clf;
  subplot(211);
  title('alphas obtenus avec LOQO');
  grid on;
  hold on;
  plot(alphaxloqo,'r');
  hold off;
  subplot(212);
  grid on;
  hold on;
  plot(alphayloqo,'b');
  hold off;
end;


%%%%%%%%%%%%%%%%%%%%%
alpha=alphaloqo;
alphax=alphaxloqo;
alphay=alphayloqo;

%%% calcul de la solution : calcul de w
ww1r = sum(-alphax.*cl1_d1 + alphay.*cl2_d1);
ww2r = sum(-alphax.*cl1_d2 + alphay.*cl2_d2);
normal = mean([ww1r ww2r]);
ww1=ww1r/normal;      % première dimension de w
ww2=ww2r/normal;      % deuxième dimension de w

bbbopt=1/(sqrt(ww1r^2 + ww2r^2));

tmpv2 = alpha.*classes';
for ii=1:nn
  cl1_v(ii) = sign( tmpv2'*lagis_rbf_gaussien([cl1_d1(ii) cl1_d2(ii)]', donnees, sig2)' + bbbopt );
  cl2_v(ii) = sign( tmpv2'*lagis_rbf_gaussien([cl2_d1(ii) cl2_d2(ii)]', donnees, sig2)' + bbbopt );
end;
bbb=bbbopt;

figure(3);
clf;
subplot(211);
title('sortie du SVM pour les points de la classe 1');
grid on;
hold on;
plot(cl1_v,'r');
hold off;
subplot(212);
title('sortie du SVM pour les points de la classe 2');
grid on;
hold on;
plot(cl2_v,'b');
hold off;


%%% au cas où
if bbb==-10000000
  fprintf(1,'b n a pas pu etre calcule correctement\n');
  bbb=-1.0
end;


%%% plot pour la route - LOQO
if (calcul_loqo==1)
  vsloqocl1 = alphaxloqo>max(alphaloqo)/20.0;
  vsloqocl2 = alphayloqo>max(alphaloqo)/20.0;
  figure(4);
  clf;
  grid on;
  hold on;
  title('les points magentas sont les vecteurs supports trouves par LOQO')
  plot(cl1_d1,cl1_d2,'*r');
  plot(cl2_d1,cl2_d2,'ob');
  for ii=1:nn
    if (vsloqocl1(ii)>0)
      plot(cl1_d1(ii),cl1_d2(ii),'+m');
    end;
    if (vsloqocl2(ii)>0)
      plot(cl2_d1(ii),cl2_d2(ii),'+m');
    end;
  end;
  xlabel('dimension 1');
  ylabel('dimension 2');
  hold off;
end;


%%% calcul et plot de la fonction de décision obtenue
nzz=500;

mind1c=mind1;
maxd1c=maxd1;
pasd1=(maxd1c-mind1c)/nzz;

mind2c=mind2;
maxd2c=maxd2;
pasd2=(maxd2c-mind2c)/nzz;

ppd1res = [];
ppd2res = [];
ppd1=mind1c;
for ii=1:nzz
  ppd2=mind2c;
  for jj=1:nzz
     clnew(ii,jj) = sign( tmpv2'*lagis_rbf_gaussien([ppd1 ppd2]', donnees, sig2)' + bbbopt );
     if ii==1
       ppd2res=[ppd2res ppd2];
     end;
    ppd2=ppd2+pasd2;
  end;
  ppd1res=[ppd1res ppd1];
  ppd1=ppd1+pasd1;
end;

%%% plot pour la route
figure(6);
clf;
title('frontiere entre les 2 classes');
hold on;
imagesc(ppd1res,ppd2res,-clnew');
% attention : "mesh" inverse l'axe des abscisses et l'axe des ordonnées
% => même comportement qu'avec "matlab", donc :-(
%mesh(ppd2res,ppd1res,clnew);
grid on;
hold on;
plot(cl1_d1,cl1_d2,'*r');
plot(cl2_d1,cl2_d2,'ob');
xlabel('dimension 1');
ylabel('dimension 2');
xlim([min(ppd1res) max(ppd1res)]);
ylim([min(ppd2res) max(ppd2res)]);
hold off;
colormap('cool');

