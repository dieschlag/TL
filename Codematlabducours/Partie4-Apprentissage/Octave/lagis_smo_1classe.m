%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la méthode SMO - pour SVM 1 classe
% Code Lagis et CentraleSupélec
%
% - Programme principal
%   Est accompagné de :
%   * examineExemple_1classe.m
%   * pas_1classe.m
%   * sage_fonc_1classe.m
%   * noyau_fonc_1classe.m
%
% Note : une version en C existe : voir "lagis_smo_1classe.c"
%
% Ce code est tiré du document :
% "Sequential Minimal Optimization for SVM" de XXX
% où le programme complet est donné en C++
%
% Arguments :
% -----------
% - "xxx"  : les données, de dimensions:
%     'N'   lignes   (nombre de dimensions)
%     'mmm' colonnes (nombre de données d'apprentissage)
% - "CCC"  : la borne supérieure de chaque "alpha"
% - "sig2" : le paramètre du noyau RBF gaussien, noyau
%     qui est le seul possible pour l'instant
% - "tdeb" : instant de commencement de l'optimisation
%
% Retours :
% ---------
%
% liste TODO
% ----------
%
% 1) 27/04/2006
%    A Finir !!!!!!!!! Voir les "????"
%
% 2) 27/04/2006
%    Décrire en détails les arguments et les retours
%
% 3) 27/04/2006
%    Trouver un moyen élégant d'étendre ceci à plus de
%    noyaux
%
% 4) 27/04/2006
%    Si c'est possible (et vraiment nécessaire), l'accélérer !!!!!!
%
% 5) 27/04/2006
%    Comprendre ce qui se passe en ce qui concerne le temps mis pour 
%    optimiser en fonction de la taille des données
%
% Stéphane Rossignol - 27/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%


function [alpha_nouv, bbb_nouv] = lagis_smo_1classe(xxx, nu, sig2, tdeb, path1);


%%% for the SVM computation, we need to add these "paths"
strversion=version;
if strcmp(version,'2.1.73')==1
  path(LOADPATH,path1);
else
  addpath(path1);
end;


%%% Initialisations
mmm = size(xxx,2);
nombreChange = 0;
examineTout = 1;
tolerance = 0.001;

%
nummm = nu*mmm;
CCC=1/nummm;
nnn = floor(nummm);
alpha = zeros(1,mmm);
%tab = randperm(mmm);
tab = [1:mmm];
alpha(tab(1:nnn)) = 1/nummm;
alpha(tab(nnn+1)) = 1 - nnn/nummm;
alpha_nouv = alpha;

% initialisation d'"erreur_cache" et de "bbb" : est-ce que je fais ça bien ?
bbb = 0;
for ii=1:mmm
   erreur_cache(ii)=0;
   if (alpha(ii)>0.0 && alpha(ii)<CCC)
      for jj=1:mmm
         erreur_cache(ii) = erreur_cache(ii) + alpha(jj)*noyau_fonc_1classe(jj,ii,xxx,sig2,path1);
      end;
      if (alpha(ii)>0 && erreur_cache(ii)>bbb)
         bbb=erreur_cache(ii);
      end;
   end;
end;
%fprintf(1,'b est initialise a %f\n',bbb);


%%%
nexaminetout = 0;
nexaminepastout = 0;

%%% Boucle principale
while (nombreChange>0 || examineTout==1)
   nombreChange=0;
   if (examineTout==1)
      nexaminetout++;
      ee=etime(clock,tdeb);
      %
      fprintf(1,'examine tout - SMO (%d - %f) (be patient !!!\n', nexaminetout,ee);
      fprintf(1,' and have a look to "data_files/examineTout.txt")\n');flush(1);
      %
      flog = fopen('LOGFILES/lagis_test_prototype.log','a');
      fprintf(flog,'examine tout - SMO (%d - %f) (be patient !!!\n', nexaminetout,ee);
      fclose(flog);

      ff=fopen('data_files/examineTout.txt','wt');
      for kk=1:mmm
%         fprintf(ff,'%d sur %d\n',kk,mmm);
         %%fprintf(1,'\nalpha avant : ');            
         %%fprintf(1,'%f ',alpha);
         %path(DEFAULT_LOADPATH,'/home/rossigno/LAGIS-CODE/SMO-1-CLASSE/');
         [change, alpha_nouv, bbb_nouv, erreur_cache_nouv] = examineExemple_1classe(kk, alpha, bbb, xxx, CCC, erreur_cache, tolerance, sig2, mmm, path1);
         %%fprintf(1,'\nalpha apres : ');
         %%fprintf(1,'%f ',alpha);
         %%fprintf(1,'\nalpha_nouv : ');
         %%fprintf(1,'%f ',alpha_nouv);
         %fprintf(1,'\nerreur_cache : ');
         %fprintf(1,'%f ',erreur_cache);
         %fprintf(1,'\nerreur_cache_nouv : ');
         %fprintf(1,'%f ',erreur_cache_nouv);
         %fprintf(1,'\n');
         %pause;
         alpha = alpha_nouv;
         bbb = bbb_nouv;
         erreur_cache = erreur_cache_nouv;
         nombreChange = nombreChange + change;
      end;
      fclose(ff);
   else
      nexaminepastout++;
      ee=etime(clock,tdeb);
      %
      fprintf(1,'examine pas tout - SMO (%d - %f)\n',nexaminepastout,ee);flush(1);
      %
      flog = fopen('LOGFILES/lagis_test_prototype.log','a');
      fprintf(flog,'examine pas tout - SMO (%d - %f)\n',nexaminepastout,ee);
      fclose(flog);

      for kk=1:mmm
         if (alpha(kk)!=0 && alpha(kk)!=CCC)
            [change, alpha_nouv, bbb_nouv, erreur_cache_nouv] = examineExemple_1classe(kk, alpha, bbb, xxx, CCC, erreur_cache, tolerance, sig2, mmm, path1);
            alpha = alpha_nouv;
            bbb = bbb_nouv;
            erreur_cache = erreur_cache_nouv;
            nombreChange = nombreChange + change;
         end;
      end;
   end;

   if (examineTout==1)
      examineTout = 0;
   elseif (nombreChange==0)
      examineTout = 1;
   end;
end;

