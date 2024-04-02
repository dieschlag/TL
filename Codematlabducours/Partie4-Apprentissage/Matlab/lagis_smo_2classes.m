%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la méthode SMO - pour SVM 2 classes
% Code Lagis et CentraleSupélec
%
% - programme principal
%
% Ce code est tiré du document :
% "Sequential Minimal Optimization for SVM" de ...
% où le programme complet est donné en C++
%
% Arguments :
% -----------
% - "xxx"  : les données, de dimensions:
%     'N'   lignes   (nombre de dimensions)
%     'mmm' colonnes (nombre de données d'apprentissage)
% - "yyy"  : les labels
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
% 1) 26/04/2006
%    Aussi vite que possible, modifier ceci pour le cas "1 classe" !!!
%
% 2) 26/04/2006
%    Décrire en détails les arguments et les retours
%
% 3) 26/04/2006
%    Trouver un moyen élégant d'étendre ceci à plus de noyaux
%
% 4) 26/04/2006
%    Si c'est possible (et vraiment nécessaire), l'accélérer
%
% 5) 26/04/2006
%    Comprendre ce qui se passe en ce qui concerne le temps mis pour 
%    optimiser en fonction de la taille des données
%
% Stéphane Rossignol - 24/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

function [alpha_new, bbb_new] = lagis_smo_2classes(xxx, yyy, CCC, sig2, tdeb);


%%% Initialisations
mmm = size(xxx,2);
numChanged = 0;
examineAll = 1;
error_cache = -yyy;
bbb = 0;
tolerance = 0.001;
alpha = zeros(mmm,1);
alpha_new = alpha;
nexaminall = 0;
nexaminnotall = 0;


%%% Boucle principale
while (numChanged>0 || examineAll==1)
   numChanged=0;
   if (examineAll==1)
      nexaminall=nexaminall+1;
      ee=etime(clock,tdeb);
      fprintf(1,"examine all - SMO (%d - %f)\n",nexaminall,ee);flush(1);
      for kk=1:mmm
         [change, alpha_new, bbb_new, error_cache_new] = examineExample(kk, alpha, bbb, xxx, yyy, CCC, error_cache, tolerance, sig2);
         %fprintf(1,"\nalpha : ");          fprintf(1,"%f ",alpha);
         %fprintf(1,"\nalpha_new : ");      fprintf(1,"%f ",alpha_new);
         %fprintf(1,"\nerror_cache : ");    fprintf(1,"%f ",error_cache);
         %fprintf(1,"\nerror_cache_new : ");fprintf(1,"%f ",error_cache_new);
         %fprintf(1,"\n");
         %pause;
         alpha = alpha_new;
         bbb = bbb_new;
         error_cache = error_cache_new;
         numChanged = numChanged + change;
      end;
   else
      nexaminnotall=nexaminnotall+1;
      ee=etime(clock,tdeb);
      fprintf(1,"examine not all - SMO (%d - %f)\n",nexaminnotall,ee);flush(1);
      for kk=1:mmm
         if (alpha(kk)~=0 && alpha(kk)~=CCC)
            [change, alpha_new, bbb_new, error_cache_new] = examineExample(kk, alpha, bbb, xxx, yyy, CCC, error_cache, tolerance, sig2);
            alpha = alpha_new;
            bbb = bbb_new;
            error_cache = error_cache_new;
            numChanged = numChanged + change;
         end;
      end;
   end;

   if (examineAll==1)
      examineAll = 0;
   elseif (numChanged==0)
      examineAll = 1;
   end;
end;

