%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la méthode SMO - pour SVM 1 classe
% Code Lagis et CentraleSupélec
%
% - "examineExemple_1classe" est appelée par le programme principal
%   "lagis_smo_1classe.m"
%
% liste TODO
% ----------
%
% - 02/05/2006
%   Description des arguments et des retours
%
% - 02/05/2006
%   Accélérer !!!!!!!!!!
%
% Stéphane Rossignol - 27/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

function [resultat, alpha_nouv, bbb_nouv, erreur_cache_nouv] = examineExemple_1classe(i1, alpha, bbb, xxx, CCC, erreur_cache, tolerance, sig2, mmm, path1)


%%% for the SVM computation, we need to add these "paths"
strversion=version;
if strcmp(version,'2.1.73')==1
  path(LOADPATH,path1);
else
  addpath(path1);
end;


%%%
resultat_tS = 0;
resultat    = 0;
alpha_nouv  = alpha;
bbb_nouv    = bbb;
erreur_cache_nouv = erreur_cache;


%%%
alpha1 = alpha(i1);

if (alpha1>0.0 && alpha1<CCC)
   E1 = erreur_cache(i1);
else
   E1 = sage_fonc_1classe(i1, alpha, bbb, xxx, sig2, mmm, path1);
end;
%fprintf(1,'\ni1=%d alpha1=%f E1=%f(%f)', i1,alpha1,E1,sage_fonc_1classe(i1, alpha, bbb, xxx, sig2, mmm, path1));

if ( (E1<-tolerance && alpha1<CCC) || (E1>tolerance && alpha1>0.0))
   % méthode 1
   tmax = 0;
   i2 = -1;
   for kk=1:mmm
      if (alpha(kk)>0 && alpha(kk)<CCC)
         E2 = erreur_cache(kk);
         temp = abs(E1 - E2);
%         fprintf(1,'\n kk=%d E2=%f temp=%f tmax=%f',kk,E2,temp,tmax);
         if (temp>tmax)
            tmax = temp;
            i2 = kk;
         end;
      end;
   end;
   if (i2>=0)
      [resultat_tS, alpha_nouv, bbb_nouv, erreur_cache_nouv] = pas_1classe(i1, i2, alpha, xxx, CCC, erreur_cache, bbb, sig2, mmm, path1);
      alpha = alpha_nouv;
      bbb   = bbb_nouv;
      erreur_cache = erreur_cache_nouv;
      if (resultat_tS==1)
         resultat=1;       %%% on sort tout de suite de "examineExample"
      end;
   end;
   %%fprintf(1,'\ni2=%d resultat_tS(1)=%d\n',i2,resultat_tS);

   %%fprintf(1,'\nCCC=%f\n',CCC);
   %%for ii=1:mmm
   %%  fprintf(1,'alpha=%f\n',alpha(ii));
   %%end;
   % méthode 2
   if (resultat==0)
      k0 = floor(mmm*rand(1))+1;
      for kk=k0:k0+mmm-1
         if kk>mmm
            i2=kk-mmm;
         else
            i2=kk;
         end;
         if (alpha(i2)>0.0 && alpha(i2)<CCC)
%            fprintf(1,'methode 2 active - i2=%d\n',i2);
            [resultat_tS, alpha_nouv, bbb_nouv, erreur_cache_nouv] = pas_1classe(i1, i2, alpha, xxx, CCC, erreur_cache, bbb, sig2, mmm, path1);
            alpha = alpha_nouv;
            bbb   = bbb_nouv;
            erreur_cache = erreur_cache_nouv;
            if (resultat_tS==1)
               resultat=1;       %%% on sort tout de suite de "examineExemple"
            end;
         end;
      end;
   end;
%   for ii=1:mmm
%     fprintf(1,'alpha_nouv=%f\n',alpha_nouv(ii));
%   end;
%   fprintf(1,'resultat_tS(2)=%d\n',resultat_tS);

   % méthode 3
   if (resultat==0)
      k0 = floor(mmm*rand(1))+1;
      for kk=k0:k0+mmm-1
         if kk>mmm
            i2=kk-mmm;
         else
            i2=kk;
         end;
%         fprintf(1,'methode 3 active - i2=%d\n',i2);
         [resultat_tS, alpha_nouv, bbb_nouv, erreur_cache_nouv] = pas_1classe(i1, i2, alpha, xxx, CCC, erreur_cache, bbb, sig2, mmm, path1);
         alpha = alpha_nouv;
         bbb   = bbb_nouv;
         erreur_cache = erreur_cache_nouv;
         if (resultat_tS==1)
            resultat=1;       %%% on sort tout de suite de "examineExample"
         end;
      end;
   end;
%   for ii=1:mmm
%     fprintf(1,'alpha_nouv=%f\n',alpha_nouv(ii));
%   end;
%   fprintf(1,'resultat_tS(3)=%d\n',resultat_tS);
end;

