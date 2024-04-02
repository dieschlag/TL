%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la méthode SMO - pour SVM 2 classes
% Code Lagis et CentraleSupélec
%
% - examineExample
%
% liste TODO
% ----------
%
% Stéphane Rossignol - 25/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

function [resultat, alpha_new, bbb_new, error_cache_new] = examineExample(i1, alpha, bbb, xxx, yyy, CCC, error_cache, tolerance, sig2)


%%%
resultat  = 0;
alpha_new = alpha;
bbb_new   = bbb;
error_cache_new = error_cache;


%%%
mmm = length(alpha);


%%%
y1 = yyy(i1);
alpha1 = alpha(i1);

if (alpha1>0.0 && alpha1<CCC)
   E1 = error_cache(i1);
else
   E1 = learned_func(i1, alpha, bbb, xxx, yyy, sig2) - y1;
end;

r1 = y1*E1;

if ( (r1<-tolerance && alpha1<CCC) || (r1>tolerance && alpha1>0.0))
   % méthode 1
   tmax = 0;
   i2 = -1;
   for kk=1:mmm
      if (alpha(kk)>0 && alpha(kk)<CCC)
         E2 = error_cache(kk);
         temp = abs(E1 - E2);
         if (temp>tmax)
            tmax = temp;
            i2 = kk;
         end;
      end;
   end;
   if (i2>=0)
      [resultat_tS, alpha_new, bbb_new, error_cache_new] = takeStep(i1, i2, alpha, xxx, yyy, CCC, error_cache, bbb, sig2);
      alpha = alpha_new;
      bbb   = bbb_new;
      error_cache = error_cache_new;
      if (resultat_tS==1)
         resultat=1;       %%% on sort tout de suite de "examineExample"
      end;
   end;

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
            [resultat_tS, alpha_new, bbb_new, error_cache_new] = takeStep(i1, i2, alpha, xxx, yyy, CCC, error_cache, bbb, sig2);
            alpha = alpha_new;
            bbb   = bbb_new;
            error_cache = error_cache_new;
            if (resultat_tS==1)
               resultat=1;       %%% on sort tout de suite de "examineExample"
            end;
         end;
      end;
   end;
   
   % méthode 3
   if (resultat==0)
      k0 = floor(mmm*rand(1))+1;
      for kk=k0:k0+mmm-1
         if kk>mmm
            i2=kk-mmm;
         else
            i2=kk;
         end;
         [resultat_tS, alpha_new, bbb_new, error_cache_new] = takeStep(i1, i2, alpha, xxx, yyy, CCC, error_cache, bbb, sig2);
         alpha = alpha_new;
         bbb   = bbb_new;
         error_cache = error_cache_new;
         if (resultat_tS==1)
            resultat=1;       %%% on sort tout de suite de "examineExample"
         end;
      end;
   end;
end;

