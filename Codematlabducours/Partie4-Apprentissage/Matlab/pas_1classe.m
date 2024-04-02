%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la méthode SMO - pour SVM 1 classe
% Code Lagis et CentraleSupélec
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

function [resultat, alpha_nouv, bbb_nouv, erreur_cache_nouv] = pas_1classe (i1, i2, alpha, xxx, CCC, erreur_cache, bbb, sig2, mmm, path1)


%%% for the SVM computation, we need to add these "paths"
strversion=version;
if strcmp(version,'2.1.73')==1
  path(LOADPATH,path1);
else
  addpath(path1);
end;


%%%
alpha_nouv = alpha;
bbb_nouv = bbb;
erreur_cache_nouv = erreur_cache;
resultat = 1;


%%%
epsi = 1e-3;


%%%
if (i1==i2)
   resultat=0;
end;


if (resultat==1)
   %%%
   alpha1 = alpha(i1);
   if (alpha1>0.0 && alpha1<CCC)
      E1 = erreur_cache(i1);
   else
      E1 = sage_fonc_1classe(i1, alpha, bbb, xxx, sig2, mmm, path1);
   end;

   alpha2 = alpha(i2);
   if (alpha2>0.0 && alpha2<CCC)
      E2 = erreur_cache(i2);
   else
      E2 = sage_fonc_1classe(i2, alpha, bbb, xxx, sig2, mmm, path1);
   end;


   %%%
   gamma = alpha1+alpha2;
   if (gamma>CCC)
      LLL = gamma-CCC;
      HHH = CCC;
   else
      LLL = 0.0;
      HHH =gamma;
   end;

   if (LLL==HHH)
      resultat=0;
   end;
end;



if (resultat == 1)
   %%%
   %k11 = noyau_fonc_1classe(i1, i1, xxx, sig2, path1);
   k12 = noyau_fonc_1classe(i1, i2, xxx, sig2, path1);
   %k22 = noyau_fonc_1classe(i2, i2, xxx, sig2, path1);
   %etamine = 2.0*k12 - k11 - k22;

   k11=1.0;k22=1.0;
   etamine = 2.0*k12 - 2.0;


   %%%
   if (etamine<0.0)
      a2 = alpha2 + (E2-E1)/etamine;
      if (a2<LLL)
         a2 = LLL;
      elseif a2>HHH
         a2 = HHH;
      end;
   else
      c1 = etamine/2.0;
      c2 = (E1-E2) - etamine*alpha2;
      Lobj = c1*LLL*LLL + c2*LLL;
      Hobj = c1*HHH*HHH + c2*HHH;

      if (Lobj>Hobj+epsi)
         a2 = LLL;
      elseif (Lobj<Hobj-epsi)
         a2 = HHH;
      else
         a2 = alpha2;
      end;
   end;

   if ( abs(a2-alpha2) < epsi*(a2 + alpha2 + epsi) )
      resultat = 0;
   end;
end;


if (resultat == 1)
   %%%
   a1 = alpha1 - (a2-alpha2);
   if (alpha1<0.0)
      a2 = a2 + a1;
      a1 = 0;
   elseif (a1>CCC)
      ttt = a1-CCC;
      a2 = a2 + ttt;
      a1 = CCC;
   end;


   %%%
   if (a1>0.0 && a1<CCC)
      bbb_nouv = bbb + E1 + (a1 - alpha1)*k11 + (a2 - alpha2)*k12;
   else
      if (a2>0.0 && a2<CCC)
         bbb_nouv = bbb + E2 + (a1 - alpha1)*k12 + (a2 - alpha2)*k22;
      else
         b1 = bbb + E1 + (a1 - alpha1)*k11 + (a2 - alpha2)*k12;
         b2 = bbb + E2 + (a1 - alpha1)*k12 + (a2 - alpha2)*k22;
 
         % méthode 1
         bbb_nouv = (b1+b2)/2.0;

         % méthode 2
         %pp1=rand(1);
         %pp2=rand(1);
         %bbb_nouv = (pp1*b1+pp2*b2)/(pp1+pp2);
      end;
   end;
   %%% v. p. 1451 de "Estimating, etc." ; cependant, rien d'intéressant ne se produit
   %%% ne se produit
   %epsibbb=1e-4;
   %if bbb_nouv>epsibbb;
   %   bbb_nouv=bbb_nouv-epsibbb;
   %end;
   
   delta_bbb = bbb_nouv - bbb;

   %fprintf(1,'%f %f -> %f\n',bbb,bbb_nouv,delta_bbb);
   %pause;


   %%%
   t1 = (a1-alpha1);
   t2 = (a2-alpha2);
   for ii=1:mmm
      if (alpha(ii)>0.0 && alpha(ii)<CCC)
         erreur_cache_nouv(ii) = erreur_cache_nouv(ii) + t1*noyau_fonc_1classe(i1,ii,xxx,sig2,path1) +  t2*noyau_fonc_1classe(i2,ii,xxx,sig2,path1) - delta_bbb;
      end;
   end;
   erreur_cache_nouv(i1) = 0;
   erreur_cache_nouv(i2) = 0;


   %%%
   alpha_nouv(i1) = a1;
   alpha_nouv(i2) = a2;
end;

