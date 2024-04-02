%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la méthode SMO - pour SVM 2 classes
% Code Lagis et CentraleSupélec
%
% liste TODO
% ----------
%
% Stéphane Rossignol - 25/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

function [resultat, alpha_new, bbb_new, error_cache_new] = takeStep (i1, i2, alpha, xxx, yyy, CCC, error_cache, bbb, sig2)


%%%
alpha_new = alpha;
bbb_new = bbb;
error_cache_new = error_cache;
resultat = 1;

mmm = length(alpha);


%%%
epsi = 0.001;


%%%
if (i1==i2)
   resultat=0;
end;


if (resultat==1)
   %%%
   alpha1 = alpha(i1);
   y1 = yyy(i1);
   if (alpha1>0.0 && alpha1<CCC)
      E1 = error_cache(i1);
   else
      E1 = learned_func(i1, alpha, bbb, xxx, yyy, sig2) - y1;
   end;

   alpha2 = alpha(i2);
   y2 = yyy(i2);
   if (alpha2>0.0 && alpha2<CCC)
      E2 = error_cache(i2);
   else
      E2 = learned_func(i2, alpha, bbb, xxx, yyy, sig2) - y2;
   end;

   sss = y1*y2;


   %%%
   if (y1==y2)
      gamma = alpha1+alpha2;
      if (gamma>CCC)
         LLL = gamma-CCC;
         HHH = CCC;
      else
         LLL = 0.0;
         HHH =gamma;
      end;
   else
      gamma = alpha1-alpha2;
      if (gamma>0.0)
         LLL = 0.0;
         HHH = CCC-gamma;
      else
         LLL = -gamma;
         HHH = CCC;
      end;
   end;

   if (LLL==HHH)
      resultat=0;
   end;
end;



if (resultat == 1)
   %%%
   k11 = kernel_func(i1, i1, xxx, sig2);
   k12 = kernel_func(i1, i2, xxx, sig2);
   k22 = kernel_func(i2, i2, xxx, sig2);
   etamine = 2.0*k12 - k11 - k22;


   %%%
   if (etamine<0.0)
      a2 = alpha2 + y2*(E2-E1)/etamine;
      if (a2<LLL)
         a2 = LLL;
      elseif a2>HHH
         a2 = HHH;
      end;
   else
      c1 = etamine/2.0;
      c2 = y2*(E1-E2) - etamine*alpha2;
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
   a1 = alpha1 - sss*(a2-alpha2);
   if (alpha1<0.0)
      a2 = a2 + sss*a1;
      a1 = 0;
   elseif (a1>CCC)
      ttt = a1-CCC;
      a2 = a2 + sss*ttt;
      a1 = CCC;
   end;


   %%%
   if (a1>0.0 && a1<CCC)
      bbb_new = bbb + E1 + y1*(a1 - alpha1)*k11 + y2*(a2 - alpha2)*k12;
   else
      if (a2>0.0 && a2<CCC)
         bbb_new = bbb + E2 + y1*(a1 - alpha1)*k12 + y2*(a2 - alpha2)*k22;
      else
         b1 = bbb + E1 + y1*(a1 - alpha1)*k11 + y2*(a2 - alpha2)*k12;
         b2 = bbb + E2 + y1*(a1 - alpha1)*k12 + y2*(a2 - alpha2)*k22;
         bbb_new = (b1+b2)/2.0;
      end;
   end;
   delta_bbb = bbb_new - bbb;


   %%%
   t1 = y1*(a1-alpha1);
   t2 = y2*(a2-alpha2);
   for ii=1:mmm
      if (alpha(ii)>0.0 && alpha(ii)<CCC)
         error_cache_new(ii) = error_cache_new(ii) + t1*kernel_func(i1,ii,xxx,sig2) + t2*kernel_func(i2,ii,xxx,sig2) - delta_bbb;
      end;
   end;
   error_cache_new(i1) = 0;
   error_cache_new(i2) = 0;


   %%%
   alpha_new(i1) = a1;
   alpha_new(i2) = a2;
end;

