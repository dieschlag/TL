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

% hypothèse : ici, on utilise toujours le noyau non linéaire (RBF gaussien)

function [sss] = learned_func(kk, alpha, bbb, xxx, yyy, sig2)

mmm = length(alpha);
sss = 0.0;

for ii=1:mmm
   if (alpha(ii)>0)
      sss = sss + alpha(ii)*yyy(ii)*kernel_func(ii,kk,xxx,sig2);
   end;
end;

sss = sss - bbb;
