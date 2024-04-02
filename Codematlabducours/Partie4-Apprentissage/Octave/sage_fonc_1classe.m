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

% hypothèse : ici, on utilise toujours le noyau non linéaire (RBF gaussien)

function [sss] = sage_fonc_1classe(kk, alpha, bbb, xxx, sig2, mmm, path1)


%%% for the SVM computation, we need to add these "paths"
strversion=version;
if strcmp(version,'2.1.73')==1
  path(LOADPATH,path1);
else
  addpath(path1);
end;


sss = 0.0;


for ii=1:mmm
   if alpha(ii)>0
      sss = sss + alpha(ii)*noyau_fonc_1classe(ii,kk,xxx,sig2,path1);
      %alpha(ii)
      %noyau_fonc_1classe(ii,kk,xxx,sig2,path1)
      %sss
      %pause
   end;
end;

sss = sss - bbb;

