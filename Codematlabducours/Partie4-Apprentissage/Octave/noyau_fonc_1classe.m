%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la m\'ethode SMO - pour SVM 1 classe
% Code Lagis et CentraleSupélec
%
% A priori, même chose que pour SMO 2 classes
%
% liste TODO
% ----------
%
% - 27/04/06
%   Ce n'est pas très souple !!!!!!!!!!!!!!!!!!
%
% St\'ephane Rossignol - 25/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%


function [ki1i2] = noyau_fonc_1classe(i1, i2, xxx, sig2, path1)


%%% for the SVM computation, we need to add these "paths"
strversion=version;
if strcmp(version,'2.1.73')==1
  path(LOADPATH,path1);
else
  addpath(path1);
end;


ki1i2 = lagis_rbf_gaussien(xxx(:,i1), xxx(:,i2), sig2);
