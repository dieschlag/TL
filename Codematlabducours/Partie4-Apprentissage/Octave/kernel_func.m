%%%%%%%%%%%%%%%%%%%%
% Optimisation avec la méthode SMO - pour SVM 2 classes
% Code Lagis et CentraleSupélec
%
% liste TODO
% ----------
%
% - 25/04/06
%   Ce n'est pas très souple !!!!!!!!!!!!!!!!!!
%
% Stéphane Rossignol - 25/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

function [ki1i2] = kernel_func(i1, i2, xxx, sig2)

ki1i2 = lagis_rbf_gaussien(xxx(:,i1), xxx(:,i2), sig2);
