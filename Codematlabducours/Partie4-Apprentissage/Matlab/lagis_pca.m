%%%%%%%%%%%%%%%%%%%%
% PCA
% Code Lagis
%
%
% Arguments :
% -----------
%
% - "data_center" doit compter 
%     * n_dimensions        lignes et 
%     * taille_echantillons colonnes
%     => par exemple, size(data_center)=> (2, 300)
%   et chaque ligne doit être centrée
% - "nvec" : nombre de vecteurs propres sur lesquels on projette
% - on est sous "octave" ou pas
%
% liste TODO
% ----------
%
% 1) 30/01/2006 -  28/03/2006 Résolu
%    Donner la possibilité de projeter sur plus de deux dimensions
%
%
% Stéphane Rossignol - 30/01/2006 ; et 2021
%
%%%%%%%%%%%%%%%%%%%%

function [projection] = lagis_pca (data_center,nvec,use_octave,use_mdavy)


%%% matrice de covariance
cov_matrix=(data_center*data_center')/length(data_center);


%%% valeurs propres et vecteurs propres
[eig_vec,eig_val] = eig(cov_matrix);
eig_val = eig_val*ones(length(cov_matrix),1);
if use_mdavy==0
  [eig_val,pos] = sort(eig_val,'descend');
else
  % parce que sur l'"octave" installé sur la machine de Manuel
  % l'argument MODE ('descend'/'ascend') n'existe pas
  [eig_valdavy,posdavy] = sort(eig_val);
  pos=posdavy(length(posdavy):-1:1);
  eigval=eig_valdavy(length(eig_valdavy):-1:1);
end;
%fprintf(1,"Valeurs propres 'eig' :\n");
%fprintf(1,"\t%f\n",eig_val);


%%% ce que l'on obtient avec svd (ajout du 22/09/06)
%[u, s, v] = svd(cov_matrix);
%size(cov_matrix)
%eig_val
%s


%%%
faitfig=0;
if faitfig==1
  figure(1)
  clf;
  grid on;
  hold on;
  plot(eig_val);
  title('valeurs propres -- PCA');
  %if use_octave==0
  %  zoom on;
  %end;
  hold off;
%  print -dpsc2 figures/donnees_lannion_valeursproprespca.ps
end;
%%eig_val/eig_val(1)
%%pause;


%%% trie des 'nvec' plus grandes valeurs propres
for ii=1:nvec
  vec = eig_vec(:,pos(ii));

  % projection 
  projection(ii,:) = vec'*data_center;
end;

