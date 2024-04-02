%%%%%%%%%%%%%%%%%%%%
% SVM 1 classe - avec noyau
%     0 : RBF gaussien
%     1 : "quadratique rationnel"
%     ...
%     De plus, des données à 'n' dimensions sont possibles
% Code Lagis et CentraleSupélec
%
% - Calcul du problème dual
% - Optimisation avec une méthode à point intérieur
%   (qui n'appartient ni au Lagis ni à CentraleSupélec !!!!!!!!!!!!!!!!!)
%
% Arguments :
% -----------
%
% Retours :
% ---------
%
% liste TODO
% ----------
%
% - 11/04/2006
%   Expliciter les arguments et les retours
%
% Stéphane Rossignol - 11/04/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

function [wwx,bbb,alpha] = lagis_svm_1classe_noyau_nd (xx,sig2,nu,use_octave,noyau)


%%% initialisations
ndimensions = size(xx,1);


%%% calcul de la solution : optimisation sous contraintes
tmatsca=clock;
mmm = size(xx,2);
for ii=1:mmm
  vec1=xx(:,ii);

  switch noyau
    case 0
      % noyau RBF gaussien
      matsca(ii,:) = lagis_rbf_gaussien(vec1, xx, sig2);
    case 1
      % noyau quadratique rationnel
      matsca(ii,:) = lagis_quadratique_rationnel(vec1, xx, sig2);
    case 2
      % noyau polynomial
      matsca(ii,:) = lagis_polynomial(vec1, xx, sig2);
    case 3
      % noyau multi-quadratique
      matsca(ii,:) = lagis_multi_quadratique(vec1, xx, sig2);
    otherwise
      fprintf(1,'Noyau inexistant\n');
  end;
end;
ematsca=etime(clock,tmatsca);
fprintf(1,'temps requis pour former H pour les SVM sur %d points : %f\n',mmm,ematsca);

%min(min(matsca))
%max(max(matsca))

c = zeros(mmm,1);
H = matsca;
A = ones(1,mmm);
b = 1;
%l = zeros(mmm,1); %'l'='c', donc pas la peine de former 'l'
u = ones(mmm,1)/(mmm*nu);

topt=clock;
[alpha,ytmp] = pr_loqo3(c, H, A, b, c, u, use_octave); % normalement : 'c' puis 'l', mais 'c=l'
eopt=etime(clock,topt);
fprintf(1,'temps requis pour optimiser les SVM sur %d points avec LOQO : %f\n',mmm,eopt);

alpha = alpha';

plotfig=0;
if plotfig==1
  figure(1);
  if use_octave
    clf;
  else
    clf;
  end;
  grid on;
  hold on;
  plot(alpha,'r');
  %zoom on; % "octave" se plaint
  hold off;
end;


%%% calcul de la solution : calcul de w
for ii=1:ndimensions
  wwx(ii) = sum(alpha.*xx(ii,:));
end;


%%% calcul de la solution : calcul de b
[valx,pos]=max(alpha);
switch noyau
  case 0
    % noyau RBF gaussien
    bbb = -alpha*lagis_rbf_gaussien(xx, xx(:,pos), sig2)';
  case 1
    % noyau quadratique rationnel
    bbb = -alpha*lagis_quadratique_rationnel(xx, xx(:,pos), sig2)';
  case 2
    % noyau polynomial
    bbb = -alpha*lagis_polynomial(xx, xx(:,pos), sig2)';
  case 3
    % noyau multi-quadratique
    bbb = -alpha*lagis_multi_quadratique(xx, xx(:,pos), sig2)';
  otherwise
    fprintf(1,'Noyau inexistant\n');
end;

