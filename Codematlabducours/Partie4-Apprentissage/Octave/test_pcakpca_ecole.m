%%%%%%%%%%%%%%%%%%%%
%
% PCA   : cas d'école pour l'état de l'art et notre compréhension
%         personnelle de l'outil
%
% K-PCA : même chose que ci-dessus
%
% PPR   : même chose que ci-dessus (26/09/2006) => c'est la K-PCA de Pad\'e-Rayleigh-Ritz, non fournie en 2021
%
% Stéphane Rossignol -- 23/01/2006 ; et 2021
%
%%%%%%%%%%%%%%%%%%%%


clear all;
close all;


% initialisations
faitfig=1; % converti les figures en .ps (pour les rapports, etc.)
utilise_kpca_c = 0; %%% K-PCA reprogrammée en C => non fournie pour la ST7 de CS

%nnn =3700; % nombre d'échantillons par classe : attention ne peut pas être trop grand !!!
%nnn =1000;
nnn = 200;
s_c = 0.2; % std (commun à toutes les classes et toutes
           % les dimensions, pour simplifier)

% classe 1
s1_x = s_c; s1_y = s_c; s1_z = s_c;
m1_x = 0.0; m1_y = 0.0; m1_z = 0.0;
cl1_x = s1_x*randn(1,nnn) + m1_x;
cl1_y = s1_y*randn(1,nnn) + m1_y;
cl1_z = s1_z*randn(1,nnn) + m1_z;

% classe 2
s2_x = s_c; s2_y = s_c; s2_z = s_c;
m2_x = 2.0; m2_y = 2.0; m2_z = 2.0;
cl2_x = s2_x*randn(1,nnn) + m2_x;
cl2_y = s1_y*randn(1,nnn) + m2_y;
cl2_z = s1_z*randn(1,nnn) + m2_z;

% classe 3
s3_x = s_c; s3_y = s_c; s3_z = s_c;
m3_x = 2.0; m3_y = 0.0; m3_z = 2.0;
cl3_x = s3_x*randn(1,nnn) + m3_x;
cl3_y = s3_y*randn(1,nnn) + m3_y;
cl3_z = s3_z*randn(1,nnn) + m3_z;


% plot pour la route
ouiplot=1;
if ouiplot==1
  figure(1);
  clf;
  plot3(cl1_x,cl1_y,cl1_z,'ro');
  hold on;
  plot3(cl2_x,cl2_y,cl2_z,'b*');
  plot3(cl3_x,cl3_y,cl3_z,'g+');
  xlabel('dimension 1');
  ylabel('dimension 2');
  zlabel('dimension 3');
  title('donnees dans espace-observation original');
  grid on;
  %%zoom on;
  hold off;
  if faitfig==1
    print -dpsc2 donnees_initiales.ps
  end;
end;


% matrice de données formée
dd_brutes = zeros(3,3*nnn);
dd_brutes(1,:) = [cl1_x cl2_x cl3_x];
dd_brutes(2,:) = [cl1_y cl2_y cl3_y];
dd_brutes(3,:) = [cl1_z cl2_z cl3_z];
clear cl1_x cl1_y cl1_z; % (GM+ -> voir lagis_kpca.m)
clear cl2_x cl2_y cl2_z; % (GM+ -> voir lagis_kpca.m)
clear cl3_x cl3_y cl3_z; % (GM+ -> voir lagis_kpca.m)


% donn\'ees centrées
dd(1,:) = dd_brutes(1,:) - mean(dd_brutes(1,:));
dd(2,:) = dd_brutes(2,:) - mean(dd_brutes(2,:));
dd(3,:) = dd_brutes(3,:) - mean(dd_brutes(3,:));
clear dd_brutes; % (GM+ -> voir lagis_kpca.m)


%%%
% PCA

projection = lagis_pca(dd,2,1,0);


% plot pour la route
if ouiplot==1
  figure(2);
  clf;
  grid on;
  hold on;
  plot(projection(1,0*nnn+1:1*nnn),projection(2,0*nnn+1:1*nnn),'ro');
  plot(projection(1,1*nnn+1:2*nnn),projection(2,1*nnn+1:2*nnn),'b*');
  plot(projection(1,2*nnn+1:3*nnn),projection(2,2*nnn+1:3*nnn),'g+');
  xlabel('dimension 1');
  ylabel('dimension 2');
  title('donnees dans espace-observation apres projection PCA');
  %%zoom on;
  hold off;
  if faitfig==1
    print -dpsc2 donnees_projetees_PCA.ps
  end;
end;
clear cov_matrix eig_val eig_vec pos projection; % (GM+ -> voir lagis_kpca.m)


%%%
% K-PCA

% estimation/choix du sig2 pour RBF_kernel
sig2 = 0.7; %0.7
%%sig2 = lagis_sig2(dd); % ATTENTION : extrêmement gourmant en temps
%%fprintf(1,'sig2=%f sig2_empirique=~%f\n',sig2,0.7); %sig_empirique est faux


%% K-PCA :
%  - avec LS_SVM, utiliser "kpca" ; avec le code du Lagis, utiliser "lagis_kpca"
%% - "lagis_kpca_prr" a été ajout\'e le 26/09/06
%for sig2=0.1:0.01:1.0
  %%%
  % calcul de la KPCA classique

  if utilise_kpca_c==0  %%%%%%%%%%%%% 1. KPCA en matlab/octave
    tt=cputime;

    %[eigval_kpca,eigvec_kpca] = kpca(dd','RBF_kernel',sig2);
    %[eigval_kpca,eigvec_kpca] = kpca(dd','lin_kernel',[]);

    [eigvec_kpca,eigval_kpca] = lagis_kpca(dd,sig2,'RBFgaussien',1);
    %[eigvec_kpca,eigval_kpca] = lagis_kpca(dd,2.0,'polynomial',1);
    ee=cputime-tt;
  else                  %%%%%%%%%%%%% 2. KPCA en matlab/octave => pas disponible pour la ST7 de CS
    % écriture des données
    ndimensions=3;
    if (3*nnn>10)
      npropres=10; % nombre de valeurs/vecteurs propres sauvés par le programme C
    else
      npropres=3*nnn;
    end;

    %versionpg = "_v1";
    versionpg = "_v2";

    ff=fopen(["fichiers_donnees/donnees_kpca_c" versionpg ".txt"],"wt");
    fprintf(ff,"%d ",3*nnn);
    fprintf(ff,"%d ",ndimensions);
    fprintf(ff,"%f ",sig2);
    fprintf(ff,"%d ",npropres);
    for ii=1:ndimensions
      fprintf(ff,"%f ",dd(ii,:));
    end;
    fclose(ff);

    % lancement du programme C
    %%tt=cputime;
    tt=clock;
    system(["lagis_kpca" versionpg]);
    ee=etime(clock,tt);
    %%ee=cputime-tt; %% ee et tt obtenus avec cputime ne servent \`a rien

    % lecture des données
    if strcmp(versionpg,"_v1")==1
      load fichiers_donnees/resultats_kpca_c_v1.txt;
      eigval_kpca = resultats_kpca_c_v1(1,:);
      eigvec_kpca = resultats_kpca_c_v1(2:3*nnn+1,:);
    elseif strcmp(versionpg,"_v2")==1
      load fichiers_donnees/resultats_kpca_c_v2.txt;
      eigval_kpca = resultats_kpca_c_v2(1,:);
      eigvec_kpca = resultats_kpca_c_v2(2:3*nnn+1,:);
    end;
  end;
  fprintf(1,'temps requis pour calculer la K-PCA : nnn=%d t=%f\n',nnn,ee);


  % "projection" sur les deux (et trois) vecteurs propres principaux
  [aaa,pos]=sort(eigval_kpca,'descend');
  vec1 = eigvec_kpca(:,pos(1));
  vec2 = eigvec_kpca(:,pos(2));
  vec3 = eigvec_kpca(:,pos(3));
  proj2d=[vec1 vec2]';
  proj3d=[vec1 vec2 vec3]';


  % essai de détermination du meilleur sig2, basé => non disponible en 2021
  % sur les convex hull
  % => est-ce utile ?
  %convex_hull1 = convhull(proj2d(1,0*nnn+1:1*nnn),proj2d(2,0*nnn+1:1*nnn));
  %cx1 = proj2d(1,convex_hull1);
  %convex_hull2 = convhull(proj2d(1,1*nnn+1:2*nnn),proj2d(2,1*nnn+1:2*nnn));
  %convex_hull3 = convhull(proj2d(1,2*nnn+1:3*nnn),proj2d(3,2*nnn+1:3*nnn));


  % plot pour la route
  figure(3);
  clf;
  grid on;
  hold on;
  plot(proj2d(1,0*nnn+1:1*nnn),proj2d(2,0*nnn+1:1*nnn),'or');
  plot(proj2d(1,1*nnn+1:2*nnn),proj2d(2,1*nnn+1:2*nnn),'*b');
  plot(proj2d(1,2*nnn+1:3*nnn),proj2d(2,2*nnn+1:3*nnn),'+g');
  xlabel('dimension 1');
  ylabel('dimension 2');
  title('donnees dans espace-observation apres projection K-PCA 2D');
  %%zoom on;
  hold off;
  print("figures/donnees_projetees_KPCA.ps","-dpsc2");

  prj3d=0;    %%% la 3D n'apporte rien, a priori
  if prj3d==1
    figure(4);
    clf;
    plot3(proj3d(1,0*nnn+1:1*nnn),proj3d(2,0*nnn+1:1*nnn),proj3d(3,0*nnn+1:1*nnn),'or');
    hold on;
    plot3(proj3d(1,1*nnn+1:2*nnn),proj3d(2,1*nnn+1:2*nnn),proj3d(3,1*nnn+1:2*nnn),'*b');
    plot3(proj3d(1,2*nnn+1:3*nnn),proj3d(2,2*nnn+1:3*nnn),proj3d(3,2*nnn+1:3*nnn),'+g');
    xlabel('dimension 1');
    ylabel('dimension 2');
    zlabel('dimension 3');
    title('donnees dans espace-observation apres projection K-PCA 3D');
    grid on;
    %%zoom on;
    hold off;
  end;


  %%%
  % calcul de la KPCA PRR
  % => elle a besoin de la KPCA classique lors de l'étape de validation
  do_prr=0;
  if do_prr==1
    nproj=2; % nbre de dimensions où l'on projette
    tt=cputime;
    [eigvec_kpca_prr,eigval_kpca_prr] = lagis_kpca_prr(dd,sig2,'RBFgaussien',nproj,aaa);
%    [eigvec_kpca,eigval_kpca] = lagis_kpca_prr(dd,2.0,'polynomial',nproj);

    ee=cputime-tt;
    fprintf(1,'temps requis pour calculer la K-PCA PRR : %f\n',ee);


    % "projection" sur les deux vecteurs propres principaux
    [aaa_prr,pos_prr]=sort(eigval_kpca_prr,'descend');
    vec1_prr = eigvec_kpca_prr(:,pos_prr(1));
    vec2_prr = eigvec_kpca_prr(:,pos_prr(2));
    proj2d_prr=[vec1_prr vec2_prr]';


    % plot pour la route
    figure(5);
    clf;
    grid on;
    hold on;
    plot(proj2d_prr(1,0*nnn+1:1*nnn),proj2d_prr(2,0*nnn+1:1*nnn),'or');
    plot(proj2d_prr(1,1*nnn+1:2*nnn),proj2d_prr(2,1*nnn+1:2*nnn),'*b');
    plot(proj2d_prr(1,2*nnn+1:3*nnn),proj2d_prr(2,2*nnn+1:3*nnn),'+g');
    xlabel('dimension 1');
    ylabel('dimension 2');
    title('donnees dans espace-observation apres projection K-PCA PRR 2D');
    %%zoom on;
    hold off;
    print("figures/donnees_projetees_KPCA_PRR.ps","-dpsc2");


    %%%
    fprintf(1,"\n\n\n");
    fprintf(1,"\nValeurs propres classiques :\n");
    fprintf(1,"%.3f ",aaa(1:length(aaa_prr)));

    fprintf(1,"\nValeurs propres Pad\'e-Rayleigh-Ritz :\n");
    fprintf(1,"%.3f ",aaa_prr);
    fprintf(1,"\n\n\n");
  end;


%  fprintf(1,'sig2=%f\n',sig2);
%  pause(1);
%end;
