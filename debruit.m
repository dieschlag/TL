% Challenge de débruitage d'un signal.
%
% Les meilleures performances que j'obtiens, comme taux de débruitage, sont
% de 19. Essayez de faire mieux que moi. Je classerai les groupes par 
% performances, et plus les performances seront bonnes, plus vous 
% améliorerez votre note pour le cours de représentations parcimonieuses. 
% 20% de la note viendra des petits challenges comme celui-ci. Bien sûr,
% lors de votre présentation finale, il faudra expliquer votre méthode
% de débruitage (comment marche votre 'debruiteur.m').
%
% Vous n'avez rien à changer au code ci-dessous.
% Ecrivez une fonction 'debruiteur.m' qui reçoit en argument le signal à
% débruiter (fait par moi ci-dessous) dans le vecteur 'signal' ; et qui 
% renvoie le signal débruité dans le vecteur 'sigrec'. Je me charge 
% ensuite de mesurer les performances de votre débruiteur, utilisé à la 
% ligne 36 du code. Je testerai votre débruiteur groupe par groupe, avec
% ce code.
% 
% S. Rossignol -- 2024

clear all;
close all;

nll=100; %%% je fais des statistiques sur 1000 signaux
rati=zeros(nll,1);
for ll=1:nll
    %%% signal : un sinus tout seul avec du bruit
    fe=16000;
    mttt=1.5;
    ttt=[0:1/fe:mttt];
    sigpur = 0.5*cos(2*pi*500*ttt);
    signal = sigpur + 0.1*randn(1,length(sigpur));
    signal = signal/max(abs(signal));

    %%%%%%%%%%%% debruiteur est votre fonction
    [sigrec]=debruiteur(signal);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    erreur1=sum((sigpur-signal).^2)/sum(sigpur.^2);
    erreur2=sum((sigpur-sigrec).^2)/sum(sigpur.^2);

    rati(ll)=erreur1/erreur2;
    fprintf(1,'%d\n',ll)
end;
ratio=mean(rati);

figure(1);
clf;
plot(rati);
hold off;

fprintf(1,'taux de débruitage : %f\n', ratio);

function signal_debruite = debruiteur(signal_bruit)
    % Décomposition en ondelettes
    niveau_decomposition = 5; % Niveau de décomposition en ondelettes
    ondelette = 'db14'; % Type d'ondelette
    [c, l] = wavedec(signal_bruit, niveau_decomposition, ondelette);

    % Filtrage des coefficients d'ondelettes (seuils universels)
    thr = thselect(c,'rigrsure');
    wtthr = wthresh(c,"s",thr);
    


    % Reconstruction du signal débruité
    signal_debruite = waverec(wtthr, l, ondelette)
end

