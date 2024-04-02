
close all



% Paramètres du signal
fundamental_freq = 200; % Fréquence fondamentale en Hz
sampling_freq = 16000; % Fréquence d'échantillonnage en Hz
duration = 1; % Durée du signal en secondes
amplitude = 1; % Amplitude du signal sinusoïdal

% Générer le temps
t = linspace(0, duration, sampling_freq * duration);

% Générer le signal sinusoïdal
sinus = amplitude * sin(2 * pi * fundamental_freq * t);

% Générer le bruit blanc
bruit = 0.1*randn(1, length(t));

% Ajouter le bruit blanc au signal sinusoïdal
signal_bruit = sinus + bruit;

% Décomposition en ondelettes
niveau_decomposition = 7; % Niveau de décomposition en ondelettes
ondelette = 'coif4'; % Type d'ondelette
[c, l] = wavedec(signal_bruit, niveau_decomposition, ondelette);

% Filtrage des coefficients d'ondelettes (seuils universels)
seuil = sqrt(2 * log(length(signal_bruit)));
c_filtre = wthresh(c, 'h', seuil);

% Reconstruction du signal débruité
signal_debruite = waverec(c_filtre, l, ondelette);

% Calcul de la puissance du bruit à partir du signal bruité
puissance_bruit_avant_debruitage = mean((signal_bruit - sinus).^2);

% Calcul de la puissance du bruit à partir du signal débruité
puissance_bruit_apres_debruitage = mean((signal_debruite - sinus).^2);

% Calcul de la puissance du signal à partir du signal bruité
puissance_signal_avant_debruitage = mean(sinus.^2);

% Calcul de la puissance du signal à partir du signal débruité
puissance_signal_apres_debruitage = mean((signal_debruite - sinus).^2 + sinus.^2);

% Calcul de la RSB avant et après le débruitage (en dB)
rsb_avant_debruitage = 10 * log10(puissance_signal_avant_debruitage / puissance_bruit_avant_debruitage);
rsb_apres_debruitage = 10 * log10(puissance_signal_apres_debruitage / puissance_bruit_apres_debruitage);

% Affichage des résultats
fprintf('RSB avant débruitage : %.2f dB\n', rsb_avant_debruitage);
fprintf('RSB après débruitage : %.2f dB\n', rsb_apres_debruitage);


% Afficher le signal original, le signal bruité et le signal débruité
figure;
subplot(3,1,1);
plot(t, sinus, 'b');
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal sinusoïdal pur');

subplot(3,1,2);
plot(t, signal_bruit, 'r');
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal bruité');

subplot(3,1,3);
plot(t, signal_debruite, 'g');
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal débruité');


% Paramètres du dossier de sauvegarde
dossier_sauvegarde = 'TL2'; % Nom du dossier de sauvegarde
chemin_dossier = fullfile(pwd, dossier_sauvegarde); % Chemin complet du dossier de sauvegarde

% Vérifier si le dossier existe et le créer s'il n'existe pas
if ~isfolder(chemin_dossier)
    mkdir(chemin_dossier);
else
    % Vider le contenu du dossier s'il existe déjà
    rmdir(chemin_dossier, 's');
    mkdir(chemin_dossier);
end

% Enregistrer les signaux dans le dossier de sauvegarde
audiowrite(fullfile(chemin_dossier, 'signal_sans_bruit.wav'), sinus, sampling_freq);
audiowrite(fullfile(chemin_dossier, 'signal_bruit.wav'), signal_bruit, sampling_freq);
audiowrite(fullfile(chemin_dossier, 'signal_debruite.wav'), signal_debruite, sampling_freq);

function signal_debruite = debruiteur(signal_bruit)
    bonjour = 'bonjour'
    % Décomposition en ondelettes
    niveau_decomposition = 7; % Niveau de décomposition en ondelettes
    ondelette = 'vade'; % Type d'ondelette
    [c, l] = wavedec(signal_bruit, niveau_decomposition, ondelette);

    % Filtrage des coefficients d'ondelettes (seuils universels)
    thr = thselect(c,'rigrsure');
    wtthr = wthresh(wt,"s",thr);
    


    % Reconstruction du signal débruité
    signal_debruite = waverec(wtthr, l, ondelette)
end



