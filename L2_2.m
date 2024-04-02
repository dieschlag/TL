close all

% Paramètres du signal
fe = 16000; % Fréquence d'échantillonnage en Hz
duree = 1; % Durée du signal en secondes
t = 0:1/fe:duree-1/fe; % Vecteur temps

% Génération des signaux sinusoidaux
f1 = 100; % Fréquence du premier sinus en Hz
f2 = 500; % Fréquence du deuxième sinus en Hz
f3 = 1000;
signal1 = cos(2*pi*f1*t); % Premier cosinus
signal2 = cos(2*pi*f2*t); % Deuxième cosinus
signal3 = cos(2*pi*f3*t)

% Somme des deux signaux sinusoidaux
signal = signal1 + signal2 + signal3;

% Calcul du spectre du signal
% Calcul de la transformée de Fourier
Y = fft(signal);
L = length(signal);
f_signal = (0:L-1)*(fe/L);

% Calcul de l'amplitude du spectre (module de la transformée de Fourier)
Pxx_signal = abs(Y/L);

% Tracer le spectre avec l'axe des abscisses représentant la fréquence
figure;
plot(f_signal, Pxx_signal, 'r', 'LineWidth', 1); % Spectre en rouge
title('Spectre du signal originel');
xlabel('Fréquence (Hz)');
ylabel('Amplitude');

level = 9;
wv = "db20";
[C,L] = wavedec(signal,level,wv);

mra = zeros(level+1,numel(signal));
for k=1:level
    mra(k,:) = wrcoef("d",C,L,wv,k);
end
figure;
for k=1:level
    nexttile
    plot(mra(k,101:612))
    title("Projection Onto Detail Subspace "+num2str(k))
end

mra(end,:) = wrcoef("a",C,L,wv,level);
mraSum = sum(mra,1);
max(abs(mraSum-signal));

% Nom du dossier
nom_dossier = 'Sep';

% Chemin complet du dossier
chemin_dossier = fullfile(pwd, nom_dossier);

% Vérifier si le dossier existe et le créer s'il n'existe pas
if ~isfolder(chemin_dossier)
    mkdir(chemin_dossier);
else
    % Vider le contenu du dossier s'il existe déjà
    delete(fullfile(chemin_dossier, '*'));
end



dossier_sauvegarde = 'Sep'; % Nom du dossier de sauvegarde
chemin_dossier = fullfile(pwd, dossier_sauvegarde); % Chemin complet du dossier de sauvegarde
liste = [];

audiowrite(fullfile(chemin_dossier, 'signal_original.wav'), signal, fe);
% Vérifier si les coefficients de la matrice sur une ligne décrivent un signal d'amplitude 1
for i = 1:size(mra, 1)
    i
    if max(abs(mra(i, :)) - 1) <= 0.1
        % Mettre à zéro les coefficients des lignes ne vérifiant pas ce critère
        mra(i, :) = 0;
    else 
        % Retrancher les coefficients à la variable signal
        liste = [liste, i]
        signal = signal - mra(i,:);
        audiowrite(fullfile(chemin_dossier, ['ondelette', num2str(i), '.wav']), mra(i,:), fe);

    end
end

% Enregistrer la variable signal au format WAV
audiowrite(fullfile(chemin_dossier, 'signal.wav'), signal, fe);

figure;
tiledlayout(level+1,1)


for k=1:level
    nexttile
    plot(mra(k,101:612))
    title("Projection Onto Detail Subspace "+num2str(k))
end

% Affichage du spectre du signal "signal" sur une nouvelle figure
figure;
[Pxx_signal, f_signal] = periodogram(signal, [], [], fe);
plot(f_signal, 10*log10(Pxx_signal), 'r', 'LineWidth', 2); % Spectre en rouge
title('Spectre du signal');
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');

% Calcul de la transformée de Fourier
Y = fft(signal);
L = length(signal);
f_signal = (0:L-1)*(fe/L);

% Calcul de l'amplitude du spectre (module de la transformée de Fourier)
Pxx_signal = abs(Y/L);

% Tracer le spectre avec l'axe des abscisses représentant la fréquence
figure;
plot(f_signal, Pxx_signal, 'r', 'LineWidth', 1); % Spectre en rouge
title('Spectre du signal après suppression');
xlabel('Fréquence (Hz)');
ylabel('Amplitude');


% Boucle pour afficher les spectres des signaux mra(i,:) pour chaque élément de la liste "liste"
for i = liste
    % Calcul de la transformée de Fourier
    Y = fft(mra(i,:));
    L = length(signal);
    f_signal = (0:L-1)*(fe/L);
    
    % Calcul de l'amplitude du spectre (module de la transformée de Fourier)
    Pxx_signal = abs(Y/L);
    
    % Tracer le spectre avec l'axe des abscisses représentant la fréquence
    figure;
    plot(f_signal, Pxx_signal, 'r', 'LineWidth', 1); % Spectre en rouge
    title(['Spectre de mra(', num2str(i), ')']);
    xlabel('Fréquence (Hz)');
    ylabel('Amplitude');

end




% nexttile
% plot(mra(end,101:612))
% title("Projection Onto Approximation Subspace")
% seuil = thselect(signal,'rigrsure');
% nbSignaux = 0;



% % Décomposition en ondelettes
% niveau_decomposition = 9; % Niveau de décomposition
% ondelette = 'db20'; % Type d'ondelette
% [c, l] = wavedec(signal, niveau_decomposition, ondelette);
% 
% % Représentation graphique des coefficients d'ondelettes
% figure;
% for i = 1:niveau_decomposition
%     subplot(niveau_decomposition,1,i);
%     plot(c(l(i)+1:l(i+1))); % Affichage des coefficients pour chaque niveau de décomposition
%     title(['Niveau de décomposition ', num2str(i)]);
%     xlabel('Indice du coefficient');
%     ylabel('Valeur du coefficient');
% end

% % Séparation des composantes basses et hautes fréquences
% division = floor(length(l)/2)
% c_basses = c(1:division); % Garder les coefficients de basses fréquences intacts
% c_hautes = c(division+1:end); % Mettre à zéro les coefficients de hautes fréquences
% 
% % Reconstruction des signaux séparés
% signal_basses = waverec(c_basses, l(1:division), ondelette);
% signal_hautes = waverec(c_hautes, l(division+1:end), ondelette);

% % Représentation graphique des signaux reconstruits
% figure;
% subplot(3,1,1);
% plot(t, signal);
% title('Signal d''origine');
% xlabel('Temps (s)');
% ylabel('Amplitude');
% subplot(3,1,2);
% plot(t, signal_basses);
% title('Composante de basses fréquences');
% xlabel('Temps (s)');
% ylabel('Amplitude');
% subplot(3,1,3);
% plot(t, signal_hautes);
% title('Composante de hautes fréquences');
% xlabel('Temps (s)');
% ylabel('Amplitude');
% 
% % Calcul de la transformée de Fourier des signaux
% spectre_basses = abs(fft(signal_basses));
% spectre_hautes = abs(fft(signal_hautes));
% 
% % Fréquences correspondant aux composantes du spectre
% frequence = (0:length(signal_basses)-1) * (fe / length(signal_basses));
% 
% % Tracé du spectre
% figure;
% subplot(2,1,1);
% plot(frequence, spectre_basses);
% xlabel('Fréquence (Hz)');
% ylabel('Amplitude');
% title('Spectre de la composante de basses fréquences');
% 
% subplot(2,1,2);
% plot(frequence, spectre_hautes);
% xlabel('Fréquence (Hz)');
% ylabel('Amplitude');
% title('Spectre de la composante de hautes fréquences');