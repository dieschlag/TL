% Paramètres
num_birds = 3; % Nombre d'oiseaux
num_recordings = 1000; % Nombre d'enregistrements par oiseau dans SONS
num_vc_recordings = 500; % Nombre d'enregistrements par oiseau dans SONS-VC
duration = 0.05; % Durée de chaque trille en secondes (50 ms)

% Initialisation des tableaux de caractéristiques
features_sons = zeros(num_birds * num_recordings, 2); % Pour SONS
features_vc = zeros(num_birds * num_vc_recordings, 2); % Pour SONS-VC

% Lecture des fichiers audio dans SONS
for bird = 1:num_birds
    for recording = 1:num_recordings
        filename = fullfile('SONS', ['oiseau', num2str(bird)], [num2str(recording), '.wav']);
        [signal, fs] = audioread(filename);
        
        % Extraction des caractéristiques (fréquence et amplitude de la fondamentale)
        [freq, ampl] = extract_features(signal, fs, duration);
        
        % Stockage des caractéristiques dans le tableau approprié
        index = (bird - 1) * num_recordings + recording;
        features_sons(index, :) = [freq, ampl];
    end
end

% Lecture des fichiers audio dans SONS-VC
for bird = 1:num_birds
    for recording = 1:num_vc_recordings
        filename = fullfile('SONS-VC', ['oiseau', num2str(bird)], [num2str(recording), '.wav']);
        [signal, fs] = audioread(filename);
        
        % Extraction des caractéristiques (fréquence et amplitude de la fondamentale)
        [freq, ampl] = extract_features(signal, fs, duration);
        
        % Stockage des caractéristiques dans le tableau approprié
        index = (bird - 1) * num_vc_recordings + recording;
        features_vc(index, :) = [freq, ampl];
    end
end

% Affichage des caractéristiques pour vérification
disp('Caractéristiques extraites pour SONS :');
disp(features_sons);

disp('Caractéristiques extraites pour SONS-VC :');
disp(features_vc);

% Fonction pour extraire les caractéristiques (fréquence et amplitude de la fondamentale)
function [freq, ampl] = extract_features(signal, fs, duration)
    % Ici vous pouvez mettre votre code d'extraction des caractéristiques, par exemple, l'analyse spectrale.
    % Pour simplifier, je vais générer des valeurs aléatoires pour l'exemple.
    freq = rand(); % Fréquence aléatoire
    ampl = rand(); % Amplitude aléatoire
end