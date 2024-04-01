
%%Question 1

% Paramètres du signal
fundamental_freq = 880; % Fréquence fondamentale en Hz
sampling_freq = 8000; % Fréquence d'échantillonnage en Hz
duration = 1; % Durée du signal en secondes

% Générer le signal harmonique composé de 4 sinus
t = linspace(0, duration, sampling_freq * duration);
signal = sin(2 * pi * fundamental_freq * t) + sin(2 * pi * 2 * fundamental_freq * t) + ...
         sin(2 * pi * 3 * fundamental_freq * t) + sin(2 * pi * 4 * fundamental_freq * t);

% Écrire le signal dans un fichier .wav
filename = 'signal_harmonique.wav';
audiowrite(filename, signal, sampling_freq);

% Lire le fichier .wav pour vérifier le son
% [y, fs] = audioread(filename);
% sound(y, fs);

% Analyse du signal
% Vous pouvez utiliser la transformée de Fourier pour analyser les fréquences du signal
% par exemple : 
N = length(signal);
frequencies = (0:N-1) * (sampling_freq / N);
signal_fft = abs(fft(signal));
plot(frequencies, signal_fft);
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('Transformée de Fourier du signal');

%% Question 2

% Paramètres du vibrato
vibrato_amplitude = 20; % Amplitude du vibrato en Hz
vibrato_frequency = 5; % Fréquence du vibrato en Hz

% Générer l'onde sinusoïdale pour le vibrato
vibrato_wave = vibrato_amplitude * sin(2 * pi * vibrato_frequency * t);

% Ajouter le vibrato au signal de départ
vibrato_signal = signal .* (1 + vibrato_wave);

% Écrire le signal avec vibrato dans un fichier .wav
vibrato_filename = 'signal_harmonique_vibrato.wav';
audiowrite(vibrato_filename, vibrato_signal, sampling_freq);

% Afficher le spectrogramme pour observer le vibrato
spectrogram(vibrato_signal, 256, 250, [], sampling_freq, 'yaxis');
title('Spectrogramme du signal avec vibrato');

% Afficher le spectre avec pspectrum
pspectrum(vibrato_signal, sampling_freq, 'spectrogram', 'OverlapPercent', 0, 'FrequencyLimits', [0, 2000]);
title('Spectrogramme du signal avec vibrato');