%On initialise le vecteur de temps et les 4 sinusoïdes
temps = linspace(0,2,16000);
sin1v = zeros(16000,1);
sin2v = zeros(16000,1);
sin3v = zeros(16000,1);
sin4v = zeros(16000,1);

fs = 8000;
for i = 1:16000
    sin1v(i) = 0.1*cos(880*temps(i)*2*pi+4*sin(2*pi*5*temps(i)));
    sin2v(i) = 0.2*cos(1760*temps(i)*2*pi+8*sin(2*pi*5*temps(i)));
    sin3v(i) = 0.3*cos(2640*temps(i)*2*pi+12*sin(2*pi*5*temps(i)));
    sin4v(i) = 0.4*cos(3520*temps(i)*2*pi+16*sin(2*pi*5*temps(i)));
end
sin1 = 0.1*cos(880*temps*2*pi);
sin2 = 0.2*cos(1760*temps*2*pi);
sin3 = 0.3*cos(2640*temps*2*pi);
sin4 = 0.4*cos(3520*temps*2*pi);
signal = sin1 + sin2 + sin4 + sin3; 
signalv = sin1v + sin2v + sin3v + sin4v;  
%spectrogram(signalv,320,80)
audiowrite('signal.wav',signal,8000); %On enregistre le signal dans un fichier .wav

sigfft = fft(signal);

% Normalisation de l'amplitude
sigfft = abs(sigfft) / length(signal);

% Fréquences correspondantes
freqs = linspace(0, fs, length(signal));

% Limiter les fréquences jusqu'à 4000 Hz
freqs = freqs(freqs <= 4000);
sigfft = sigfft(1:length(freqs));

% Affichage du spectre
figure;
plot(freqs, sigfft);
title('Spectre du signal');
xlabel('Fréquence (Hz)');
ylabel('Amplitude');

bbas=fir2(300,[0 0.5 0.5 1],[1 1 0 0]);
figure(2)
freqz(bbas,1,8000);
hold off



aa=1;
bhaut=fir2(300,[0 0.5 0.5 1],[0 0 1 1]);

figure(3); %%% réponse en fréquence du filtre
freqz(bhaut,1,8000);
hold off


sigfilb=filtfilt(bbas,aa,signalv);


sigfilh=filtfilt(bhaut,aa,signalv);

sousechb = zeros(8000,1);
for i = 1:8000
    sousechb(i) = sigfilb(2*i);
end



sousechh = zeros(8000,1);
for i = 1:8000
    sousechh(i) = sigfilh(2*i);
end



sigfilbb=filtfilt(bbas,aa,sousechb);


sousechbb = zeros(4000,1);
for i = 1:4000
    sousechbb(i) = sigfilbb(2*i);
end



sigfilbh = filtfilt(bhaut,aa,sousechb);


sousechbh = zeros(4000,1);
for i = 1:4000
    sousechbh(i) = sigfilbh(2*i);
end
fftsschbh = abs(fft(sousechbh));


sigfilhb =filtfilt(bbas,aa,sousechh); 



sousechhb = zeros(4000,1);
for i = 1:4000
    sousechhb(i) = sigfilhb(2*i);
end
fftsschhb = abs(fft(sousechhb));



sigfilhh =filtfilt(bhaut,aa,sousechh); 

sousechhh = zeros(4000,1);
for i = 1:4000
    sousechhh(i) = sigfilhh(2*i);
end
fftsschhh = abs(fft(sousechhh));

% freqbb = freqhil(sousechbb,2000,1400);
% disp(freqbb)
% 
% freqbh = freqhil(sousechbh,2000,1400);
% disp(freqbh)
% 
% freqhb = freqhil(sousechhb,2000,1400);
% disp(freqhb)
% 
% freqhh = freqhil(sousechhh,2000,1400);
% disp(freqhh)
tabfbb = freqhil(sousechbb,fs/4);
tabfbh = freqhil(sousechbh,fs/4);
tabfhb = freqhil(sousechhb,fs/4);
tabfhh = freqhil(sousechhh,fs/4);
figure(4)
plot(temps(101:2000)/2000,tabfbb(101:2000));
hold on
plot(temps(101:2000)/2000,tabfbh(101:2000));
hold on
plot(temps(101:2000)/2000,tabfhb(101:2000))
hold on
plot(temps(101:2000)/2000,tabfhh(101:2000))

title("Evolution de la fréquence au cours du temps avec un vibrato d'amplitude 20Hz")
xlabel("Temps en secondes")
ylabel("Fréquence en Hz")
hold off


%--------------------------------Fonctions---------------------------------
function plot_fft(tab,titre) % Fonction utilisée pour tracer les différents graphes
fftab = abs(fft(tab));
plot(fftab);
title(titre)
end


function freq = freqhil(tab,fs)
%trouve la fréquence du signal représenté par un tableau, en utilisant la
%transformation de Hilbert
hil = hilbert(tab);
phases = unwrap(angle(hil));
freqinst = zeros(length(phases)-1);
for ii=1:length(phases)-1
    freqinst(ii) = (phases(ii+1)-phases(ii))/2/pi*fs;
end

    freq=freqinst;
end