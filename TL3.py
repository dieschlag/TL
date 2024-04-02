import os
import numpy as np
from scipy.io import wavfile
from scipy.signal import find_peaks
from sklearn.neighbors import KNeighborsClassifier
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

def extract_features(file_path):
    sample_rate, data = wavfile.read(file_path)
    
    # Extraction de la fréquence fondamentale et de l'amplitude moyenne
    f0 = extract_fundamental_frequency(data, sample_rate)
    mean_amplitude = np.mean(np.abs(data))
    
    return [f0, mean_amplitude]

def extract_fundamental_frequency(data, sample_rate):
    # Analyse spectrale pour trouver la fréquence fondamentale
    _, frequencies, _ = np.fft.fft(data, sample_rate)
    peaks, _ = find_peaks(np.abs(frequencies))
    
    fundamental_frequency = frequencies[peaks[0]]
    
    return fundamental_frequency

def load_data(directory):
    features = []
    labels = []
    
    for bird_directory in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, bird_directory)):
            for file_name in os.listdir(os.path.join(directory, bird_directory)):
                file_path = os.path.join(directory, bird_directory, file_name)
                features.append(extract_features(file_path))
                labels.append(bird_directory)
    
    return np.array(features), np.array(labels)

def train_classifier(features, labels, classifier_type='knn'):
    if classifier_type == 'knn':
        classifier = KNeighborsClassifier(n_neighbors=3)
    elif classifier_type == 'gmm':
        classifier = GaussianMixture(n_components=3)
    else:
        raise ValueError("Invalid classifier type. Choose 'knn' or 'gmm'.")
    
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    classifier.fit(features_scaled, labels)
    
    return classifier

# Chargement des données d'entraînement depuis le répertoire SONS
train_features, train_labels = load_data('Documents/CENTRALE 2A/ST7/TL/SONS')

# Entraînement du classifieur KNN
knn_classifier = train_classifier(train_features, train_labels, classifier_type='knn')

# Entraînement du classifieur GMM
gmm_classifier = train_classifier(train_features, train_labels, classifier_type='gmm')
