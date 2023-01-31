clear all
close all
clc


%% Suppression du bruit provoqué par les mouvements 
figure;
%Question 1:


load("ecg.mat"); 


x=ecg;


%Question 2:

fs=500;

N = length(x); 
ts=1/fs;
t=(0:N-1)*ts;
subplot(1,3,1)
plot(t,x);
title("le signal ECG ");


%Question 3:

y = fft(x);
f =(0:N-1)*(fs/N);

f2 = (-N/2:N/2-1)*(fs/N);
subplot(1,3,2)
plot(f2,fftshift(abs(y))) 
title("spectre Amplitude du signal ECG")

filtre_haut = ones(size(x));

fc = 0.5;
index_h = ceil(fc*N/fs); 
filtre_haut(1:index_h)=0; 
filtre_haut(N-index_h+1:N)=0; 

filtre=filtre_haut.*y;
ecg1=ifft(filtre,"symmetric");

%Question 4:
subplot(1,3,3)
plot(t,ecg1);
title("signal filtre ecg1")


 %% Suppression des interférences des lignes électriques 50Hz
% Filtre Notch

subplot(2,3,4)
plot(f2,fftshift(abs(fft(ecg1))))
title("TF de ecg1")




%Question 5:
figure;
notch_filter = ones(size(x)); 

fc2 = 50;

index_h2 = ceil((fc2*N)/fs)+1; 

notch_filter(index_h2)= 0;
notch_filter(N-index_h2+2)= 0; 

filtre2=notch_filter.*filtre;

%Question 6:

ecg2=ifft(filtre2,"symmetric");
subplot(1,3,1)
plot(f2,fftshift(abs(filtre2)));
title("Spectre de ECG2")
subplot(1,3,2)
plot(t,ecg2);
title("Signal  ECG2")

% subplot(211)
% plot(t,x)
% xlim([0.5 1.5])
% subplot(212)
% plot(t,ecg2)
% xlim([0.5 1.5])


%% Amélioration du rapport signal sur bruit



%Question 7:

figure;
filtre_bas = zeros(size(x));
fc3 = 37;
index_h3 = ceil(fc3*N/fs);
filtre_bas(1:index_h3)=1;
filtre_bas(N-index_h3+1:N)=1; 
ecg3_freq =  filtre_bas.*fft(ecg2);

ecg3 = ifft(ecg3_freq,"symmetric");
 %Question 8:
 
 
subplot(1,3,1)
plot(t,x)
xlim([0.5 1.5])
title('signal de depart ecg')
subplot(1,3,2)
plot(t,ecg3)
xlim([0.5 1.5])
title('signal ecg3')

%% Identification de la fréquence cardiaque avec la fonction dautocorrélation

%Question 9:

figure;
[acf,lags] = xcorr(ecg3,ecg3); 
stem(lags/fs,acf)
