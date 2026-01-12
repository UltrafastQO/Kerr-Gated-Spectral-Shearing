%% Processing and plotting two XPM module results - fine control over spectral bandwidth
% This code can be used to process the data collected with two XPM
% modules in operation. Data collected here consist of 2D timing histograms
% between the idler photon (detector label 2) and signal photon (detector
% label 3). Histograms show correlations between click times on detector 2
% and 3, after a trigger pulse (from the master laser) is measured at
% detector 1. 

clear all; clc; close all;
cd '/Users/kfenwick/Desktop/Spectral Shearing/Github Repo/Figure4'; % required for colormap, USER MUST CHANGE

x = 1:1:20000; % range of bins along signal axis in collected timing histograms
timestamp = x*20; % binwidth is 20ps, convert unitless x into timestamp [ps]
wavelength = -timestamp/1029+1744; % convert timestamp into wavelength

c=3e8; % speed of light in vacuum [m/s]
freq = c./(wavelength*10^(-9))*10^(-12)-c./(1550*10^(-9))*10^(-12); % frequency shift corresponding to each wavelength [Hz]

gaussEqn = 'a*exp(-4*log(2)*(x-b)^2/c^2)+d';
% excludePoints = freq > 8; % USER MUST CHANGE CUTOFF VALUE, for datasets where pump noise is more present, exclusion for high frequencies is necessary to obtain a good fit

startPoints = [10 0 2 0.1]; % USER MUST CHANGE VALUES, [a b c d] in 'gaussEqn'

targetBW = ["10"; "20"; "30"; "40"; "50"; "60"; "70"; "80"; "90"; "100"];
power1 = ["10"; "0"; "2_7"; "5_7"; "5"; "6"; "8"; "11"; "15"; "20"];
power2 = ["10"; "0"; "0"; "0"; "5"; "6"; "8"; "11"; "15"; "20"];

for k = 1:1:10
    for i = 1:1:20
        filename = strcat('/Users/kfenwick/Desktop/Spectral Shearing/Github Repo/Figure4/Bandwidth Manipulation/switch1_', power1(k), 'mW_switch2_', power2(k),'mW_BW_', targetBW(k), 'nm/scan1delay', num2str(i), '.mat');
        load(filename)
        
        corrs2D = corrs123+corrs124;
        y = squeeze(smoothdata(sum(corrs2D(115:160,:,:)),"movmean",20)); % sum over 1-photon contributions from idler
        
        spectra(:,i) = y;  % add current spectrum to matrix of spectra
        
        [peakWavVal,peakWavInd] = max(y); % find peak of current spectrum
        peakWav(i) = wavelength(peakWavInd); % wavelength of peak 
        peakWavFreq(i) = freq(peakWavInd); % frequency of peak 
        
        ind = find(y>0.5*max(y)); % find range between which the spectrum is greater than half its maximum value
        bandwidthLambda(i) = abs(wavelength(ind(end))-wavelength(ind(1))); % use these indices to compute wavelength bandwidth (FWHM)
        bandwidthFreq(i) = abs(freq(ind(end))-freq(ind(1))); % use these indices to compute frequency bandwidth (FWHM)
        
%         f1 = fit(freq.',spectra(:,i),gaussEqn,'Start', startPoints,'Exclude', excludePoints); % use this if wanting to exclude some points in the fit (for noisier data)
        f1 = fit(freq.',spectra(:,i),gaussEqn,'Start', startPoints);
%         plot(f1,freq,spectra(:,i)) % If visualizing the fit is desired, recommended if tuning the fit parameters (a, b, c, d)!
%         drawnow
        
        fitParams = confint(f1,0.95);
        
        fittedShift(i) = (fitParams(1,2)+fitParams(2,2))/2;
        fittedShiftErr(i) = fitParams(2,2)-fittedShift(i);
        
        fittedBW(i) = (fitParams(1,3)+fitParams(2,3))/2;
        fittedBWErr(i) = fitParams(2,3)-fittedBW(i);
        
        fittedRelShift(i) = fittedShift(i)/2.54;

    end
    
    relFreqShift = peakWavFreq./bandwidthFreq;
    
    avgSpec(:,k) = mean(spectra,2);
    errSpec(:,k) = std(spectra,0,2);
    avgBWwav(k) = mean(bandwidthLambda);
    stdBWwav(k) = std(bandwidthLambda);
    avgBWfreq(k) = mean(bandwidthFreq);
    stdBWfreq(k) = std(bandwidthFreq);
    avgPeakWavFreq(k) = mean(peakWavFreq);
    errPeakWavFreq(k) = std(peakWavFreq);
    
    avgFittedShift(k) = mean(fittedShift);
    stdFittedShift(k) = std(fittedShift);
    avgFittedBW(k) = mean(fittedBW);
    stdFittedBW(k) = std(fittedBW);
    avgfittedRelShift(k) = mean(fittedRelShift);
    stdfittedRelShift(k) = sqrt(std(fittedRelShift)^2+0.02^2);
    
    avgRelFreqShift(k) = mean(relFreqShift);
    stdRelFreqShift(k) = std(relFreqShift);
    
end

%% bandwidth expanded spectra

figure('Position',[600 100 600 300]);
set(gcf,'color','w');
a=axes();
set(a,'FontSize', 18, 'LineWidth',3);

target = -50:10:50;

colours = ['#4B0082',"#000000", "#8B008B","#9370DB","#9932CC","#BA55D3","#DA70D6","#EE82EE","#DDA0DD","#D8BFD8"];

ii=0;
for i = 1:1:10
    plot(freq,10^3*(avgSpec(:,i)./sum(avgSpec(:,i))).','-','Color',colours(i),'LineWidth',3)
    hold on
    ii=ii-1;
end

xlabel('Freq shift (THz)')
ylabel('Target BW (nm)')
xlim([-15 15])
legend('C1','Input','E1','E2','E3','E4','E5','E6','E7','E8','location','east outside')

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth', 2);
