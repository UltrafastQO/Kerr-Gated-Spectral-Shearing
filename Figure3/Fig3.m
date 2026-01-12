%% Processing and plotting two XPM module results
% This code can be used to process the data collected with two XPM
% modules in operation. Data collected here consist of 2D timing histograms
% between the idler photon (detector label 2) and signal photon (detector
% label 3). Histograms show correlations between click times on detector 2
% and 3, after a trigger pulse (from the master laser) is measured at
% detector 1. 

clear all; clc; close all;
cd '/Users/kfenwick/Desktop/Spectral Shearing/Github Repo/Figure3'; % required for colormap, USER MUST CHANGE


x = 1:1:20000; % range of bins along signal axis in collected timing histograms
timestamp = x*20; % binwidth is 20ps, convert unitless x into timestamp [ps]
wavelength = -timestamp/1029+1744; % convert timestamp into wavelength

c=3e8; % speed of light in vacuum [m/s]
freq = c./(wavelength*10^(-9))*10^(-12)-c./(1550*10^(-9))*10^(-12); % frequency shift corresponding to each wavelength [Hz]

power1 = [5 5 10 20]; % powers used in XPM module 1 [mW]; needed for data file labelling
power2 = [5 15 10 20]; % powers used in XPM module 2 [mW]; needed for data file labelling

for k = 1:1:4
    for i = 1:1:21 % loop over XMP pump delays
        %%%USER MUST CHANGE FILENAME DIRECTORY
        filename = strcat('/Users/kfenwick/Desktop/Spectral Shearing/Github Repo/Figure3/switch1_', num2str(power1(k)),'mW_switch2_', num2str(power2(k)), 'mW/scan1delay', num2str(i), '.mat');
        load(filename)
        
        delay1 = delayStage(i);
        delay2 = powerStage(i);
        
        tau1(i) = (delay1)*10^(-3)*2/3e8*10^12+10.66667-0.2; % delay in ps
        tau2(i) = -((delay2)*10^(-3)*2/3e8*10^12-0.266667)-0.2; % delay in ps
        
        corrs2D = corrs123+corrs124;
        
        y = squeeze(smoothdata(sum(corrs2D(115:160,:,:)),"movmean",20)); % sum over 1-photon contributions from idler
        corrSpec(:,:,i) = y; % add current spectrum to matrix of spectra
        
        [peakWavVal,peakWavInd] = max(y); % find peak of current spectrum
        peakWav(i,:) = freq(peakWavInd); % frequency of peak
        
        for j = 1:1:21
            ind = find(y(:,j)>0.5*max(y(:,j))); % find range between which the spectrum is greater than half its maximum value
            bandwidth(j,i) = abs(freq(ind(end))-freq(ind(1))); % use these indices to compute frequency bandwidth (FWHM)
        end
    end
    
    t1 = 3:1:21; % tau1 time indices needed for snake plot (i.e., when tau1=tau2, gives the same phase shift)
    t2 = 20:-1:2; % tau2 time indices needed for snake plot (i.e., when tau1=tau2, gives the same phase shift)

    for i = 1:1:length(t1)
        spectra(:,i)=corrSpec(:,t1(i),t2(i)); % store spectra when tau1=tau2
    end
    
    figure(3*k-2)
    set(gcf,'color','w');
    a=axes();
    box off;
    set(a,'FontSize', 18, 'LineWidth',3);
    
    %%%plot shifted spectra as a function of either wavelength or frequency
    surf(wavelength,-tau1(3:1:21),spectra.'/integrationtime,'EdgeColor','none')
%     surf(freq,-tau1(3:1:21),spectra.'/integrationtime,'EdgeColor','none')
    view(2)
    colorbar
    xlabel('Wavelength (nm)')
%     xlabel('Frequency (THz)')
    ylabel('Pump #1 Delay (ps)')
    ylim([-1.8 1.8])
    %%%choose below xlim depending on whether you are plotting wavelength or frequency
    xlim([1400 1700])
%     xlim([(c./(1700*10^(-9))*10^(-12)-c./(1550*10^(-9))*10^(-12)) (c./(1400*10^(-9))*10^(-12)-c./(1550*10^(-9))*10^(-12))])
    colormap(slanCM('viridis'));
    caxis([0 1])
    
    set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth', 2);
    
    figure(3*k-1)
    set(gcf,'color','w');
    a=axes();
    box off;
    set(a,'FontSize', 18, 'LineWidth',3);
    
    surf(tau1,tau2,peakWav.','EdgeColor','none')
    view(2)
    colorbar
    xlabel('\Delta\tau_1 (ps)')
    ylabel('\Delta\tau_2 (ps)')
    xlim([-1.8 1.8])
    ylim([-1.8 1.8])
    colormap(slanCM('fusion'));
    caxis([-10 10])
    
    set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth', 2);
    
    figure(3*k)
    set(gcf,'color','w');
    a=axes();
    box off;
    set(a,'FontSize', 18, 'LineWidth',3);
    
    surf(tau1,tau2,(bandwidth-mean(bandwidth(1:3,1))).','EdgeColor','none')
    view(2)
    colorbar
    xlabel('\Delta\tau_1 (ps)')
    ylabel('\Delta\tau_2 (ps)')
    xlim([-1.8 1.8])
    ylim([-1.8 1.8])
    colormap(slanCM('viola'));
    
    caxis([-10 10])
    
    set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth', 2);
end