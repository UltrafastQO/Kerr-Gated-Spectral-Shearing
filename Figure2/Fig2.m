%% Processing and plotting single XPM module results
% This code can be used to process the data collected with a single XPM
% module in operation. Data collected here consist of 2D timing histograms
% between the idler photon (detector label 2) and signal photon (detector
% label 3). Histograms show correlations between click times on detector 2
% and 3, after a trigger pulse (from the master laser) is measured at
% detector 1. 

clear all; clc; close all;
cd '/Users/kfenwick/Desktop/Spectral Shearing/Github Repo/Figure2'; % required for colormap


x = 1:1:20000; % range of bins along signal axis in collected timing histograms
timestamp = x*20; % binwidth is 20ps, convert unitless x into timestamp [ps]
wavelength = -timestamp/1029+1744; % convert timestamp into wavelength

c=3e8; % speed of light in vacuum [m/s]
freq = c./(wavelength*10^(-9))*10^(-12)-c./(1550*10^(-9))*10^(-12); % frequency shift corresponding to each wavelength [Hz]

for k = 1:1:4
    power = k*5;
    for i = 1:1:71 % loop over XMP pump delays
        filename = strcat('/Users/kfenwick/Desktop/Spectral Shearing/Github Repo/Figure2/SNAKE_switch2_', num2str(power), 'mW/scan1delay', num2str(i), '.mat');
        load(filename)
        
        corrs2D = corrs123;
        y = squeeze(smoothdata(sum(corrs2D(115:160,:,:)),"movmean",20)); % sum over 1-photon contributions from idler
        spectra(:,i) = y; % add current spectrum to matrix of spectra
        
        [peakWavVal,peakWavInd] = max(y); % find peak of current spectrum
        peakWav(i,k) = wavelength(peakWavInd); % wavelength of peak 
        peakWavFreq(i,k) = freq(peakWavInd); % frequency of peak 
        
        ind = find(y>0.5*max(y)); % find range between which the spectrum is greater than half its maximum value
        bandwidthLambda(i,k) = abs(wavelength(ind(end))-wavelength(ind(1))); % use these indices to compute wavelength bandwidth (FWHM)
        bandwidthFreq(i,k) = abs(freq(ind(end))-freq(ind(1))); % use these indices to compute frequency bandwidth (FWHM)
    end
    tau = 2*powerStage*10^(-3)./c*10^(12)-0.13333;
    
    figure(k)
    set(gcf,'color','w');
    a=axes();
    box off;
    set(a,'FontSize', 18, 'LineWidth',3);
    
    %%%plot shifted spectra as a function of either wavelength or frequency
    surf(wavelength,tau,spectra.'/integrationtime,'EdgeColor','none')
    % surf(freq,tau,spectra.'/integrationtime,'EdgeColor','none')
    view(2)
    colorbar
    xlabel('Wavelength (nm)')
%     xlabel('Frequency (THz)')
    ylabel('Pump #1 Delay (ps)')
    %%%choose below xlim depending on whether you are plotting wavelength or frequency
    xlim([1450 1650])
    % xlim([(c./(1650*10^(-9))*10^(-12)-c./(1550*10^(-9))*10^(-12)) (c./(1450*10^(-9))*10^(-12)-c./(1550*10^(-9))*10^(-12))])
    ylim([-2 2])
    colormap(slanCM('viridis')); % requires being in Figure2 folder
    caxis([0 50/120])
    
    set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth', 2);
end

%%
figure(5);
f = gcf; % Get the current figure handle
f.Position(3:4) = [1200, 400];
set(gcf,'color','w');
a=axes();
box off;
set(a,'FontSize', 18, 'LineWidth',3);

subplot(1,2,1)
yyaxis left
plot(tau,peakWavFreq(:,1),'o-','Color',[0.90 0.36 0.05],'MarkerFaceColor',[0.90 0.36 0.05],'LineWidth',2,'MarkerSize',10)
hold on
plot(tau,peakWavFreq(:,2),'s-','Color',[0.37 0.43 0.76],'MarkerFaceColor',[0.37 0.43 0.76],'LineWidth',2,'MarkerSize',10)
hold on
plot(tau,peakWavFreq(:,3),'d-','Color',[0.95 0.59 0.00],'MarkerFaceColor',[0.95 0.59 0.00],'LineWidth',2,'MarkerSize',10)
hold on
plot(tau,peakWavFreq(:,4),'^-','Color',[0.65 0.25 0.69],'MarkerFaceColor',[0.65 0.25 0.69],'LineWidth',2,'MarkerSize',10)
hold on
ylim([-10 7])
ylabel('Frequency Shift (THz)')
yyaxis right
ylabel('Center Wavelength (\mu m)')
ylim([(c./(c./(1550*10^(-9))+7*10^(12)))*10^9-2 (c./(c./(1550*10^(-9))-10*10^(12)))*10^9-2])
xlim([-2 2])
xlabel('Delay 1 (ps)')

ax = gca; % Get handle to current axes (which is 'right' after the last yyaxis call)
ax.YDir = 'reverse'; 
ax.YAxis(1).Color = 'k'; % Sets the left y-axis color to black
ax.YAxis(2).Color = 'k'; % Sets the right y-axis color to black

legend('P_{XPM_1}=25 nJ','P_{XPM_1}=50 nJ','P_{XPM_1}=75 nJ','P_{XPM_1}=100 nJ','Location','SouthEast')

subplot(1,2,2)
yyaxis left
plot(tau,bandwidthFreq(:,1)-mean(bandwidthFreq(1:5,1)),'o-','Color',[0.90 0.36 0.05],'MarkerFaceColor',[0.90 0.36 0.05],'LineWidth',2,'MarkerSize',10)
hold on
plot(tau,bandwidthFreq(:,2)-mean(bandwidthFreq(1:5,2)),'s-','Color',[0.37 0.43 0.76],'MarkerFaceColor',[0.37 0.43 0.76],'LineWidth',2,'MarkerSize',10)
hold on
plot(tau,bandwidthFreq(:,3)-mean(bandwidthFreq(1:5,3)),'d-','Color',[0.95 0.59 0.00],'MarkerFaceColor',[0.95 0.59 0.00],'LineWidth',2,'MarkerSize',10)
hold on
plot(tau,bandwidthFreq(:,4)-mean(bandwidthFreq(1:5,4)),'^-','Color',[0.65 0.25 0.69],'MarkerFaceColor',[0.65 0.25 0.69],'LineWidth',2,'MarkerSize',10)
hold on
xlim([-2 2])
ylim([-1 23])
xlabel('Delay 1 (ps)')
ylabel('\Delta Bandwidth (THz)')
yyaxis right
ylabel('\Delta Bandwidth (\mu m)')
legend('P_{XPM_1}=25 nJ','P_{XPM_1}=50 nJ','P_{XPM_1}=75 nJ','P_{XPM_1}=100 nJ','Location','NorthEast')

ax = gca; % Get handle to current axes (which is 'right' after the last yyaxis call)
ax.YAxis(1).Color = 'k'; % Sets the left y-axis color to black
ax.YAxis(2).Color = 'k'; % Sets the right y-axis color to black

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth', 2);