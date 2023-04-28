close all; clc; clear all;
format long e;
% scenario parameters
scenario.carrier_freq = 2e9; % carrier frequency in Hz
scenario.bandwidth = 250e3; % bandwidth in Hz
scenario.path_loss_exponent = 2.5; %
scenario.noise_spectral_density = -174; % dBm
scenario.harvesting_threshold = -35; % dBm %Da vedere
scenario.phase_off = 1; %Offset, ambiguità di fase
scenario.packet_length = 1;
scenario.shadowing = 1;
scenario.HOM = 2; %Higher Modulation Order
scenario.niter = 1e5; %Numero iterazioni
scenario.maxerrs = 1e4; %Numero massimo di errori
sigma_lognormaldB = [3];
sigma_lognormalNeper = sigma_lognormaldB*log(10)/10;

% LINK1
l1.sigma_nakagamif = 1;
l1.sigma_nakagamib = 1; %stessa nakagami per RX1-tag1 e RX1-tag2?
l1.mf_nakagami = 4;
l1.mb_nakagami = 1;
% l1.sigma_lognormalf = sigma_lognormalNeper;
% l1.sigma_lognormalb = sigma_lognormalNeper;

scenario.link1 = l1;

% LINK1
l2.sigma_nakagamif = 1;
l2.sigma_nakagamib = 1; %stessa nakagami per RX2-tag1 e RX2-tag2?
l2.mf_nakagami = 4;
l2.mb_nakagami = 1;
% l2.sigma_lognormalf = sigma_lognormalNeper;
% l2.sigma_lognormalb = sigma_lognormalNeper;

scenario.link2 = l2;

alpha = 10;% alpha >= 1, path gain
scenario.alpha=alpha;

Ptx0 = 0.2;
scenario.ptx0=Ptx0;

N0 = db2pow(scenario.noise_spectral_density - 30); %Per averla in Watt/Hz
Pn = N0*scenario.bandwidth; %Potenza del rumore in Watt
scenario.pn=Pn;

%assegno le coordinate cartesiane ai due tag e le metto nella struct
% coordinate_tag1=[0 4; 0 3 ; 0 3];
% coordinate_tag2=[0 -4; 0 -3 ; 0 -3];
coordinate_tag1=[0 2; 0 2];
coordinate_tag2=[0 -2; 0 -2];
parametri_scenario.coordinate_tag1=coordinate_tag1;
parametri_scenario.coordinate_tag2=coordinate_tag2;
%assegno le coordinate polari ai due tag e le metto nella struct
% coordinate_rx1=[3 pi; 3 pi; 3 pi];
% coordinate_rx2=[3 0; 3 0; 3 0];
coordinate_rx1=[39.95 pi; 39.95 pi];
coordinate_rx2=[39.95 0; 39.95 0];
parametri_scenario.coordinate_rx1=coordinate_rx1;
parametri_scenario.coordinate_rx2=coordinate_rx2;

%Chiamo la function scenario che restituisce gli angoli Theta12_i, per le
%matrice di correlazioni
setup = setup_scenario(parametri_scenario);
theta_12_1=setup.Theta_12_1;
theta_12_2=setup.Theta_12_2;

%Chiamo la function calcolo_snr che restituisce il snr medio in una struct
av_rx_snr = calcolo_snr(setup, scenario);
snr_1=av_rx_snr.snr1;
snr_2=av_rx_snr.snr2;


gamma = 1; %é il coefficiente di riflessione di uno dei due tag, normalizzato ad 1 (in questo caso uguali)
av_rx_sir = pow2db(1./gamma); %Il primo, ha Reflection coefficient=1, Ricorda SIR>=1


estim_err_var = 0; %Ho una perfetta stima del canale, canale flat, stimo una volta e addio
HOM = [2]; %Higher Modulation Order
l_var = length(snr_1); %Lo uso per ciclare sul numero di scenari da simulare
l_param = length(HOM);

BER1_SIC_MULTI = zeros(l_param,l_var);
BER2_SIC_MULTI = zeros(l_param,l_var);
CONDBER1_SIC = zeros(l_param,l_var);
SNRCOND1_SIC = zeros(l_param,l_var);
CONDBER2_SIC = zeros(l_param,l_var);
SNRCOND2 = zeros(l_param,l_var);
POW_RATE1 = zeros(l_param,l_var);
POW_RATE2 = zeros(l_param,l_var);
POW_DIFF = zeros(l_param,l_var);
CUM_BIT_SIC_MULTI = zeros(l_param,l_var);
BER1_SIC = zeros(l_param,l_var);
BER2_SIC = zeros(l_param,l_var);
CUM_BIT_SIC = zeros(l_param,l_var);
sigma = sigma_lognormalNeper;

A=5.0130e-05;
B=0.99994987;


for i=1:l_param
    
    scenario.HOM = HOM(i);
    for j=1:l_var
    %Sigmab è la matrice di correlazione    
        Sigmab_tag1 = sigma^2*[(A+B) (A*cos(theta_12_1(i))+B);...
        (A*cos(theta_12_1(i))+B) (A+B)]; 
        %l1.sigma_lognormalf = sigma;
        l1.Sigmab = Sigmab_tag1;

        Sigmab_tag2 = sigma^2*[(A+B) (A*cos(theta_12_2(i))+B);...
        (A*cos(theta_12_2(i))+B) (A+B)]; 
        %l2.sigma_lognormalf = sigma;
        l2.Sigmab = Sigmab_tag2;
        scenario.link1 = l1;
        scenario.link2 = l2;
       
        %Prendo i snr per lo scenario i-esimo da simulare
        av_rx_snr1=snr_1(i,:);
        av_rx_snr2=snr_2(i,:);

        [ber1_div,ber2_div,snr1,snr2,ber1,ber2,powering_rate1,powering_rate2,cum_bit_div,cum_bit] = multi_Shadowing_SIC(scenario,scenario.alpha,av_rx_snr1,av_rx_snr2,av_rx_sir,estim_err_var);
        BER1_SIC_MULTI(i,j) = ber1_div;
        BER2_SIC_MULTI(i,j) = ber2_div;
        BER1_SIC(i,j) = ber1;
        BER2_SIC(i,j) = ber2;
        POW_RATE1(i,j) = powering_rate1;
        POW_RATE2(i,j) = powering_rate2;
        SNRCOND1(i,j) = snr1;
        SNRCOND2(i,j) = snr2;
        CUM_BIT_SIC_MULTI(i,j) = cum_bit_div;
        CUM_BIT_SIC(i,j) = cum_bit;
    end
end