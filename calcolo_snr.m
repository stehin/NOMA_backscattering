%La function restituisce i 4*N av_rx_snr come struct, composta da 4
%vettori: 
% av_rx_snr1 (componente 1: snr11 SNR ricevuto da RX1 quando backscattera TAG1,
% componente 2: snr12 SNR ricevuto da RX1 quando backscattera TAG2)
% av_rx_snr2 (componente 1: snr21 SNR ricevuto da RX2 quando backscattera TAG1,
% componente 2: snr22 SNR ricevuto da RX2 quando backscattera TAG2),
function [av_rx_snr] = calcolo_snr(setup, scenario)
    av_rx_snr1=zeros(length(setup.df1),2);
    av_rx_snr2=zeros(length(setup.df1),2);
    for i=1:length(setup.df1)
        df=[setup.df1(i) setup.df2(i)]; %Distanze forward,Hf1 e poi Hf2
        
        d11 = setup.d11(i); % dist(Tag1,RX1) all'iterazione i
        d12 = setup.d12(i); % dist(Tag2,RX1) all'iterazione i
        d21 = setup.d21(i); % dist(Tag1,RX2) all'iterazione i
        d22 = setup.d22(i); % dist(Tag2,RX2) all'iterazione i

        drx1 = [d11 d12]; % [dist(Tag1,RX1) dist(Tag2,RX1)]
        drx2 = [d21 d22]; % [dist(Tag1,RX2) dist(Tag2,RX2)]
        
        wvlngt = freq2wavelen(scenario.carrier_freq);
        Plf = (wvlngt./(4*pi*df)).^scenario.path_loss_exponent; %Calcola il path loss forward

        Plb1 = (wvlngt./(4*pi*drx1)).^scenario.path_loss_exponent; %Calcola il path loss backward dai tag a RX1
        Plb2 = (wvlngt./(4*pi*drx2)).^scenario.path_loss_exponent; %Calcola il path loss backward dai tag a RX2
        Pl1 = Plf.*Plb1; % Path losses at RX1
        Pl2 = Plf.*Plb2; % Path losses at RX2

        av_rx_snr_1 = pow2db(scenario.alpha*scenario.ptx0*Pl1/scenario.pn); %SNRs at RX1 in dB
        av_rx_snr_2 = pow2db(scenario.alpha*scenario.ptx0*Pl2/scenario.pn); %SNRs at RX2 in dB

        av_rx_snr1(i,:)=av_rx_snr_1;
        av_rx_snr2(i,:)=av_rx_snr_2;
    end
    av_rx_snr.snr1=av_rx_snr1;
    av_rx_snr.snr2=av_rx_snr2;
