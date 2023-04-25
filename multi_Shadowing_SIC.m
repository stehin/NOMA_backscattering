function [ber1_div,ber2_div,SNR1,SNR2,ber1,ber2,powering_rate1,powering_rate2,cumulative_bits_div,cumulative_bits] = multi_Shadowing_SIC(scenario,alpha,av_rx_snr1,av_rx_snr2,av_rx_sir,estim_err_var)

    % System parameters
    
    %fc = scenario.carrier_freq; % carrier frequency [Hz]
    
    %wave_length = freq2wavelen(fc); % wavelength [m] MOVE THIS
    
    %k0 = scenario.path_loss_exponent; % path-loss exponent
    
    W = scenario.bandwidth; % available bandwidth [Hz]
    
    NodBm = scenario.noise_spectral_density;
    No = db2pow(NodBm - 30); % noise power spectral density [W/Hz]
    
    Pn = No*W; % noise power
    
    %gamma_pdBm = scenario.harvesting_threshold;
    %gamma_p = gamma_pdBm -30; % tags' sensitivity threshold [dBW]
    
    n_iter = scenario.niter;
    max_num_err = scenario.maxerrs;
    monostatic = scenario.monostatic;
    
    l_pkt = scenario.packet_length;
    %GAMMA = [gamma_min gamma_max];
    
    phase_flag = scenario.phase_off;
    
    if monostatic
        
    %     %db1 = df1;
    %     PLfs_f1 = (wave_length./(4*pi*df1)).^k0; % free space attenuation vector
    %     %PLfs_b1 = (wave_length./(4*pi*db1)).^k0; % free space attenuation vector
    %     
    %     PLfs1 = PLfs_f1.^2;
    %     
    %     % Link 2 pareameters
    %     
    %     %db2 = df2;
    %     PLfs_f2 = (wave_length./(4*pi*df2)).^k0; % free space attenuation vector
    %     %PLfs_b2 = (wave_length./(4*pi*db2)).^k0; % free space attenuation vector
    %     
    %     PLfs2 = PLfs_f2.^2;
        
        sigma_nakagami1f = scenario.link1.sigma_nakagamif;
        
        sigma_nakagami2f = scenario.link2.sigma_nakagamif;
        
        mf_nakagami1 = scenario.link1.mf_nakagami;
        
        mf_nakagami2 = scenario.link2.mf_nakagami;
        
        pdf1 = makedist('Nakagami','mu',mf_nakagami1,'omega',sigma_nakagami1f);
        pdf2 = makedist('Nakagami','mu',mf_nakagami2,'omega',sigma_nakagami2f);
        
        
    else
        
    %     % Link 1 pareameters
    %     
    %     db1 = dtx_rx - df1;
    %     PLfs_f1 = (wave_length./(4*pi*df1)).^k0; % free space attenuation vector
    %     PLfs_b1 = (wave_length./(4*pi*db1)).^k0; % free space attenuation vector
    %     
    %     PLfs1 = PLfs_f1*PLfs_b1;
    %     
    %     
    %     % Link 2 pareameters
    %     
    %     db2 = dtx_rx - df2;
    %     PLfs_f2 = (wave_length./(4*pi*df2)).^k0; % free space attenuation vector
    %     PLfs_b2 = (wave_length./(4*pi*db2)).^k0; % free space attenuation vector
    %     
    %     PLfs2 = PLfs_f2*PLfs_b2;
        
        sigma_nakagami1f = scenario.link1.sigma_nakagamif;
        sigma_nakagami1b = scenario.link1.sigma_nakagamib;
        
        
        sigma_nakagami2f = scenario.link2.sigma_nakagamif;
        sigma_nakagami2b = scenario.link2.sigma_nakagamib;
        
        
        mf_nakagami1 = scenario.link1.mf_nakagami;
        mb_nakagami1 = scenario.link1.mb_nakagami;
        
        mf_nakagami2 = scenario.link2.mf_nakagami;
        mb_nakagami2 = scenario.link2.mb_nakagami;
        
        pdf1 = makedist('Nakagami','mu',mf_nakagami1,'omega',sigma_nakagami1f); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
        pdb1 = makedist('Nakagami','mu',mb_nakagami1,'omega',sigma_nakagami1b); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
        pdf2 = makedist('Nakagami','mu',mf_nakagami2,'omega',sigma_nakagami2f); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
        pdb2 = makedist('Nakagami','mu',mb_nakagami2,'omega',sigma_nakagami2b); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
        if(scenario.shadowing)
    %     sigma_lognormal1f = scenario.link1.sigma_lognormalf;
        Sigma_b1 = scenario.link1.Sigmab;
    %     sigma_lognormal2f = scenario.link2.sigma_lognormalf;
        Sigma_b2 = scenario.link2.Sigmab;
    %     lnf1 = makedist('Normal','mu',0,'sigma',sigma_lognormal1f);
        %lnb1 = makedist('Normal','mu',0,'sigma',sigma_lognormal1b);
    %     lnf2 = makedist('Normal','mu',0,'sigma',sigma_lognormal2f);
        %lnb2 = makedist('Normal','mu',0,'sigma',sigma_lognormal2b);
        end
        
    end
    
    %forward_pl_vector = [PLfs_f1 PLfs_f2]; % path-loss for the forward links
    %PL_vector = [PLfs1 PLfs2]; % composite path-loss
    
    errors_div = zeros(1,2);
    errors = zeros(1,2);
    backscattered_bits = zeros(1,2);
    conditioned_errors = zeros(1,2);
    succ_bit_div = zeros(1,2);
    succ_bit = zeros(1,2);
    bit_both_powered = 0;
    bit_not_powered = zeros(1,2);
    bit_count = 0;
    
    %cond_snr1 = [];
    %cond_snr2 = [];
    if(scenario.HOM > 2)
        MQAM = scenario.HOM;
        d = sqrt(6/(MQAM-1));
        n0 =  (sqrt(MQAM)/2-0.5);
        count=1;
        %Genero ala costellazione di una MQAM e la salvo come vettore riga
        for nn=0:(sqrt(MQAM)-1)
            for mm = 0:(sqrt(MQAM)-1)
                sn(count,1) = (nn - n0)*d + 1i*(mm - n0 )*d;
                count = count + 1;
            end
        end
        sn = sn';
    else
        %Genero ala costellazione binaria
        sn = [-1 1];
    end
    for i=1:n_iter
        
        if monostatic
            
            % generate channel gain for the forward link
            af1 = random(pdf1,1,1);
    
            phif1 = phase_flag*2*pi*rand(1,1);
            
            hf1 = af1.*exp(1i*phif1);
    
           
            
            % generate channel gain for the reverse link
            hb1 = hf1;
            
            h1 = hf1*hb1;
            
            % generate channel gain for the forward link
            af2 = random(pdf2,1,1);
            phif2 = phase_flag*2*pi*rand(1,1);
            
            hf2 = af2.*exp(1i*phif2);
            
            % generate channel gain for the reverse link
            hb2 = hf2; % hb2 = 1 for test
            
            h2 = hf2*hb2;
            
        else
            
    % %         % generate channel gain for the forward link
    % %         %af1 = gamrnd(mf_nakagami1, sigma_nakagami1f/mf_nakagami1,1,1);
             af1 = random(pdf1,1,1); %Ho deciso che segue la Nakagami
             phif1 = phase_flag*2*pi*rand(1,1); %Fase uniformemente distribuita slide 59
    % %         
    % %         %hf1 = sqrt(af1).*exp(1i*phif1);
             hf1 = af1.*exp(1i*phif1); %slide 59, def dell'inviluppo del segnale ricevuto
    % % 
    % %           % generate channel gain for the forward link
    % %         %af2 = gamrnd(mf_nakagami2, sigma_nakagami2f/mf_nakagami2,1,1);
             af2 = random(pdf2,1,1);
    % %         
             phif2 = phase_flag*2*pi*rand(1,1);
    % %         
    % %         %hf2 = sqrt(af2).*exp(1i*phif2);
              hf2 = af2.*exp(1i*phif2);
    %           hf2 = hf1;
    % %         
    % %         
    % %         % generate channel gain for the reverse link (dal tag al
                %ricevitore. Ho due tag e due ricevitori, 4 backward links)
    % %         %ab1 = gamrnd(mb_nakagami1, sigma_nakagami1b/mb_nakagami1,1,1);
             ab11 = random(pdb1,1,1);
             ab12 = random(pdb1,1,1);
    % %        
             phib11 = phase_flag*2*pi*rand(1,1);
             phib12 = phase_flag*2*pi*rand(1,1);
    % %         
    % %         %hb1 = sqrt(ab1).*exp(1i*phib1);
             hb11 = ab11.*exp(1i*phib11);
             hb12 = ab12.*exp(1i*phib12);
    %           hb11 = hf1;
    %           hb12 = hf1;
    % %         
             h11 = hf1*hb11;
             h12 = hf2*hb12;
    % % 
    % %          h1 = hb1;
    % %         
    % %       
    % %         
    % %         % generate channel gain for the reverse link
    % %         %ab2 = gamrnd(mb_nakagami2, sigma_nakagami2b/mb_nakagami2,1,1);
            ab21 = random(pdb2,1,1);
            ab22 = random(pdb2,1,1);
    % %        
             phib21 = phase_flag*2*pi*rand(1,1);
             phib22 = phase_flag*2*pi*rand(1,1);
    % %         
    % %         %hb2 = sqrt(ab2).*exp(1i*phib2);
             hb21 = ab21.*exp(1i*phib21);
             hb22 = ab22.*exp(1i*phib22);
    %           hb21 = hb11;
    %           hb22 = hb12;
    % %         
             h21 = hf1*hb21;
             h22 = hf2*hb22;
    % %     
    % %         
    % %         %h2 = hb2;
    % % 
            %Aggiunta Shadowing, solo nel backward
            if(scenario.shadowing)
    %              X_b1 = mvnrnd([0 0],Sigma_b1);
                 X_b1 = mvnrnd([0 0],Sigma_b1); %genera una v.a. gaussiana bidimensionale
    %              X_b2 = mvnrnd([0 0],Sigma_b2);
                X_b2 = mvnrnd([0 0],Sigma_b2);
    %              x_s1f = random(lnf1,1,1);
    %              x_s2f = random(lnf2,1,1);
    %              x_s1f = X_s(1);
    %              x_s2f = X_s(2);
                 x_s11b = X_b1(1);
                 x_s11 = x_s11b;
    %             x_s11 = x_s11b; 
                x_s12b = X_b2(1);
                x_s12 =  x_s12b;
    %             x_s12 = x_s12b;
                x_s21b = X_b1(2);
                 x_s21 =  x_s21b;
    %            x_s21 = x_s21b;
                x_s22b = X_b2(2);
                 x_s22 = x_s22b;
    %             x_s22 = x_s22b;
                Shadowing1 = [x_s11 x_s12];
                Shadowing1 = exp(Shadowing1);
                Shadowing2 = [x_s21 x_s22];
                Shadowing2 = exp(Shadowing2);
    %             Shadowing2 = Shadowing1;
        end
        end
        % generate modulated symbols (BPSK)
        
    %     x = randi([0 1],l_pkt,2); % soruce bits
    %     
    %     xmod = 2*x - 1; % modulated symbols
        xmod = [datasample(sn,1) datasample(sn,1)]; %prendo due simboli a caso della costellazione
    %     xmod = qammod(x,M,'UnitAveragePower',true);
        forward_ch_vector = [hf1 hf2]; % channel gains of the forward links
        %forward_ch_vector = [1 1];
    %     h11 = 3;
    %         h12 = 3;
    %         h21 = 3;
    %         h22 = 3;
      
        if(scenario.shadowing)
            CH_vector1 = [h11 h12].*Shadowing1; % channel gains of the composite links
            %Lo shadowing agisce sia lungo Re che lungo Im
            CH_vector2 = [h21 h22].*Shadowing2;
        else
            
            CH_vector1 = [h11 h12]; % channel gains of the composite links
            CH_vector2 = [h21 h22];
        end
        
        %Prx_sensors = Ptx*[forward_pl_vector(1)*abs(forward_ch_vector(1)).^2 forward_pl_vector(2)*abs(forward_ch_vector(2)).^2];
        %Prx_sensors = Ptx*[abs(forward_ch_vector(1)).^2 abs(forward_ch_vector(2)).^2];
        if(scenario.powering) % powering condition 
        powered = alpha.*(abs(forward_ch_vector).^2) >= 1;
        
    %     disp( Prx_sensors - db2pow(gamma_p)./forward_pl_vector)
    %     
    %     disp(abs(forward_ch_vector).^2)
    %     pause
      
        powered_index = find(powered); %Per sapere la posizione di quelli che hanno restituito
        %un esito True dopo il check.
        else
            powered_index = [1 2];
        end
        
    %                         disp(powered)
    %                         pause
        
        if isempty(powered_index) %Nessun bit è trasmesso, non si attiva il diodo
            
            %disp('if 1')
            
            
            % detection with noise only
            %errors(1) = errors(1) + randi([0,1],1,1);
            %errors(2) = errors(2) + randi([0,1],1,1);
            
            %Incremento il numero di bit non "trasmessi"
            bit_not_powered(1) = bit_not_powered(1) + 1;
            bit_not_powered(2) = bit_not_powered(2) + 1;
            
            
        elseif length(powered_index) == 1 %Un solo tag si è attivato
            
            %disp('if 2')
            %errors(not(powered)) = errors(not(powered)) + randi([0,1],1,1);
            
            %incremento il numero di bit backscatterati per il tag che si è
            %attivato
            backscattered_bits(powered_index) = backscattered_bits(powered_index) + 1;
            %Incrementa il numero di bit non backscatterati per il tag che non 
            %si è attivato, utilizzando l'indicizzazione logica
            bit_not_powered(not(powered)) = bit_not_powered(not(powered)) + 1;
            
            %Considero solo il canale utilizzato dal tag che si è attivato
            CH_vector1_ = CH_vector(powered_index);
            %Considero solo il simbolo che il tag che si è attivato
            %backscattera
            xmod_ = xmod(:,powered_index);
            
            %PL_vector_ = PL_vector(powered_index);
            
            % AWGN
            w1 = sqrt(Pn)*(randn(l_pkt,1)+1i*randn(l_pkt,1))/sqrt(2); %rand genera una 
            %matrice l_pkt x 1 di v.a. con distribuzione gaussiana
            %standardizzata. Lo genero per entrambe le direzioni
            w2 = sqrt(Pn)*(randn(l_pkt,1)+1i*randn(l_pkt,1))/sqrt(2);
            
            % etimation error (nel nostro caso stima perfetta del canale,
            % estim_err_var=0
            err = sqrt(estim_err_var)*(randn(1,1)+1i*randn(1,1))/sqrt(2);
            
            %Trovo la potenza ricevuta, senza il contributo del rumore
            rx_power1 = db2pow(av_rx_snr1).*Pn;
            rx_power2 = db2pow(av_rx_snr1).*Pn;
            
            % received signal
            y1 = sqrt(rx_power1).*CH_vector1_.*xmod_ + w1;
            y2 = sqrt(rx_power2).*CH_vector2_.*xmod_ + w2;
            
            v = [-1 1];
            v = v(ones(1,l_pkt),:);
            
            %Eseguo il SIC su entrambi i ricevitori, per un unico simbolo,
            %quello backscatterato
            tilde_x1 = y - sqrt(rx_power1)*(CH_vector1_ + err)*v;
            tilde_x1 = abs(tilde_x1).^2;
            
            [~,ix] = min(tilde_x1,[],2);
            
            xdec1 = ix - 1;
    
            tilde_x2 = y - sqrt(rx_power2)*(CH_vector2_ + err)*v;
            tilde_x2 = abs(tilde_x2).^2;
            
            [~,ix] = min(tilde_x2,[],2);
    
            xdec2 = ix - 1;
            
            %xdec = 1;
            
            %estim_channel = CH_vector + err;
            
            %         if real(y./estim_channel) < 0
            %
            %             xdec = -1;
            %
            %         end
            
            %Aggiorna le statistiche sui bit decodificati correttamente, o sui
            %bit sbagliati
            if not(isequal(xdec1,x(:,powered_index))|isequal(xdec2,x(:,powered_index)))%xdec ~= x(powered_index)
                
                errors_div(powered_index) =  errors_div(powered_index) + 1;
               
            else
                
                succ_bit_div(powered_index) = succ_bit_div(powered_index) + 1;
                
            end
            
        else
            
            %disp('if 3')
            bit_both_powered = bit_both_powered + 1; %Aggiorno il numero di casi in cui
            %entrambi i bit sono backscatterati
            
            backscattered_bits = backscattered_bits + 1;
            
            % etimation error
            %Nel nostro caso perfetta stima di canale
            err = sqrt(estim_err_var)*(randn(1,2)+1i*randn(1,2))/sqrt(2);
            
    %         rx_power1 = db2pow(av_rx_snr1)*Pn;
    %         rx_power2 = db2pow(av_rx_snr2)*Pn;
            
            %POW_vector = Ptx*[PL_vector(1)*abs(CH_vector(1)+err(1)).^2 PL_vector(2)*abs(CH_vector(2)+err(2)).^2]; % backscattered powered vector
            
    %         pow_vector1 = [rx_power1 rx_power1];
    %         pow_vector2 = [rx_power2 rx_power2];
    
    %         pow_vector1 = rx_power1; 
    %         pow_vector2 = rx_power2;
               
            pow_vector1 = av_rx_snr1;
            pow_vector2 = av_rx_snr2;
    
    
            
    %         Ps1 = Ptx*GAMMA(1)*PL_vector(1)*abs(CH_vector(1)).^2;
    %         cond_snr1 = cat(2,cond_snr1,Ps1/Pn);
            
    %         Ps2 = Ptx*GAMMA(2)*PL_vector(2)*abs(CH_vector(2)).^2;
    %         cond_snr2 = cat(2,cond_snr2,Ps2/Pn);
            
            % Ordering the channel gains in order to assign the reflection
            % coefficients
            [~,index1] = sort(abs(CH_vector1).^2,'descend');
            
            %index = [1 2];
            
            %Assigning the smallest reflection coefficient(Gamma2) to the user
            %experiencing the worst channel conditions
            pow_vector1(index1(2)) = db2pow(pow_vector1(index1(2)) - av_rx_sir)*Pn;
            pow_vector2(index1(2)) = db2pow(pow_vector2(index1(2)) - av_rx_sir)*Pn;
            pow_vector1(index1(1)) = db2pow(pow_vector1(index1(1)))*Pn;
            pow_vector2(index1(1)) = db2pow(pow_vector2(index1(1)))*Pn;
            [~,index2] = sort(pow_vector2.*abs(CH_vector2).^2,'descend');
            
            y1 = sqrt(pow_vector1(index1(1))).*CH_vector1(index1(1)).*xmod(:,index1(1)) + sqrt(pow_vector1(index1(2))).*CH_vector1(index1(2)).*xmod(:,index1(2));
            y2 = sqrt(pow_vector2(index2(1))).*CH_vector2(index2(1)).*xmod(:,index2(1)) + sqrt(pow_vector2(index2(2))).*CH_vector2(index2(2)).*xmod(:,index2(2));
            
    %         disp(y)
    %         pause
                
            % AWGN
            w1 = sqrt(Pn)*(randn(l_pkt,1)+1i*randn(l_pkt,1))/sqrt(2);
            w2 = sqrt(Pn)*(randn(l_pkt,1)+1i*randn(l_pkt,1))/sqrt(2);
            y1 = y1 + w1; % add the noise
            y2 = y2 + w2;
            %         disp(y)
    %         pause
            
            % RX1 Decoding
            
            %estim_channel = CH_vector(index(1)) + err(index(1));
            
            v = [-1 1];
            v = v(ones(1,l_pkt),:);
            
            %Decoding of the first packet at RX1
    %         tilde_x11 = y1 - sqrt(pow_vector1(index1(1)))*(CH_vector1(index1(1))+err(index1(1)))*v;
    %         tilde_x11 = abs(tilde_x11).^2;
             tilde_x11 = y1 - sqrt(pow_vector1(index1(1)))*(CH_vector1(index1(1))+err(index1(1))).*sn;
             tilde_x11 = abs(tilde_x11).^2;
    %         disp(tilde_x)
    %         pause
            
            [~,ix] = min(tilde_x11,[],2);
            
    %         disp(ix)
    %         pause
            
    %         xdec11 = ix - 1;
              xdec11 = sn(ix);
            
           %Interference Cancellation at RX1
            y1 = y1 - sqrt(pow_vector1(index1(1)))*(CH_vector1(index1(1))+err(index1(1)))*xdec11;
            
            %estim_channel = CH_vector(index(2)) + err(index(2));
            
    %         v = [-1 1];
    %         v = v(ones(1,l_pkt),:);
            
            %Decoding of the second packet at RX1
    %         tilde_x12 = y1 - sqrt(pow_vector1(index1(2)))*(CH_vector1(index1(2))+err(index1(2)))*v;
    %         tilde_x12 = abs(tilde_x12).^2;
            tilde_x12 = y1 - sqrt(pow_vector1(index1(2)))*(CH_vector1(index1(2))+err(index1(2))).*sn;
            tilde_x12 = abs(tilde_x12).^2;
            [~,ix] = min(tilde_x12,[],2);
            
    %         xdec12 = ix - 1;
            xdec12 = sn(ix);
             % RX2 Decoding
            
            %estim_channel = CH_vector(index(1)) + err(index(1));
            
            %Decoding of the first packet at RX2
    %         tilde_x21 = y2 - sqrt(pow_vector2(index2(1)))*(CH_vector2(index2(1))+err(index2(1)))*v;
    %         tilde_x21 = abs(tilde_x21).^2;
            tilde_x21 = y2 - sqrt(pow_vector2(index2(1)))*(CH_vector2(index2(1))+err(index2(1))).*sn;
    %         disp(tilde_x)
    %         pause
            tilde_x21 = abs(tilde_x21).^2;
            [~,ix] = min(tilde_x21,[],2);
            
    %         disp(ix)
    %         pause
            
    %         xdec21 = ix - 1;
            xdec21 = sn(ix);
            
           %Interference Cancellation at RX2
            y2 = y2 - sqrt(pow_vector2(index2(1)))*(CH_vector2(index2(1))+err(index2(1)))*xdec21;
            
            %estim_channel = CH_vector(index(2)) + err(index(2));
            
    %         v = [-1 1];
    %         v = v(ones(1,l_pkt),:);
            
            %Decoding of the second packet at RX2
    %         tilde_x22 = y2 - sqrt(pow_vector2(index2(2)))*(CH_vector2(index2(2))+err(index2(2)))*v;
    %         tilde_x22 = abs(tilde_x22).^2;
            tilde_x22 = y2 - sqrt(pow_vector2(index2(2)))*(CH_vector2(index2(2))+err(index2(2))).*sn;
            tilde_x22 = abs(tilde_x22).^2;
            [~,ix] = min(tilde_x22,[],2);
            
    %         xdec22 = ix - 1;
            xdec22 = sn(ix);
            
            
    %%    Symbol Comparison
    
    % Single RX
        if not(isequal(xdec11,xmod(:,index1(1)))) % first bit in error
                
                errors(index1(1)) = errors(index1(1)) + 1;
                
                conditioned_errors(index1(1)) = conditioned_errors(index1(1)) + 1;
                
            else
                
                succ_bit(index1(1)) = succ_bit(index1(1)) + 1;
                
        end
        if not(isequal(xdec12,xmod(:,index1(2))))%xdec ~= x(index(2))
                
                errors(index1(2)) = errors(index1(2)) + 1;
                
                conditioned_errors(index1(2)) = conditioned_errors(index1(2)) + 1;
                
            else
                
                succ_bit(index1(2)) = succ_bit(index1(2)) + 1;
                
            end
    % Diversity
    
        if(index1==index2) %Stessi coefficienti assegnati ai tag per i due canali
    %     sir1(i) = (abs(CH_vector1(index1(1)))^2)/abs(CH_vector1(index1(2)))^2;
    %     sir2(i) = (abs(CH_vector1(index2(1)))^2)/abs(CH_vector1(index2(2)))^2;
            if not(isequal(xdec11,xmod(:,index1(1)))|isequal(xdec21,xmod(:,index1(1)))) % first bit in error
                %Entrambi i ricevitori decodificano male il primo bit
                
                errors_div(index1(1)) = errors_div(index1(1)) + 1;
                
                conditioned_errors(index1(1)) = conditioned_errors(index1(1)) + 1;
                
            else
                
                succ_bit_div(index1(1)) = succ_bit_div(index1(1)) + 1;
                
            end
            
           
            
            if not(isequal(xdec12,xmod(:,index1(2)))|isequal(xdec22,xmod(:,index1(2)))) % second bit in error
                %Entrambi i ricevitori decodificano male il secondo bit
                errors_div(index1(2)) = errors_div(index1(2)) + 1;
                
                conditioned_errors(index1(2)) = conditioned_errors(index1(2)) + 1;
                
            else
                
                succ_bit_div(index1(2)) = succ_bit_div(index1(2)) + 1;
                
            end
        else %Il primo tag ha canale migliore su uno e peggiore sull'altro
        
            if not(isequal(xdec11,xmod(:,index1(1)))|isequal(xdec22,xmod(:,index2(2)))) % first bit in error
                %Entrambi i ricevitori decodificano male il primo bit
                errors_div(index1(1)) = errors_div(index1(1)) + 1;
                
                conditioned_errors(index1(1)) = conditioned_errors(index1(1)) + 1;
                
            else
                
                succ_bit_div(index1(1)) = succ_bit_div(index1(1)) + 1;
                
            end
            
           
            
            if not(isequal(xdec12,xmod(:,index1(2)))|isequal(xdec21,xmod(:,index2(1)))) % second bit in error
                
                %Entrambi i ricevitori decodificano male il secondo bit
                errors_div(index1(2)) = errors_div(index1(2)) + 1;
                
                conditioned_errors(index1(2)) = conditioned_errors(index1(2)) + 1;
                
            else
                
                succ_bit_div(index1(2)) = succ_bit_div(index1(2)) + 1;
                
            end
    
        end
    
            
        end
        
        bit_count = bit_count + 1;
        
        if errors_div(1) == max_num_err || errors_div(2) == max_num_err
            break
        end
        
        %     disp(errors(i,:))
        %     pause
    end
    
    ber1_div = errors_div(1)/backscattered_bits(1);%/bit_count;
    ber2_div = errors_div(2)/backscattered_bits(2);%/bit_count;
    
    ber1 = errors(1)/backscattered_bits(1);%/bit_count;
    ber2 = errors(2)/backscattered_bits(2);%/bit_count;
    
    conditioned_ber1 = conditioned_errors(1)/bit_both_powered;
    conditioned_ber2 = conditioned_errors(2)/bit_both_powered;
    
    powering_rate1 = bit_not_powered(1)/bit_count;
    powering_rate2 = bit_not_powered(2)/bit_count;
    
    %SNR = 10*log10(mean(snr));
    %Ps1 = Ptx*PL_vector(1)*GAMMA(1);
    
    
    SNR1 = 1; %pow2db(mean(cond_snr1));
    SNR2 = 1; %pow2db(mean(cond_snr2));
    
    POW_DIFF = 1; %abs(SNR1 - SNR2);
    
    cumulative_bits_div = (succ_bit_div(1) + succ_bit_div(2))/(2*bit_count);%(succ_bit(1) + succ_bit(2))/(2*bit_count-bit_not_powered(1)-bit_not_powered(2));
    cumulative_bits = (succ_bit(1) + succ_bit(2))/(2*bit_count);%(succ_bit(1) + succ_bit(2))/(2*bit_count-bit_not_powered(1)-bit_not_powered(2));
    
    end
    
    %Cicla sui 100000 bit, per i diversi 5 parametri del channel gain.
    %Restituisce la media dei bit decodificati con successo dal primo
    %ricevitore. 
    %Restituisce la media dei bit decodificati con successo da uno dei due o da
    %entrambi i ricevitori.