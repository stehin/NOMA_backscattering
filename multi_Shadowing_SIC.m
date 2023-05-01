function [ber1_div,ber2_div,SNR1,SNR2,ber1,ber2,powering_rate1,powering_rate2,cumulative_bits_div,cumulative_bits] = multi_Shadowing_SIC(scenario,alpha,av_rx_snr1,av_rx_snr2,av_rx_sir,estim_err_var)

    % System parameters
    W = scenario.bandwidth; % available bandwidth [Hz]
    NodBm = scenario.noise_spectral_density; %noise power spectral density [dBm/Hz]
    No = db2pow(NodBm - 30); % noise power spectral density [W/Hz]
    Pn = No*W; % noise power [W]
    n_iter = scenario.niter;
    max_num_err = scenario.maxerrs;
    l_pkt = scenario.packet_length;
    phase_flag = scenario.phase_off;
    
    
        
    sigma_nakagami1f = scenario.link1.sigma_nakagamif; %sigma della Nakagami
    sigma_nakagami1b = scenario.link1.sigma_nakagamib; 
    
    
    sigma_nakagami2f = scenario.link2.sigma_nakagamif;
    sigma_nakagami2b = scenario.link2.sigma_nakagamib;
    
    
    mf_nakagami1 = scenario.link1.mf_nakagami; %m della Nakagami
    mb_nakagami1 = scenario.link1.mb_nakagami;
    
    mf_nakagami2 = scenario.link2.mf_nakagami;
    mb_nakagami2 = scenario.link2.mb_nakagami;
    
    pdf1 = makedist('Nakagami','mu',mf_nakagami1,'omega',sigma_nakagami1f); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
    pdb1 = makedist('Nakagami','mu',mb_nakagami1,'omega',sigma_nakagami1b); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
    pdf2 = makedist('Nakagami','mu',mf_nakagami2,'omega',sigma_nakagami2f); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
    pdb2 = makedist('Nakagami','mu',mb_nakagami2,'omega',sigma_nakagami2b); %Crea la PDF con distribuzione Nakagami, con i parametri passati.
    
    %se è presente lo shadowing recupero le matrici di correlazione
    if(scenario.shadowing)
        Sigma_b1 = scenario.link1.Sigmab;
        Sigma_b2 = scenario.link2.Sigmab;
    end
    
    %Inizializzazione a zero delle variabili da calcolare
    errors_div = zeros(1,2);
    errors = zeros(1,2);
    backscattered_bits = zeros(1,2);
    conditioned_errors = zeros(1,2);
    succ_bit_div = zeros(1,2);
    succ_bit = zeros(1,2);
    bit_both_powered = 0;
    bit_not_powered = zeros(1,2);
    bit_count = 0;
    
    
    if(scenario.HOM > 2)
        MQAM = scenario.HOM;
        d = sqrt(6/(MQAM-1));
        n0 =  (sqrt(MQAM)/2-0.5);
        count=1;
        %Genero la costellazione di una MQAM e la salvo come vettore riga
        for nn=0:(sqrt(MQAM)-1)
            for mm = 0:(sqrt(MQAM)-1)
                sn(count,1) = (nn - n0)*d + 1i*(mm - n0 )*d;
                count = count + 1;
            end
        end
        sn = sn';
    else
        %Genero la costellazione binaria
        sn = [-1 1];
    end
    for i=1:n_iter
        
        % generate channel gain for the forward link1
        af1 = random(pdf1,1,1); %Modulo del forward link1, ho deciso che segue la Nakagami
        phif1 = phase_flag*2*pi*rand(1,1); %Fase uniformemente distribuita slide 59
        hf1 = af1.*exp(1i*phif1); %slide 59, def dell'inviluppo del segnale ricevuto
        
        % generate channel gain for the forward link2
        af2 = random(pdf2,1,1);    
        phif2 = phase_flag*2*pi*rand(1,1);
        hf2 = af2.*exp(1i*phif2);
        
        % generate channel gain for the reverse link (dal tag al ricevitore. 
        % Ho due tag e due ricevitori, 4 backward links)
        ab11 = random(pdb1,1,1); %Modulo del channel gain Backward11, 
        % generato con stessa Nakagami pdb1, usata anche che per ab12
        ab12 = random(pdb1,1,1); %Generati con stessa Nakagami pdb1
        
        phib11 = phase_flag*2*pi*rand(1,1);  %Fase del channel gain Backward11, 
        %generato con distribuzione uniforme, nulla se phase_offs=0
        phib12 = phase_flag*2*pi*rand(1,1);
        
        hb11 = ab11.*exp(1i*phib11);
        hb12 = ab12.*exp(1i*phib12);
        
        % generate channel gain for the reverse link
        ab21 = random(pdb2,1,1);
        ab22 = random(pdb2,1,1);
        phib21 = phase_flag*2*pi*rand(1,1);
        phib22 = phase_flag*2*pi*rand(1,1);
        hb21 = ab21.*exp(1i*phib21);
        hb22 = ab22.*exp(1i*phib22);
        
        %Unisco i contributi di forward e backward link
        h11 = hf1*hb11;
        h12 = hf2*hb12;
        h21 = hf1*hb21;
        h22 = hf2*hb22;
    
        %Aggiunta Shadowing, solo nel backward
        if(scenario.shadowing)
            X_b1 = mvnrnd([0 0],Sigma_b1); %genera una v.a. gaussiana bidimensionale,
            %a media nulla, a partire dalla matrice di correlazione Sigma_b
            X_b2 = mvnrnd([0 0],Sigma_b2);
            x_s11b = X_b1(1); 
            x_s11 = x_s11b; %nel forward non c'è shadowing
            x_s12b = X_b2(1);
            x_s12 =  x_s12b; %nel forward non c'è shadowing
            x_s21b = X_b1(2);
            x_s21 =  x_s21b; %nel forward non c'è shadowing
            x_s22b = X_b2(2);
            x_s22 = x_s22b; %nel forward non c'è shadowing
            Shadowing1 = [x_s11 x_s12];
            Shadowing1 = exp(Shadowing1);
            Shadowing2 = [x_s21 x_s22];
            Shadowing2 = exp(Shadowing2);
        end
        % generate modulated symbols (BPSK)
        xmod = [datasample(sn,1) datasample(sn,1)]; %prendo due simboli a caso della costellazione
        %forward_ch_vector = [hf1 hf2]; % channel gains of the forward links ?
    
        %Inserisco i channel gains che tengono conto di entrambi i forward e 
        %backward link n due vettori CH_vector1 (RX1) e CH_vector2 (RX2) con o
        %senza shadowing a seconda dello scenario.
        if(scenario.shadowing)
            CH_vector1 = [h11 h12].*Shadowing1; % channel gains of the composite links
            CH_vector2 = [h21 h22].*Shadowing2;
        else
            
            CH_vector1 = [h11 h12]; % channel gains of the composite links
            CH_vector2 = [h21 h22];
        end
        
        %powered_index = [1 2]; %entrambi i tag backscatterano ?
        bit_both_powered = bit_both_powered + 1; %Aggiorno il numero di casi in cui
        %entrambi i bit sono backscatterati
        
        backscattered_bits = backscattered_bits + 1;
        
        %estimation error
        %Nel nostro caso perfetta stima di canale
        err = sqrt(estim_err_var)*(randn(1,2)+1i*randn(1,2))/sqrt(2);
        pow_vector1 = av_rx_snr1; %SNR al RX1, è un vettore
        pow_vector2 = av_rx_snr2; %SNR al RX2, è un vettore
            
        % Ordering the channel gains in order to assign the reflection
        % coefficients
        %Il primo elemento di index1, indica quale tra tag1 e tag2 è il tag con
        %realizzazione migliore del canale composto fino ad RX1
        [~,index1] = sort(abs(CH_vector1).^2,'descend');
        
        %Assigning the smallest reflection coefficient(Gamma2) to the user
        %experiencing the worst channel conditions
        pow_vector1(index1(2)) = db2pow(pow_vector1(index1(2)) - av_rx_sir)*Pn;%a chi ha guadagno minore, toglie SIR
        pow_vector2(index1(2)) = db2pow(pow_vector2(index1(2)) - av_rx_sir)*Pn;
        pow_vector1(index1(1)) = db2pow(pow_vector1(index1(1)))*Pn; %potenza del segnale ricevuto da RX1 dal tag che ha guadagno maggiore
        pow_vector2(index1(1)) = db2pow(pow_vector2(index1(1)))*Pn;
    
        %Il primo elemento di index2, indica quale tra tag1 e tag2 è il tag con
        %realizzazione migliore del canale composto fino ad RX2
        [~,index2] = sort(pow_vector2.*abs(CH_vector2).^2,'descend');
        
        y1 = sqrt(pow_vector1(index1(1))).*CH_vector1(index1(1)).*xmod(:,index1(1)) + sqrt(pow_vector1(index1(2))).*CH_vector1(index1(2)).*xmod(:,index1(2));
        y2 = sqrt(pow_vector2(index2(1))).*CH_vector2(index2(1)).*xmod(:,index2(1)) + sqrt(pow_vector2(index2(2))).*CH_vector2(index2(2)).*xmod(:,index2(2));
    
        % AWGN
        w1 = sqrt(Pn)*(randn(l_pkt,1)+1i*randn(l_pkt,1))/sqrt(2);
        w2 = sqrt(Pn)*(randn(l_pkt,1)+1i*randn(l_pkt,1))/sqrt(2);
        y1 = y1 + w1; % add the noise
        y2 = y2 + w2;
    
        % RX1 Decoding
        %Decoding of the first packet at RX1
        tilde_x11 = y1 - sqrt(pow_vector1(index1(1)))*(CH_vector1(index1(1))+err(index1(1))).*sn;
        tilde_x11 = abs(tilde_x11).^2;
        [~,ix] = min(tilde_x11,[],2);
        xdec11 = sn(ix);
    
        %Interference Cancellation at RX1
        y1 = y1 - sqrt(pow_vector1(index1(1)))*(CH_vector1(index1(1))+err(index1(1)))*xdec11;
    
        %Decoding of the second packet at RX1
        tilde_x12 = y1 - sqrt(pow_vector1(index1(2)))*(CH_vector1(index1(2))+err(index1(2))).*sn;
        tilde_x12 = abs(tilde_x12).^2;
        [~,ix] = min(tilde_x12,[],2);
        xdec12 = sn(ix);
        
        
        % RX2 Decoding
        %Decoding of the first packet at RX2
        tilde_x21 = y2 - sqrt(pow_vector2(index2(1)))*(CH_vector2(index2(1))+err(index2(1))).*sn;
        tilde_x21 = abs(tilde_x21).^2;
        [~,ix] = min(tilde_x21,[],2);
        xdec21 = sn(ix);
            
        %Interference Cancellation at RX2
        y2 = y2 - sqrt(pow_vector2(index2(1)))*(CH_vector2(index2(1))+err(index2(1)))*xdec21;
        
        %Decoding of the second packet at RX2
        tilde_x22 = y2 - sqrt(pow_vector2(index2(2)))*(CH_vector2(index2(2))+err(index2(2))).*sn;
        tilde_x22 = abs(tilde_x22).^2;
        [~,ix] = min(tilde_x22,[],2);
        xdec22 = sn(ix);
            
            
        %Symbol Comparison
    
        %Single RX
        if not(isequal(xdec11,xmod(:,index1(1)))) % first bit in error    
            errors(index1(1)) = errors(index1(1)) + 1;
            conditioned_errors(index1(1)) = conditioned_errors(index1(1)) + 1;
        else
            succ_bit(index1(1)) = succ_bit(index1(1)) + 1;   
        end
        if not(isequal(xdec12,xmod(:,index1(2)))) %second bit in error
            errors(index1(2)) = errors(index1(2)) + 1;
            conditioned_errors(index1(2)) = conditioned_errors(index1(2)) + 1;
        else
            succ_bit(index1(2)) = succ_bit(index1(2)) + 1;
        end
    
        %Diversity
        if(index1==index2) %Su entrambi i canali uno dei due tag ha realizzazioni migliori/peggiori
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
        bit_count = bit_count + 1;
     %esco dal ciclo se raggiungo il numero massimo di errori su uno dei due bit    
        if errors_div(1) == max_num_err || errors_div(2) == max_num_err
            break
        end
    end
    
    ber1_div = errors_div(1)/backscattered_bits(1);%/bit_count;
    ber2_div = errors_div(2)/backscattered_bits(2);%/bit_count;
    
    ber1 = errors(1)/backscattered_bits(1);%/bit_count;
    ber2 = errors(2)/backscattered_bits(2);%/bit_count;
    
    % conditioned_ber1 = conditioned_errors(1)/bit_both_powered;
    % conditioned_ber2 = conditioned_errors(2)/bit_both_powered;
    
    powering_rate1 = bit_not_powered(1)/bit_count;
    powering_rate2 = bit_not_powered(2)/bit_count;
    
    SNR1 = 1; %pow2db(mean(cond_snr1));
    SNR2 = 1; %pow2db(mean(cond_snr2));
    
    %POW_DIFF = 1; %abs(SNR1 - SNR2);
    
    cumulative_bits_div = (succ_bit_div(1) + succ_bit_div(2))/(2*bit_count);%(succ_bit(1) + succ_bit(2))/(2*bit_count-bit_not_powered(1)-bit_not_powered(2));
    cumulative_bits = (succ_bit(1) + succ_bit(2))/(2*bit_count);%(succ_bit(1) + succ_bit(2))/(2*bit_count-bit_not_powered(1)-bit_not_powered(2));
    
    end
    
    %Cicla sui 100000 bit, per i diversi 5 parametri del channel gain.
    %Restituisce la media dei bit decodificati con successo dal primo
    %ricevitore. 
    %Restituisce la media dei bit decodificati con successo da uno dei due o da
    %entrambi i ricevitori.