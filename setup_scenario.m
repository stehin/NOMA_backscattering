%La function restituisce gli N angoli sottesi ai due vettori (A e B), la 
% distanza d tra i due punti E ed F (Tag, RX) e la distanza h tra i due punti F e O
%(origine, dove si trova il carrier emitter) (Tag, CE). Al suo interno calcola le diverse
%distanze e gli N angoli richiamando trigo(F,E). Riceve in input una struct
%contenente  un vettore per le coordinate (cartesiane (x,y)) del tag1, un 
% vettore per le coordinate (cartesiane (x,y)) di tag2, un vettore per le 
% coordinate (polari (rho, theta)) di RX1 e un vettore per le coordinate 
% (polari (rho, theta)) di RX2. Restituisce in output una struct
% contenente un vettore per gli angoli Theta_12_1, un vettore per gli
% angoli Theta_12_2, un vettore contenente le distanze distanze tag1-CE
% (df1), un vettore contenente le distanze distanze tag2-CE (df2), un 
% vettore contenente le distanze distanze tag1-RX1 (d11), un 
% vettore contenente le distanze distanze tag1-RX2 (d21), un 
% vettore contenente le distanze distanze tag2-RX1 (d12), un 
% vettore contenente le distanze distanze tag2-RX2 (d22).
function [setup] = setup_scenario(parametri_scenario)
    %Ho un'iterazione del ciclo per ogni diverso scenario. Al termine di ogni
    %iterazione i vettori d11, d21, d12, d22, Theta_12_1 e Theta_12_2 sono
    %aggiornati.
    
        %Inizializzo i vettori con valori nulli
        d11=zeros(1,length(parametri_scenario.coordinate_tag1));
        d12=zeros(1,length(parametri_scenario.coordinate_tag1));
        d21=zeros(1,length(parametri_scenario.coordinate_tag1));
        d22=zeros(1,length(parametri_scenario.coordinate_tag1));
        Theta_12_1=zeros(1,length(parametri_scenario.coordinate_tag1));
        Theta_12_2=zeros(1,length(parametri_scenario.coordinate_tag1));
        df1=zeros(1,length(parametri_scenario.coordinate_tag1));
        df2=zeros(1,length(parametri_scenario.coordinate_tag1));
    
        for i=1:length(parametri_scenario.coordinate_tag1)
            pos_tag1=parametri_scenario.coordinate_tag1(i,:);
            pos_tag2=parametri_scenario.coordinate_tag2(i,:);
            pos_rx1=parametri_scenario.coordinate_rx1(i,:);
            pos_rx2=parametri_scenario.coordinate_rx2(i,:);
            [Theta_12_1_l, d_11, df_1a] = trigo(pos_tag1, pos_rx1);
            [Theta_12_1_r, d_21, df_1b] = trigo(pos_tag1, pos_rx2);
            [Theta_12_2_l, d_12, df_2b] = trigo(pos_tag2, pos_rx1);
            [Theta_12_2_r, d_22, df_2a] = trigo(pos_tag2, pos_rx2);
            Theta_12_1_tot=Theta_12_1_l+Theta_12_1_l;
            Theta_12_2_tot=Theta_12_2_l+Theta_12_2_l;
            d11(i)=d_11;
            d12(i)=d_12;
            d21(i)=d_21;
            d22(i)=d_22;
            df1(i)=df_1a; %con df_1b è lo stesso
            df2(i)=df_2a; %con df_2b è lo stesso
            Theta_12_1(i)=Theta_12_1_tot;
            Theta_12_2(i)=Theta_12_2_tot;
    
        end
        setup.d11=d11;
        setup.d21=d21;
        setup.d12=d12;
        setup.d22=d22;
        setup.df1=df1;
        setup.df2=df2;
        setup.Theta_12_1= Theta_12_1;
        setup.Theta_12_2= Theta_12_2;
    