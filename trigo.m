%La function restituisce l'angolo sotteso ai due segmenti, la distanza d
%tra i due punti pos_rx e pos_tag (Tag, RX) e la distanza h tra i due punti 
% pos_tag e O (origine, dove si trova il carrier emitter) (Tag, CE). 
% Vuole in input le coordinate polari di pos_rx (RX), r quelle cartesiane 
% di pos_tag (tag).
function [theta_rad, d, h] = trigo(pos_tag, pos_rx)
    %Converto pos_rxda coordinate polari in coordinate cartesiane
    temp=[pos_rx(1)*cos(pos_rx(2)) pos_rx(1)*sin(pos_rx(2))];
    pos_rx=temp; 
    O=[0 0];
    A=O-pos_tag;
    B=pos_rx-pos_tag;
    mod_A=sqrt(dot(A,A));
    mod_B=sqrt(dot(B,B));
    theta_rad=acos(dot(A,B)/(mod_A*mod_B));
    %theta_grad=theta_rad*360/(2*pi);
    h=mod_A;
    d=mod_B;
end