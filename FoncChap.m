function [phi,phi1deriv,phi2deriv,phiderivst] = FoncChap(k,Xgauss)
%FoncChap donne les fonctions de formes dans l'element de référence
%
%
%        N5=(0,1)
%               |\
%               | \
%               |  \
%               |   \
%               |    \
%            N6 |     \ N4
%               |      \
%               |  TR   \
%               |        \
%               |_ _ _ _ _\
%        N1=(0,0)    N2    (1,0)=N3
%
%
%       
% Sorties:
%         phi: tableau dont les colonnes contiennent les fonctions chapeaux 
%              phi1, phi2 et phi3... phi6 dans cet ordre
%         phi1deriv: tableau dont les colonnes contiennent les gradiants  
%              des fonctions chapeaux phi1, phi2 et phi3... phi6 dans cet ordre
%         phi2deriv: tableau dont les colonnes contiennent les dérivées  
%              secondes des fonctions chapeaux phi1, phi2 et phi3... phi6 dans cet ordre
%         phi1deriv: tableau dont les colonnes contiennent les derivée croisées  
%              des fonctions chapeaux phi1, phi2 et phi3... phi6 dans cet ordre
%

s = Xgauss(:,1); 
t = Xgauss(:,2); 

if k==2
    phi=[2*s.^2+4*s.*t+2*t.^2-3*s-3*t+1,-4*s.^2-4*s.*t+4*s,2*s.^2-s,4*s.*t,2*t.^2-t,-4*s.*t-4*t.^2+4*t];
    phi1deriv=[4*s+4*t-3,-8*s-4*t+4,4*s-1,4*t,zeros(size(s,1),1),-4*t;
        4*s+4*t-3,-4*s,zeros(size(s,1),1),4*s,4*t-1,-4*s-8*t+4];
    phi2deriv=[4,-8,4,0,0,0;4,0,0,0,4,-8];
    phiderivst=[4,-4,0,4,0,-4];
end   
end
    


