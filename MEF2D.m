function [K,M] = MEF2D(f,elem,k,X,T,b,Ngauss,h,D,pho,H)
%MEF2D r�sout la solution de l'�quation biharmonique dans le 1er instant
% moyennant la m�thode des �l�ments finis de type Galerkin Pk
%
% Sorties:
%          Bb:  Matrice globale constitu�e de K et M
%          Fb:  Matrice globale second membre
%          U:  Le vecteur solution
%
% --------Quadrature de Gauss-----------
[Xgauss,Wgauss]=Quadrature(elem, Ngauss);
%
%---------Fonctions chapeaux------------
if elem==1
[phi,grad,deriv2,derivst] = FoncChap(k,Xgauss);
elseif elem==0
[phi,grad,deriv2,derivst] = FoncChapQ(k,Xgauss);
end
%
%-------- Assemblage--------------------
[ K,M ] = Assemblage2D(elem,k,X,T,f,phi,grad,deriv2,derivst,Xgauss,Wgauss,Ngauss,h,pho,H);
%[ B,F ] = GlobalAssembly(K,M,FU,h,D,pho,H);
%

end