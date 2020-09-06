function [Fb,U] = MEF2Ddyna(Bb,f,elem,k,X,T,b,Ngauss,h,Upre1,Upre2,pho,H,temps)
%MEF2Ddyna r�sout la solution de l'�quation biharmonique dans les instants ult�rieurs
% moyennant la m�thode des �l�ments finis de type Galerkin Pk
% La r�solution des instants ult�rieurs ne requiert pas toute la proc�dure
% �lement finis mais seulement le calcul du vecteur second membre du nouvel
% instant et l'inversion de la matrice puisque le sch�ma adopt� est un
% sch�ma implicite

Nn=size(X,1); %nombre des noeuds
Nt=size(T,1);  % nombre des �l�ments
FU=zeros(Nn,1);% initialisation de F comme matrice creuse

for ie = 1:Nt
    Tie =T(ie,:);
    [Fe] = SMelem(f,elem,k,X,T,ie,Ngauss,h,Upre1,Upre2,pho,H,temps);
    FU(Tie)=FU(Tie)+Fe;
    clear Fe; 
end 
FG=zeros(Nn,1);
F=[FG;FU];
%--------Conditions aux bords-----------
%On met "0" dans les lignes et les colonnes qui correspondent aux noeuds
%bords dans K et F
F(ones(1,size(b,2))*size(X,1)+b)=0;
%On met "1" sur la diagonale de K correspodant aux aux noeuds bords 
%K(b,b)=speye(length(b),length(b)); % cas o� K est sparse
% K et F apr�s conditions aux bords
%
Fb=F;
U=Bb\Fb;
end