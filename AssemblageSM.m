function F = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps)
%Assemblage2D fonction pour l'assemblage des matrices K,M et F élement par élement 
%
%  Assemblage des matrices élémentaires "Ke" dans la matrice globale K
%  Assemblage des matrices élémentaires "Me" dans la matrice globale M
%  Assemblage des seconds membres élémentaires "Fe" dans le second membre global F
%  L'assembalge s'effectue élement par élement
% 
% Sorties:
%          K: Matrice de raideur après assemblage
%          M: Matrice de masse après assemblage
%          F: Matrice second membre après assemblage
%     

Nn=size(X,1);     % nombre des noeuds
Nt=size(T,1);     % nombre des éléments
%K=sparse(zeros(Nn,Nn));% initialisation de K comme matrice creuse
F=zeros(Nn,1);   % initialisation de F

 for i = 1:Nn
    F(i)=f(X(i,1),X(i,2),temps);
 end 
 
 F(b)=0;
%K=sparse(K);
end

