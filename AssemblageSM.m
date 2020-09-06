function F = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps)
%Assemblage2D fonction pour l'assemblage des matrices K,M et F �lement par �lement 
%
%  Assemblage des matrices �l�mentaires "Ke" dans la matrice globale K
%  Assemblage des matrices �l�mentaires "Me" dans la matrice globale M
%  Assemblage des seconds membres �l�mentaires "Fe" dans le second membre global F
%  L'assembalge s'effectue �lement par �lement
% 
% Sorties:
%          K: Matrice de raideur apr�s assemblage
%          M: Matrice de masse apr�s assemblage
%          F: Matrice second membre apr�s assemblage
%     

Nn=size(X,1);     % nombre des noeuds
Nt=size(T,1);     % nombre des �l�ments
%K=sparse(zeros(Nn,Nn));% initialisation de K comme matrice creuse
F=zeros(Nn,1);   % initialisation de F

 for i = 1:Nn
    F(i)=f(X(i,1),X(i,2),temps);
 end 
 
 F(b)=0;
%K=sparse(K);
end

