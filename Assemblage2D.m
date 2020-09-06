function [ K,M] = Assemblage2D(elem,k,X,T,f,phi,grad,deriv2,derivst,Xgauss,Wgauss,Ngauss,h,pho,H)
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
K=zeros(Nn,Nn);   % initialisation de K 
M=zeros(Nn,Nn);   % initialisation de M

 for ie = 1:Nt
    Tie =T(ie,:);
    [Ke,Me] =MatElem2D(elem,k,X,T,ie,phi,grad,deriv2,derivst,Xgauss,Wgauss);
    K(Tie,Tie)=K(Tie,Tie)+Ke;
    M(Tie,Tie)=M(Tie,Tie)+Me;    
    clear Ke Me; 
 end 
%K=sparse(K);
end

