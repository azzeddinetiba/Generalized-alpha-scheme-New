function [ K,M] = Assemblage2D(elem,k,X,T,f,phi,grad,deriv2,derivst,Xgauss,Wgauss,Ngauss,h,pho,H)
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

