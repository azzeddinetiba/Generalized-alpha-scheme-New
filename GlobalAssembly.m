function [ B,F ] = GlobalAssembly(K,M,FU,h,D,pho,H)
%GlobalAssembly La fonction qui crée les matrices globales B.U=F
%
% Création de la matrice globale constituée de K, M et M/h^2
% Création de la matrice globale second membre à partir de l'assemblage
% 
% Sorties:
%          B: Matrice Globale 
%          F: Matrice second membre globale tel que B.U=F
%          
% 
FG=zeros(size(K,1),1);

B=[D*K M;
   pho*H*M/(h^2) -K];

F=[FG;FU];
end

