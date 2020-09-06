function [X,T,b] = Maillage2d(x1,x2,y1,y2,npx,npy,elem,k)
%Maillag2d crée un maillage uniforme sur un domaine 
% rectangulaire [x1,x2]x[y1,y2]
%
% Entrées:    
%   x1,x2,y1,y2:    coordonnées des coins 
%   npx,npy:        nombres de noeuds dans chaque direction 
% Sorties:   
%          X:  table des coordonnées
%          T:  table de connectivité
%          b:  noeuds sur le bords
    if elem==1

    if k==1 
    
X = zeros((npx)*(npy),2);
% Nombre des éléments dans chaque direction
nx = npx-1; ny = npy-1;
%

xs = linspace(x1,x2,npx)'; 
unos = ones(npx,1);
%
% Coordonnées des noeuds
%
yys = linspace(y1,y2,npy);
for i=1:npy
    ys = yys(i)*unos; 
    posi = (i-1)*(npx)+1:i*(npx); 
    X(posi,:)=[xs,ys];
end
% P1
% Connectivité
T = zeros(nx*ny,3);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+(npx)];
                T(ielem+1,:) = [inode+1   inode+1+npx   inode+npx];
            end
        end
%Noeuds sur le bord stockées dans "b" comme suit (bas,gauche,droite,haut)
b=[1:npx,npx+1:npx:npx*npy,2*npx:npx:npx*npy,npx*npy-npx+2:npx*npy-1];

    elseif k==2
        
        if rem(npx,2)==0 
            npx=npx+1;
        end
        if rem(npy,2)==0 
            npy=npy+1;
        end
        
b=[1:npx,npx+1:npx:npx*npy,2*npx:npx:npx*npy,npx*npy-npx+2:npx*npy-1];


        X = zeros((npx)*(npy),2);
% Nombre des éléments dans chaque direction
%

xs = linspace(x1,x2,npx)'; 
unos = ones(npx,1);
%
% Coordonnées des noeuds
%
yys = linspace(y1,y2,npy);
for i=1:npy
    ys = yys(i)*unos; 
    posi = (i-1)*(npx)+1:i*(npx); 
    X(posi,:)=[xs,ys];
end

% P2
% Connectivité
       nx=(npx-1)/2; 
       ny=(npy-1)/2;
       
        N=nx*ny*2;
        n=ny;
        
T = zeros(N,6);
        for i=1:n
            for j=((2*i-2)*npx+1):2:npx+((2*i-2)*npx+1)-3
                J=j-(i-1)*npx-i+1;
                
           T(J,:)=[j j+1 j+2 j+npx+1 j+2*npx j+npx];
           T(J+1,:)=[j+2 j+npx+2 j+2*npx+2 j+2*npx+1 j+2*npx j+npx+1];
            end
        end
        
    end
    end


if elem==0
    if k==1
    X = zeros((npx)*(npy),2);
% Nombre des éléments dans chaque direction
nx = npx-1; ny = npy-1;
%

xs = linspace(x1,x2,npx)'; 
unos = ones(npx,1);
%
% Coordonnées des noeuds
%
yys = linspace(y1,y2,npy);
for i=1:npy
    ys = yys(i)*unos;
    posi = (i-1)*(npx)+1:i*(npx);
    X(posi,:)=[xs,ys];
end
% P1
% Connectivité
T = zeros(nx*ny,4);

for j=1:ny
for i=1:nx
    I=npx*j-nx+i-1;
    T(I-j+1,1:2)=[I,I+1];
    T(I-j+1,3:4)=[I+npx+1,I+npx];
end
end

%Noeuds sur le bord stockées dans "b" comme suit (bas,gauche,droite,haut)
b=[1:npx,npx+1:npx:npx*npy,2*npx:npx:npx*npy,npx*npy-npx+2:npx*npy-1];

    
elseif k==2
        X = zeros((npx)*(npy),2);
% Nombre des éléments dans chaque direction
nx =npx-1; ny =npy-1;
%
hx = (x2-x1)/nx;
hy = (y2-y1)/ny;
xs = linspace(x1,x2,npx)'; 
unos = ones(npx,1);
%
% Coordonnées des noeuds
%
yys = linspace(y1,y2,npy);
for i=1:npy
    ys = yys(i)*unos;
    posi = (i-1)*(npx)+1:i*(npx);
    X(posi,:)=[xs,ys];
end

XP=zeros((2*npx-1)*(2*npy-1),2);

for j=1:2:2*npy-1
    J=(2*npx-1)*j-(2*npx-1-1);
    JPR=(j+1)/2;
    XP(J:2:J-1+(2*npx-1),:)=X(JPR*npx-npx+1:JPR*npx,:);
    XP(J+1:2:J-1+(2*npx-2),:)=X(JPR*npx-npx+1:JPR*npx-1,:)+(hx/2)*[ones(npx-1,1),zeros(npx-1,1)];
    if j==2*npy-1
        break;
    end
    XP(J-1+(2*npx-1)+1:J-1+(2*npx-1)+1+(2*npx-1-1),:)=XP(J:J-1+(2*npx-1),:)+(hy/2)*[zeros(2*npx-1,1),ones(2*npx-1,1)]; 
end
%%%%%%  :") %%%%%%

    X=XP;
    
% P1
% Connectivité
T = zeros(nx*ny,9);

for j=1:ny
for i=1:nx
    I=npx*j-nx+i-1;
    Fill=2*i-1+(2*(2*npx-1))*(j-1);
    T(I-j+1,1:3)=[Fill,Fill+1,Fill+2];
    T(I-j+1,4:6)=[Fill+(2*npx-1),Fill+(2*npx-1)+1,Fill+(2*npx-1)+2];
    T(I-j+1,7:9)=[Fill+2*(2*npx-1),Fill+2*(2*npx-1)+1,Fill+2*(2*npx-1)+2];
end
end

%Noeuds sur le bord stockées dans "b" comme suit (bas,gauche,droite,haut)
npx=2*npx-1; npy=2*npy-1;
b=[1:npx,npx+1:npx:npx*npy,2*npx:npx:npx*npy,npx*npy-npx+2:npx*npy-1];
    end
end

end   


