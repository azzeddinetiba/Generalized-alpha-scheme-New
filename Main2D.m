clear;
close all ;
home;
%------------------------
% Les entrées
%------------------------
valid=input('Choisir (Validation du code 1 / Votre propre problème 2) :  ');
t=input('Choisir l instant final :  ');
h=input('Choisir le pas du temps :  ');
if valid==2
Young=input('Young Modulus (Pa) :  ');
Poisson=input('Poisson :  ');
pho=input('La densité (Kg/m3) :  ');
H=input('L hauteur de la plaque (m) :  ');
D=(Young*H^3)/(12*(1-Poisson^2));
F=input('entrez la pression appliquée (Pa) : ','s');
f=str2func(['@(x,y,t) ' F]);
elseif valid==1
    D=1;
    pho=1;
    H=1;
end
Ngauss=input('Choisir le nombre de points de Gauss : ');
npx=input('Entrer nombre de points suivant (Ox) :  ');
npy=input('Entrer nombre de points suivant (Oy) :  ');
AFFICHAGET=input('Affichage de visualisation temporelle avec ou sans maillage y/n : ','s');
phoalpha=input('Choisir la valeur de rho alpha :  ');
alphaf2=phoalpha/(phoalpha+1);
alpham2=(2*phoalpha-1)/(phoalpha+1);
Beta2=0.25*(1-alpham2+alphaf2)^2;
Gamma2=0.5-alpham2+alphaf2;
alpha=input('alpha :  ');
Beta=0.25*(1+alpha)^2;
Gamma=0.5+alpha;
%------------------------------
% Les données du maillage
if valid==2
x1=input('Entrer la coordonée 1 de laxe X :  ');
x2=input('Entrer la coordonée 2 de laxe X :  ');
y1=input('Entrer la coordonée 1 de laxe Y :  ');
y2=input('Entrer la coordonée 2 de laxe Y :  ');
elseif valid==1
x1=0; x2=1;
y1=0; y2=1;
end
elem=1;
k=2;
%--------------------------------
%Compteur du temps de calcul
tic  

% maillage
[X,T,b] = Maillage2d(x1,x2,y1,y2,npx,npy,elem,k);
%--------------------------------
% Initialisation
NT=floor(t/h);
SOL=zeros(size(X,1),NT);
G=zeros(size(X,1),NT);
[Xgauss,Wgauss]=Quadrature(elem, Ngauss);
[phi,grad,deriv2,derivst] = FoncChap(k,Xgauss);
Pi=pi;

if valid==1
    f=@(x,y,t) 2*exp(-t)*sin(Pi*x)^2*sin(Pi*y)^2-4*t*exp(-t)*sin(Pi*x)^2*sin(Pi*y)^2+t^2*exp(-t)*sin(Pi*x)^2*sin(Pi*y)^2+24*t^2*exp(-t)*sin(Pi*x)^2*Pi^4*sin(Pi*y)^2-16*t^2*exp(-t)*cos(Pi*x)^2*Pi^4*sin(Pi*y)^2+8*t^2*exp(-t)*cos(Pi*x)^2*Pi^4*cos(Pi*y)^2-16*t^2*exp(-t)*sin(Pi*x)^2*cos(Pi*y)^2*Pi^4;
     %f=@(x,y,t) -225*sin(15*t)*sin(Pi*x)^2*sin(Pi*y)^2+24*sin(15*t)*sin(Pi*x)^2*Pi^4*sin(Pi*y)^2-16*sin(15*t)*cos(Pi*x)^2*Pi^4*sin(Pi*y)^2+8*sin(15*t)*cos(Pi*x)^2*Pi^4*cos(Pi*y)^2-16*sin(15*t)*sin(Pi*x)^2*cos(Pi*y)^2*Pi^4;
end

[K,M] = MEF2D(f,elem,k,X,T,b,Ngauss,h,D,pho,H);
[a,v,d] = HHTmethod(K,M,f,NT,h,alphaf2,alpham2,Beta2,Gamma2,elem,k,X,T,b,Ngauss,D,pho,H,alpha,Beta,Gamma);




for i=1:NT
    t=h*(i-1);
    U=d(:,i);
%------------------------------
% Affichage et comparaison
if valid==1
% Solution approchée
figure(1); clf; 
colormap hsv;
trisurf(T,X(:,1),X(:,2),U,'edgecolor','k','facecolor','interp');
view(2),axis([x1 x2 y1 y2]),axis equal; colorbar; title(['Solution approchee à t=',num2str(t),' s']);
%Solution exacte
figure(2); clf;
colormap hsv;
uee=(t^2)*exp(-t)*(sin(pi*X(:,1)).^2).*sin(pi*X(:,2)).^2;
trisurf(T,X(:,1),X(:,2),uee,'edgecolor','k','facecolor','interp');
view(2),axis([x1 x2 y1 y2]),axis equal; colorbar;  title(['Solution exacte à t=',num2str(t),' s']);
%--------------------------------
% Calcul de l'erreur en norme L2
error_L2=L2_Error(U,uee,X,T,Ngauss,elem,k);
error_H1=H1_Error(U,uee,X,T,Ngauss,elem,k,error_L2);
fprintf('at t= %d s \t\t %4d \t\t %20.16e \t\t %20.16e\n\n',t, npx*npy,error_L2,error_H1);
end
end
toc
%--------------------------------

figure('name', 'Deplacement transversal de la plaque');
    colormap hsv;    
    
    if AFFICHAGET=='y'
    for i=1:size(d,2)
    trisurf(T,X(:,1),X(:,2),d(:,i),'edgecolor','k','facecolor','interp');
    colorbar;
    axis off; title( sprintf('t = %f s', (i-1)*h) );
    pause(0.2)
    end
    
    
    elseif AFFICHAGET=='n'
    if rem(npx,2)==0
        npx=npx+1;
    end
    if rem(npy,2)==0
        npy=npy+1;
    end
    for i=1:size(d,2)
    post=d(:,i);
    AFF=zeros(npx,npy);  
         for j=1:npy      
             AFF(:,j)=post(npx*(j-1)+1:npx*j);
         end
    surf(linspace(x1,x2,npx),linspace(y1,y2,npy),AFF','FaceColor','interp','EdgeColor','none');
    colorbar; 
    light
    lighting gouraud
    material dull
    title( sprintf('t = %f s', (i-1)*h) );
    pause(0.1)
    end
    end
    
   