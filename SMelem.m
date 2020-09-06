function [ Fe ] = SMelem(f,elem,k,X,T,ie,Ngauss,h,pho,H,temps)
%SMelem calcul du second membre élémentaire Fe
%

[Xgauss,Wgauss]=Quadrature(elem, Ngauss);
%fonctions chapeaux
Ngeom=6;


        phi = FoncChap(k,Xgauss);

    
x=X(T(ie,1:Ngeom),:);
s = Xgauss(:,1); 
t = Xgauss(:,2); 
    
Fe=zeros(Ngeom,1); %initialiser les tailles de Me et Fe 

%les coefficients abcdef tel que x=a*s²+b*t²+c*st+d*s+e*t+f
%abcdef a deux colonnes la première de x vers (s,t) et la deuxième
%de y vers (s,t)
abcdef=[2*x(1,:)-4*x(2,:)+2*x(3,:);
        2*x(1,:)+2*x(5,:)-4*x(6,:);
        4*x(1,:)-4*x(2,:)+4*x(4,:)-4*x(6,:);
        -3*x(1,:)+4*x(2,:)-x(3,:);
        -3*x(1,:)-x(5,:)+4*x(6,:);
        x(1,:)];
    


% matrice jacobienne de taille [2*Nombre de points de Gauss,2]
J=[(2*abcdef(1,1)*s+abcdef(3,1)*t+abcdef(4,1)),(2*abcdef(2,1)*t+abcdef(3,1)*s+abcdef(5,1));
    (2*abcdef(1,2)*s+abcdef(3,2)*t+abcdef(4,2)),(2*abcdef(2,2)*t+abcdef(3,2)*s+abcdef(5,2))];

st=[s.^2 t.^2 s.*t s t];
%
%On calcule Xc=X1+J*S aux noeuds de Gauss dont le nombre est Ngauss
  Xc=zeros(Ngauss,2); 
  for i=1:Ngauss
      Xc(i,:)=abcdef(6,:)+st(i,:)*abcdef(1:5,:);
  end
  

     
%on calcule Fe1 Fe2 et Fe2...Fe6 les composantes de Fe
  for i=1:6
      s=0;
      for j=1:Ngauss
      s=s+f(Xc(j,1),Xc(j,2),temps);
      end
      Fe(i)=s;
  end
end



