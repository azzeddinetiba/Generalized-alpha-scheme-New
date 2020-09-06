function ERR = H1_Error( U,uee,X,T,Ngauss,elem,k,error_L2 )
%H1_Error calcule l'errer dans l'espace de Sobolev H1

[Xgauss,Wgauss]=Quadrature(elem, Ngauss);
[~,grad,~,~] = FoncChap(k,Xgauss);

ERR=0;
for j=1:Ngauss
    interm=0;
for ie=1:size(T,1)
    
x=X(T(ie,:),:); %Les noeuds réels de l'élement °ie
s = Xgauss(:,1); 
t = Xgauss(:,2); 

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


 for i=1:size(T,2)
     interm=interm+(grad([j j+Ngauss],i)*(U(T(ie,i))-uee(T(ie,i))))'*(grad([j j+Ngauss],i)*(U(T(ie,i))-uee(T(ie,i))));
 end
end
    ERR=ERR+interm*abs(det(J([j j+Ngauss],:)))*Wgauss(j);
end

ERR=sqrt(ERR+(error_L2)^2);
end

