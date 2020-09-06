function [a,v,d] = HHTmethod(K,M,f,NT,h,alphaf2,alpham2,Beta2,Gamma2,elem,k,X,T,b,Ngauss,D,pho,H,alpha,Beta,Gamma)

d=zeros(size(K,1),NT);
v=zeros(size(K,1),NT);
a=zeros(size(K,1),NT);
%--------Conditions aux bords-----------
%On met "0" dans les lignes et les colonnes qui correspondent aux noeuds
%bords dans K et F
K(b,:)=0; K(:,b)=0;
M(b,:)=0; M(:,b)=0;
%On met "1" sur la diagonale de K correspodant aux aux noeuds bords 
%K(b,b)=speye(length(b),length(b)); % cas où K est sparse
K(b,b)=eye(length(b),length(b));
M(b,b)=eye(length(b),length(b));
% K et F après conditions aux bords
%
temps=0;
F = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps);
a(:,1)=(1/(pho*H))*M\F; % Initialisation de l'accélération à t=0
K=D*(K-M)*(M\K)+D*K;
M=pho*H*M;


%Genenralized alpha scheme (doesn't work)
%for i=2:NT
%temps=(i-1)*h*(1-alphaf2)+alphaf2*(i-2)*h;
%F = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps);
%a(:,i)=((1-alpham2)*M+(h^2)*Beta2*(1-alphaf2)*K)\(F-alpham2*M*a(:,i-1)-K*alphaf2*d(:,i-1)-K*(1-alphaf2)*(d(:,i-1)+h*v(:,i-1)+(h^2)*(1-Beta2)*a(:,i-1)));
%v(:,i)=v(:,i-1)+h*(Gamma2*a(:,i)+(1-Gamma2)*a(:,i-1));
%d(:,i)=d(:,i-1)+h*v(:,i-1)+(h^2)*(Beta2*a(:,i)+(0.5-Beta2)*a(:,i-1));
%end

%----------------%% HHT Method %%---------------------%
for i=2:NT
temps1=(i-2)*h;
F1 = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps1);
temps2=(i-1)*h;
F2 = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps2);
a(:,i)=(M+(h^2)*Beta*(1-alpha)*K)\((1-alpha)*F2+alpha*F1-K*d(:,i-1)-h*(1-alpha)*K*v(:,i-1)-(h^2)*(1-alpha)*(0.5-Beta)*K*a(:,i-1));
a(:,i)=(M+(h^2)*Beta*(1-alpha)*K)\(F1-K*d(:,i-1)-h*(1-alpha)*K*v(:,i-1)-(h^2)*(1-alpha)*(0.5-Beta)*K*a(:,i-1));
v(:,i)=v(:,i-1)+h*(Gamma*a(:,i)+(1-Gamma)*a(:,i-1));
d(:,i)=d(:,i-1)+h*v(:,i-1)+(h^2)*(Beta*a(:,i)+(0.5-Beta)*a(:,i-1));
end
d=0.5*d;
%----------------%% HHT Method %%---------------------%
%for i=2:NT
%temps1=(i-2)*h;
%F1 = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps1);
%temps2=(i-1)*h;
%F2 = AssemblageSM(elem,k,X,T,b,f,Ngauss,h,pho,H,temps2);
%d(:,i)=2*(eye(size(d,1),size(d,1))+(h^2)*Beta*(1-alpha)*(M\K))\(d(:,i-1)+h*v(i-1)+(h^2)*(0.5-Beta)*a(:,i-1)+(h^2)*Beta*(1-alpha)*(M\F2)+(h^2)*Beta*alpha*(M\F1)-(h^2)*Beta*alpha*(M\(K*d(:,i-1))));
%a(:,i)=M\(F2-K*d(:,i));
%v(:,i)=v(:,i-1)+h*(Gamma*a(:,i)+(1-Gamma)*a(:,i-1));
%end
end