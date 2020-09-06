function [X,T,b] = MAILLAGENONSTRUCT(fd, fh, h0, bbox, p_fix)
%MAILLAGENONCONSTRUCT crée un maillage non structuré basé sur distmesh
% distmesh un code de maillage par G.Strang
% La partie suivante est la transformation de distmesh vers un maillage P2
% et création du vecteur b des noeuds de bord

[X,T,~,b]=distmesh( fd, fh, h0, bbox, p_fix);
T=[T(:,1) zeros(size(T,1),1) T(:,2) zeros(size(T,1),1) T(:,3) zeros(size(T,1),1)];
r=0;
for i=1:size(T,1)
    for j=1:size(X,1)
    if  X(j,:) == 0.5*(X(T(i,1),:)+X(T(i,3),:))
    r=1;
    T(i,2)=j;
    end
    end
    
    if r==0 
    X(size(X,1)+1,:)=0.5*(X(T(i,1),:)+X(T(i,3),:));
    T(i,2)=size(X,1);
    end
    
    
    
    r=0;
        
    
    for j=1:size(X,1)
    if X(j,:) == 0.5*(X(T(i,3),:)+X(T(i,5),:))
    r=1;
    T(i,4)=j;
    end
    end
    
    if r==0
    X(size(X,1)+1,:)=0.5*(X(T(i,3),:)+X(T(i,5),:));
    T(i,4)=size(X,1);
    end
        
    r=0;
    
    for j=1:size(X,1)
    if  X(j,:) == 0.5*(X(T(i,5),:)+X(T(i,1),:))
    r=1;
    T(i,6)=j;
    end
    end
    
    if r==0
    X(size(X,1)+1,:)=0.5*(X(T(i,5),:)+X(T(i,1),:));
    T(i,6)=size(X,1);
    end
    
    r=0;
   
    
end

 indice=1;
    b=zeros(1,size(b,2));
    for i=1:size(X,1)
    if X(i,1)==1 || X(i,2)==1 || X(i,1)==0 || X(i,2)==0
        b(:,indice)=i; 
        indice=indice+1;
    end
    end
