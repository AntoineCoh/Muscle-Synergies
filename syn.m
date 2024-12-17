function [W,C,vaf] = syn(M,k)

maxiter=100; % A MODIFIER
[d,n]=size(M);

tabvaf=zeros(1,maxiter);
tabw=zeros(d,k*maxiter);
tabc=zeros(k*maxiter,length(M));

for i=1:maxiter
    [W,C,vaf] = NMF(M,k);
    tabvaf(1,i)=vaf;
    tabw(1:d,k*(i-1)+1:k*i)=W;
    tabc(k*(i-1)+1:k*i,1:n)=C;
end

[vafmax,ind]=max(tabvaf);
vaf=vafmax;
W=tabw(1:d,k*(ind-1)+1:k*ind);
C=tabc(k*(ind-1)+1:k*ind,1:n);

Wn=zeros(d,k);
for i=1:k
    Wn(:,i)=W(:,i)/norm(W(:,i)); % Normalisation
end
W=Wn;

end