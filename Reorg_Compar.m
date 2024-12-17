function [Wo2,Co2]=Reorg_Compar(W1,W2,C2)
% Réorganise le vecteur de synergies W2 par rapport à W1
% Le vecteur de coefficients C2 est réorganisé en conséquences

[d,k]=size(W1);
sim=zeros(k,k); % Matrice carrée contenant les coef. de similarité

% Normalisation des synergies
Wn1=zeros(d,k);
Wn2=zeros(d,k);
for i=1:k
    Wn1(:,i)=W1(:,i)/norm(W1(:,i));
    Wn2(:,i)=W2(:,i)/norm(W2(:,i));
end

% Produit scalaire des vecteurs normalisés
for i=1:k
    for j=1:k                            
        sim(i,j)=dot(Wn1(:,i),Wn2(:,j));
    end
end

%% Synergies de W2 ordonnées selon W1

Wo2=zeros(d,k);
Co2=zeros(k,length(C2));

for i=1:k
    coefmax=max(sim,[],'all');
    [row,col]=find(sim==coefmax);
    Wo2(:,row)=W2(:,col);
    Co2(row,:)=C2(col,:);
    sim(:,col)=0; sim(row,:)=0;
end

%% Synergies ordonnées par similarité décroissante

% Wo1=zeros(d,k);
% Wo2=zeros(d,k);
% 
% for i=1:k
%     coefmax=max(sim,[],'all');
%     [row,col]=find(sim==coefmax);
%     Wo2(:,i)=W2(:,col);
%     Wo1(:,i)=W1(:,row);
%     sim(:,col)=0; sim(row,:)=0;
% end
end