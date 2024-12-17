
function [siltab] = Silhou (X1,X2)
% Coefficient de silhouette d'un partitionnement 
% de données contenant deux clusters X1 et X2

% Remplacer la fonction 'dist' par 'distc' pour les profils d'activation

[n1,~]=size(X1);    % Taille du cluster 1
[n2,~]=size(X2);    % Taille du cluster 2

A1=zeros(n1-1,n1);
B1=zeros(n2,n1);
A2=zeros(n2-1,n2);
B2=zeros(n1,n2);
sil=zeros(1,n1+n2);

for i=1:n1
    jind=1;
    for j=1:n1
        if j==i
            continue
        else
            A1(jind,i)=dist(X1(i,:),X1(j,:));
            jind=jind+1;
        end
    end
    for j=1:n2
        B1(j,i)=dist(X1(i,:),X2(j,:));
    end                                % Silhouette score des points de X1
    sil(1,i)=((mean(B1(:,i))-mean(A1(:,i)))/max(mean(B1(:,i)),mean(A1(:,i))));
end

for i=1:n2
    jind=1;
    for j=1:n2
        if j==i
            continue
        else
            A2(jind,i)=dist(X2(i,:),X2(j,:));
            jind=jind+1;
        end
    end
    for j=1:n1
        B2(j,i)=dist(X2(i,:),X1(j,:));
    end                                % Silhouette score des points de X2
    sil(1,i+n1)=((mean(B2(:,i))-mean(A2(:,i)))/max(mean(B2(:,i)),mean(A2(:,i))));
end

siltab=[mean(sil) median(sil) min(sil) max(sil) std(sil) 1-mean(A1,"all") 1-mean(A2,"all")];
                                                       % indicateurs de
                                                       % similarité
                                                       % intra-condition
end