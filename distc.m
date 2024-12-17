function [distance] = distc(x1,x2)

% Calcule la distance entre deux profils d'activation x1 et x2.

k=length(x1)/100;
tab=zeros(1,k);

for i=1:k
    %tab(1,i)=max(xcorr(x1(:,1+(i-1)*100:100*i),x2(:,1+(i-1)*100:100*i),'normalized'));
    tab(1,i)=min(corrcoef(x1(:,1+(i-1)*100:100*i),x2(:,1+(i-1)*100:100*i)),[],'all');
end

distance=1-mean(tab);

end