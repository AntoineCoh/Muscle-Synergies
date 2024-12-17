function [distance] = dist(x1,x2)
% Calcule la distance entre deux sets de k synergies x1 et x2.
% Cette distance est définie comme étant la moyenne des distances
% cosinus entre les k synergies alignées des deux sets.a

if mod(length(x1),3)==0
    k=3;
elseif length(x1)<=8
    k=1;
else
    k=2;
end
nbm=length(x1)/k;

tab=zeros(1,k);

for i=1:k
    tab(1,i)=dot(x1(:,1+(i-1)*nbm:nbm*i),x2(:,1+(i-1)*nbm:nbm*i));
end

distance=1-mean(tab);

end