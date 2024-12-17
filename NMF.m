function [W,C,vaf] = NMF(M,k)

[d,n]=size(M);
[W]=rand(d,k);              % Initialisation
[C]=rand(k,n);              % Valeurs aléatoires entre 0 et 1
it=1;                       % Indice d'itération
tabvaf=[1;200];             % Tableau de stockage des VAF
condition=true;

while condition==1

    t1=(W')*M;
    t2=(W')*W*C;
    C=C.*t1./t2;            % Mise à jour de C
    
    t1=M*(C');
    t2=W*C*(C');
    W=W.*t1./t2;            % Mise à jour de W

    [E]=M-W*C;              % Erreur de reconstruction
    Ef=sumsqr(E);
    Mf=sumsqr(M);
    tabvaf(:,it)=(1-Ef/Mf)*100; % Calcul de VAF

    % 5 itérations consécutives sans augmentation de 0,01% de VAF
    if it>5 & tabvaf(:,it)-tabvaf(:,it-4)<0.01 
        condition=0;        % Fin de boucle
        vaf=tabvaf(1,it);
    end
    it=it+1;                % itération suivante
end
end