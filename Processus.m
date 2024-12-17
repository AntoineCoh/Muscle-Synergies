
%%              Extraction et stockage des synergies      
%%

Nbpart=19;
CVaf=cell(16,Nbpart);           % Données stockées dans des Cell Arrays
Wall=cell(16,Nbpart);           % Structure du stockage :
Call=cell(16,Nbpart);           % 1-6 int, 7-12 exoint
Vall=cell(16,Nbpart);           % 13-14 mean, 15-16 std
Mfall=cell(2,Nbpart);
Meall=cell(2,Nbpart);           % Stockage pour vérification
Tall=cell(2,Nbpart);
tbnorm=zeros(Nbpart,8);
load("Mball.mat");

for l=1:2
for p=1:Nbpart
    % if p==5 ||                 Participants à exclure
    %     continue
    % end

for tp=1:2                      % Type d'essai : 1-INT   2-EXOINT

if l==1
    if tp==1
        type="int";
    else
        type="exoint";
    end

addpath("./Participants/P"+p+"/"+type);
load('Donnees.mat')

%%                             Timestamps
%%
addpath("./Participants/P"+p+"/"+type);

[C,D] = butter(6,2/500,'low');
AZf=filtfilt(C,D,AZ);
AXf=filtfilt(C,D,AX);
AZf=-AZf;
wYf=filtfilt(C,D,wY);

if p==2 || p==12 || (p==5 && tp==2)
    timestamps=zeros(1,12);
    ind=1;
    for i=1:length(locsmean)
        if locsmean(i)==1
            timestamps(ind)=i;
            timestamps(ind+1)=i+1000;
            ind=ind+2;
        end
    end

else

[pks,pos] = findpeaks(AZf,"MinPeakDistance",length(AZf)/6-2000);
i=1;
while size(pos)~=6
    if pks(1,i)<6
        pos(:,i)=[];
        pks(:,i)=[];
    else
        i=i+1;
    end
end
locs(1,:)=pos;          % 1e ligne AZ

[pks,pos] = findpeaks(AXf,"MinPeakDistance",length(AXf)/6-2000);
i=1;
while size(pos)~=6
    if pks(1,i)<-6
        pos(:,i)=[];
        pks(:,i)=[];
    else
        i=i+1;
    end
end
locs(2,:)=pos;          % 2e ligne AX

if tp==1 && (p==4 || p==11 || p==15 || p==18)         % Données avec
    for j=1:6                                         % AX corrompu
        for i=locs(1,j)-500:locs(1,j)+500             % ou imprécis
            if sign(wYf(i))~=sign(wYf(i+1))
                locs(3,j)=i;    % 3e ligne VY
            elseif i==length(AZf)-1
                break
            end
        end
    end
    locsmean=round((locs(1,:)+locs(3,:))/2);
elseif tp==1 && (p==10 || p==1)                       % Données avec
    for j=1:6                                         % AZ corrompu
        for i=locs(2,j)-500:locs(2,j)+500             % ou imprécis
            if sign(wYf(i))~=sign(wYf(i+1))
                locs(3,j)=i;
            elseif i==length(AXf)-1
                break
            end
        end
    end
    locsmean=round((locs(2,:)+locs(3,:))/2);
else
    locsmean=round(mean(locs,1));                     % Données correctes
    for j=1:6
        for i=locsmean(1,j)-500:locsmean(1,j)+500
            if sign(wYf(i))~=sign(wYf(i+1))
                locs(3,j)=i;
            elseif i==length(AXf)-1
                break
            end
        end
    end
    locsmean=round(mean(locs,1));
end


for j=1:6
    c=0;
    i=locsmean(1,j)+500;
    while c~=1
        if sign(wYf(i)+.2)~=sign(wYf(i+1)+.2)   % wY redescend sous le
            locs(4,j)=i;                        % seuil de 0.2 rad/s
            c=1;
        elseif i==length(AXf)-1
            break
        end
        i=i+1;
    end
end

timestamps=zeros(1,12);
for i=1:6
    timestamps(1,2*i-1)=locsmean(1,i);
    timestamps(1,2*i)=locs(4,i);
end
if timestamps(1,12)==0
    timestamps(1,12)=length(wYf);
end
end
%                                  Correction manuelle de certains indices

if p==2 && tp==1
    timestamps=timestamps+1000;
end
if p==3 && tp==1
    timestamps(1,12)=23000;
end
if p==4 && tp==2
    timestamps(1,6)=12100; timestamps(1,10)=23000; timestamps(1,12)=27100;
end
if p==6 && tp==2
    timestamps(1,1)=4400;
end
if p==8 && tp==1
    timestamps(1,1)=1800; timestamps(1,3)=6170;
end
if p==8 && tp==2
    timestamps(1,12)=35100;
end
if p==10 && tp==1
    timestamps(1,1)=890;   timestamps(1,2)=2510;  timestamps(1,3)=5970;
    timestamps(1,4)=7650;  timestamps(1,5)=11140; timestamps(1,6)=12590;
    timestamps(1,7)=17870; timestamps(1,8)=19700; timestamps(1,9)=23810;
    timestamps(1,10)=25480;timestamps(1,11)=29590;timestamps(1,12)=31300;
end
if p==11 && tp==1
    timestamps(1,8)=17000;
end
if p==14 && tp==1
    timestamps(1,2)=3100; timestamps(1,4)=7000; timestamps(1,6)=11300;
    timestamps(1,7)=24060; timestamps(1,8)=25600; timestamps(1,12)=33200;
end
if p==14 && tp==2
    timestamps(1,7)=timestamps(1,8); timestamps(1,8)=25800;
    timestamps(1,9)=timestamps(1,10); timestamps(1,10)=29200;
end
if p==15 && tp==2
    timestamps(1,4)=10870;
end
if p==16 && tp==2
    timestamps(1,1)=3300;
end
if p==17 && tp==2
    timestamps(1,1)=3100; timestamps(1,2)=4530; timestamps(1,10)=28200;
end
if p==18 && tp==1
    timestamps(1,4)=8300;
end
if p==19 && tp==1
    timestamps(1,2)=3100;
end
if p==19 && tp==2
    timestamps(1,2)=3000; timestamps(1,6)=10200; timestamps(1,8)=14100;
end


if p==1 && tp==2
   timestamps=timestamps-650;
end
if p==2 && tp==1
    timestamps(1,1)=timestamps(1,1)-700;
    timestamps(1,2)=timestamps(1,2)-700;
    timestamps(1,3)=timestamps(1,3)-300;
    timestamps(1,4)=timestamps(1,4)-300;
    timestamps(1,5)=timestamps(1,5)-100;
    timestamps(1,6)=timestamps(1,6)-100;
    timestamps(1,7)=timestamps(1,7)-400;
    timestamps(1,8)=timestamps(1,8)-400;
    timestamps(1,9)=timestamps(1,9)-400;
    timestamps(1,10)=timestamps(1,10)-400;
    timestamps(1,11)=timestamps(1,11)-600;
    timestamps(1,12)=timestamps(1,12)-600;
end
if p==2 && tp==2
    for i=1:12
        if mod(i,2)==1
            timestamps(1,i)=timestamps(1,i)-400;
        else
            timestamps(1,i)=timestamps(1,i)-200;
        end
    end
    timestamps(1,7)=timestamps(1,7)-100;
    timestamps(1,8)=timestamps(1,8)-100;
end
if p==3 && tp==1
   timestamps=timestamps-50;
end
if p==3 && tp==2
   timestamps=timestamps-400;
   timestamps(1,8)=17000;
end
if p==5 && tp==2
   timestamps=timestamps-400;
end
if p==6 && tp==1
   timestamps=timestamps-500;
end
if p==6 && tp==2
   timestamps=timestamps-300;
end
if p==7 && tp==1
   timestamps=timestamps-500;
end
if p==7 && tp==2
   timestamps=timestamps-600;
end
if p==8 && tp==1
   timestamps=timestamps-400;
end
if p==9 && tp==1
   timestamps=timestamps-100;
end
if p==10 && tp==1
   timestamps=timestamps-200;
   timestamps(1,8)=timestamps(1,8)-300;
end
if p==11 && tp==1
   timestamps=timestamps-100;
end
if p==11 && tp==2
   timestamps=timestamps-300;
   timestamps(1,4)=timestamps(1,4)-100;
end
if p==12 && tp==1
    timestamps(1,5)=timestamps(1,5)+100; timestamps(1,6)=timestamps(1,6)+100;
    timestamps(1,7)=timestamps(1,7)+100; timestamps(1,8)=timestamps(1,8)+100;
    timestamps(1,9)=timestamps(1,9)+200; timestamps(1,10)=timestamps(1,10)+200;
    timestamps(1,11)=timestamps(1,11)+200; timestamps(1,12)=timestamps(1,12)+200;
end
if p==12 && tp==2
    timestamps(1,1)=timestamps(1,1)-100; timestamps(1,2)=timestamps(1,2)-100;
    timestamps(1,3)=timestamps(1,3)-200; timestamps(1,4)=timestamps(1,4)-200;
    timestamps(1,5)=timestamps(1,5)-200; timestamps(1,6)=timestamps(1,6)-200;
    timestamps(1,7)=timestamps(1,7)-200; timestamps(1,8)=timestamps(1,8)-200;
    timestamps(1,9)=timestamps(1,9)-200; timestamps(1,10)=timestamps(1,10)-200;
end
if p==13 && tp==1
   timestamps=timestamps-300;
end
if p==13 && tp==2
   timestamps(1,2)=timestamps(1,2)-100;
end
if p==14 && tp==1
   timestamps=timestamps-200;
end
if p==15 && tp==1
   timestamps=timestamps-700;
end
if p==15 && tp==2
   timestamps=timestamps-400;
end
if p==16 && tp==2
   timestamps=timestamps+100;
end
if p==17 && tp==1
   timestamps(1,2)=timestamps(1,2)-300;
end
if p==18 && tp==2
   timestamps=timestamps+100;
end


if timestamps(1,12)>length(Mball{tp,p})
    timestamps(1,12)=length(Mball{tp,p});
end

Tall{tp,p}=timestamps;
clear locs

%%                          EMG Pre-processing
%%
[B,A] = butter(6, [30,350]/500,'bandpass');
[C,D] = butter(4,15/500,'low');

M=Mball{tp,p};
[nbm,col]=size(M);
Mf=zeros(nbm,col);
Me=zeros(nbm,col);

for i=1:nbm
    Mf(i,:)=filtfilt(B,A,M(i,:));           % Bandpass filtering
    if p==3 && tp==1 && i==2
        for j=1:length(Mf(i,:))
            if Mf(i,j)>.2 || Mf(i,j)<-.2
                Mf(i,j)=.2;
            end
        end
    end
    if p==6 && tp==1 && i==8
        for j=1:length(Mf(i,:))
            if Mf(i,j)>.3 || Mf(i,j)<-.3
                Mf(i,j)=.3;
            end
        end
    end  
    %Mf(i,:)=eemd_it(Mf(i,:));              % Denoising par eemd-it
    Me(i,:)=abs(Mf(i,:));                   % Rectify
    Me(i,:)=filtfilt(C,D,Me(i,:));          % Lowpass
    if tp==1
        maxtemp=zeros(1,6);
        for u=1:6
            if p==6 && u==2
                continue
            else
                Mtemp=Fenetre(Me,Tall{tp,p}(1,2*u-1),Tall{tp,p}(1,2*u));
                maxtemp(1,u)=max(Mtemp(i,:));
            end
        end
        tbnorm(p,i)=max(maxtemp);
        Me(i,:)=abs(Me(i,:))/tbnorm(p,i);     % Enveloppe normalisée
    else
        Me(i,:)=abs(Me(i,:))/tbnorm(p,i);
    end
end

Mfall{tp,p}=Mf;
Meall{tp,p}=Me;
%disp("(p:"+p+",tp:"+tp+")")

tiledlayout(2,1); t1=nexttile; hold on; xline(timestamps,'--'); yline(0)
plot(AXf,'r'); plot(-AZf,'b'); plot(4*wYf,'y','LineWidth',2);
legend(t1,{'','','','','','','','','','','','','','AX','AZ','wY'})
t2=nexttile;
hold on; xline(timestamps,'--')
area(Me(2,:),'FaceColor','g','EdgeColor',"#77AC30","FaceAlpha",.3);
area(Me(4,:),'FaceColor',"#D95319",'EdgeColor',"#A2142F","FaceAlpha",.3)
area(Me(5,:),'FaceColor','m','EdgeColor',"#7E2F8E","FaceAlpha",.3)
legend(t2,{'','','','','','','','','','','','','Trap.Med.','Longiss.','Biceps'})

%pause
close

%%                            Calculs VAF
%%

for i=1:6

Mtemp=Fenetre(Me,Tall{tp,p}(1,2*i-1),Tall{tp,p}(1,2*i));

[nbm,~]=size(Mtemp);
VAFk=zeros(nbm,1);

for ki=1:nbm
[~,~,vaftemp] = NMF(Mtemp,ki);
VAFk(ki,:)=vaftemp;
end

CVaf{i+6*(tp-1),p}=VAFk;
end
                                    % matrice 3D
CVafz=cat(3,CVaf{1+6*(tp-1),p},CVaf{2+6*(tp-1),p},CVaf{3+6*(tp-1),p},...
            CVaf{4+6*(tp-1),p},CVaf{5+6*(tp-1),p},CVaf{6+6*(tp-1),p});
stdCVaf=std(CVafz,0,3);             % std sur dim 3
CVafm=mean(CVafz,3);                % mean sur dim 3
CVaf{12+tp,p}=CVafm;
CVaf{14+tp,p}=stdCVaf;

%%                              Synergies
else
figure
Ctemp=cell(6,2);
Wtemp=cell(6,2);
k=tabk(p,1);

if tp==1
for i=1:6

Mtemp=Fenetre(Meall{tp,p},Tall{tp,p}(1,2*i-1),Tall{tp,p}(1,2*i));
if p==3
    Mtemp=Fenetre(Meall{tp,p},Tall{tp,p}(1,2*i-1),Tall{tp,p}(1,2*i),2);
elseif p==6
    Mtemp=Fenetre(Meall{tp,p},Tall{tp,p}(1,2*i-1),Tall{tp,p}(1,2*i),8);
end
[W,C,vaf] = syn(Mtemp,k);           % Extraction des synergies

if tp==2 && i==1
    Wztemp=cat(3,Wall{1,p},Wall{2,p},Wall{3,p},...
                 Wall{4,p},Wall{5,p},Wall{6,p}); 
    Wtemp=mean(Wztemp,3);
    [W,C]=Reorg_Compar(Wtemp,W,C);
elseif i==2                                     % Réorganisation
    [W,C]=Reorg_Compar(Wall{1+6*(tp-1),p},W,C);   % des vecteurs W
elseif i==3
    Wztemp=cat(3,Wall{1+6*(tp-1),p},Wall{2+6*(tp-1),p});
    Wtemp=mean(Wztemp,3);
    [W,C]=Reorg_Compar(Wtemp,W,C);
elseif i==4
    Wztemp=cat(3,Wall{1+6*(tp-1),p},Wall{2+6*(tp-1),p},Wall{3+6*(tp-1),p});
    Wtemp=mean(Wztemp,3);
    [W,C]=Reorg_Compar(Wtemp,W,C);
elseif i==5
    Wztemp=cat(3,Wall{1+6*(tp-1),p},Wall{2+6*(tp-1),p},...
                 Wall{3+6*(tp-1),p},Wall{4+6*(tp-1),p});
    Wtemp=mean(Wztemp,3);
    [W,C]=Reorg_Compar(Wtemp,W,C);
elseif i==6
    Wztemp=cat(3,Wall{1+6*(tp-1),p},Wall{2+6*(tp-1),p},Wall{3+6*(tp-1),p},...
                 Wall{4+6*(tp-1),p},Wall{5+6*(tp-1),p});
    Wtemp=mean(Wztemp,3);
    [W,C]=Reorg_Compar(Wtemp,W,C);
end
if tp==2 && i~=1
    Wztemp2=cat(3,Wall{1,p},Wall{2,p},Wall{3,p},...
                 Wall{4,p},Wall{5,p},Wall{6,p}); 
    Wtemp2=mean(Wztemp2,3);                       % Vérif. de l'orga.
    [W2,~]=Reorg_Compar(Wtemp2,W,C);             % des W inter-essais
    if W2~=W
        disp("Vérifier W"+i+" (p:"+p+",tp:"+tp+")")
        % if dot((reshape(Wtemp2,[],1)),(reshape(W2,[],1)))...
        %  > dot((reshape(Wtemp2,[],1)),(reshape(W2,[],1)))
        % % if mean(dot(Wtemp2,W2,1))>mean(dot(Wtemp,W,1))
        %     W=W2;
        %     disp("Switch W"+i+" (p:"+p+",tp:"+tp+")")
        % end
    end
end
Wall{i+6*(tp-1),p}=W;                           % Stockage dans
Call{i+6*(tp-1),p}=C;                           % les Cell Arrays
Vall{i+6*(tp-1),p}=vaf;
xq=linspace(1,length(C),100);                   % Downsampling de
Ctemp{i,2}=(interp1(C',xq))';                   % C à 100 points

% cms = colormap(cool(6)); %essayer cool, summer ou copper
% tl=tiledlayout(3,2*k);
% for e=1:6
%     W=Wall{e+(tp-1)*6,p};
%     for j=1:k
%         nexttile;
%         b=bar((W(:,j)/norm(W(:,j))),'b'); b.FaceColor = cms(7-i,:);
%         nt=gca;
%         nt.YLim=[0;1];
%     end
% end
%sgtitle('Synergies extraites des 6 essais du participant 4 sans exosquelette')
% tl.Padding = 'compact'; tl.TileSpacing = 'compact';
% cb = colorbar; cb.Layout.Tile = 'east';
% cb.Label.String = 'n° mouvement'; cb.Label.FontSize = 11;
% ech=linspace(0,1,7)-0.08333333; ech(:,1)=[];
% cb.Ticks=ech; cb.TickLabels=6:-1:1;

end

elseif tp==2

tabsim=zeros(1,6);
for i=1:6
Mtemp=Fenetre(Meall{tp,p},Tall{tp,p}(1,2*i-1),Tall{tp,p}(1,2*i));

if p==3
    Mtemp=Fenetre(Meall{tp,p},Tall{tp,p}(1,2*i-1),Tall{tp,p}(1,2*i),2);
elseif p==6
    Mtemp=Fenetre(Meall{tp,p},Tall{tp,p}(1,2*i-1),Tall{tp,p}(1,2*i),8); 
end

[W,C,vaf] = syn(Mtemp,k);
Vall{i+6*(tp-1),p}=vaf;
Wztemp=cat(3,Wall{1,p},Wall{2,p},Wall{3,p},...
             Wall{4,p},Wall{5,p},Wall{6,p}); 
Wcomp=mean(Wztemp,3);
[W,C]=Reorg_Compar(Wcomp,W,C);
tabsim(1,i)=mean(dot(Wcomp,W,1));
Wtemp{i,1}=W;
Ctemp{i,1}=C;
end

temp=sort(tabsim,'descend');

for i=1:6
    index=find(tabsim==temp(1,i));
    W=Wtemp{index,1};
    C=Ctemp{index,1};
    if i==1
        Wcomp=W;
        W1=W;
    elseif i==2
        [W,C]=Reorg_Compar(Wcomp,W,C);
        W2=W;
        Wcompz=cat(3,W1,W2);
        Wcomp=mean(Wcompz,3);
    elseif i==3
        [W,C]=Reorg_Compar(Wcomp,W,C);
        W3=W;
        Wcompz=cat(3,W1,W2,W3);
        Wcomp=mean(Wcompz,3);
    elseif i==4
        [W,C]=Reorg_Compar(Wcomp,W,C);
        W4=W;
        Wcompz=cat(3,W1,W2,W3,W4);
        Wcomp=mean(Wcompz,3);
    elseif i==5
        [W,C]=Reorg_Compar(Wcomp,W,C);
        W5=W;
        Wcompz=cat(3,W1,W2,W3,W4,W5);
        Wcomp=mean(Wcompz,3);
    elseif i==6
        [W,C]=Reorg_Compar(Wcomp,W,C);
    end
    Wall{index+6*(tp-1),p}=W;
    Call{index+6*(tp-1),p}=C;
    xq=linspace(1,length(C),100);
    Ctemp{index,2}=(interp1(C',xq))';
end

end
                                % matrice 3D
Wz=cat(3,Wall{1+6*(tp-1),p},Wall{2+6*(tp-1),p},Wall{3+6*(tp-1),p},...
         Wall{4+6*(tp-1),p},Wall{5+6*(tp-1),p},Wall{6+6*(tp-1),p}); 
stdW=std(Wz,0,3);               % std sur dim 3
Wm=mean(Wz,3);                  % mean sur dim 3
Wall{12+tp,p}=Wm;
Wall{14+tp,p}=stdW;
                                % matrice 3D
Vz=cat(3,Vall{1+6*(tp-1),p},Vall{2+6*(tp-1),p},Vall{3+6*(tp-1),p},...
         Vall{4+6*(tp-1),p},Vall{5+6*(tp-1),p},Vall{6+6*(tp-1),p}); 
stdV=std(Vz,0,3);               % std sur dim 3
Vm=mean(Vz,3);                  % mean sur dim 3
Vall{12+tp,p}=Vm;
Vall{14+tp,p}=stdV;

                                % matrice 3D
Cz=cat(3,Ctemp{1,2},Ctemp{2,2},Ctemp{3,2},...
         Ctemp{4,2},Ctemp{5,2},Ctemp{6,2});
stdC=std(Cz,0,3);               % std sur dim 3
Cm=mean(Cz,3);                  % mean sur dim 3
Call{12+tp,p}=Cm;
Call{14+tp,p}=stdC;

                                % Indicateurs
disp("std moyen : "+mean(stdW,'all')+" (p:"+p+",tp:"+tp+")")
%pause
close
end
end
end
%%                 Plot Courbes VAF et choix de K
%%
if l==1

tl=tiledlayout(1,2);
cm = colormap(turbo(Nbpart));
for i=1:2
    nexttile
    for x=1:Nbpart
    plot(CVaf{12+i,x},'Color',cm(x,:),'LineWidth',2)
    hold on
    end
    ylim([80;100])
    yline(90,'--'); xline(2,'--'); xline(3,'--');
    hold off               
    ylabel('VAF (%)') 
    xlabel('Nombre de synergies')
    if i==1                             % Chaque courbe affichée
        title('Sans exosquelette')      % représente la moyenne des
    else                                % courbes VAF obtenues sur les
        title('Avec exosquelette')      % 6 essais pour chaque participant
    end
end
tl.Padding = 'compact'; tl.TileSpacing = 'compact';
cb = colorbar; cb.Layout.Tile = 'east';
cb.Label.String = 'n° participant';
ech=linspace(0,1,20)-0.02631579; ech(:,1)=[];
cb.Ticks=ech; cb.TickLabels=1:Nbpart;
title(tl,'VAF en fonction du nombre de synergies')
close

% prompt="Valeur de k : ";
% k=input(prompt);                        % Sélection de K
%k=3;
load('tabk.mat')

end
%%                              Sauvegarde
%%
save Wall.mat Wall
save Call.mat Call
save Vall.mat Vall
save CVaf.mat CVaf
save Tall.mat Tall
save Mfall.mat Mfall
save Meall.mat Meall
end
clearvars -except Wall Call Vall CVaf Tall Mfall Meall tbnorm