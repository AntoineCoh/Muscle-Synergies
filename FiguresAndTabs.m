%% Run sections

%%                  k synergies W des 6 essais
load('Wall.mat')
for p=6           % N° Participant
for tp=1:2         % Type : 1-int  2-exoint
figure
[~,k]=size(Wall{1+(tp-1)*6,p});
cms = colormap(cool(6));
tl=tiledlayout(3,2*k);
for i=1:6
    W=Wall{i+(tp-1)*6,p};
    for j=1:k
        if i==1 && j==1 && tp==1 && p==17
            ax1 = nexttile;
        else
        nexttile;
        end
        b=bar((W(:,j)/norm(W(:,j))),'b'); b.FaceColor = cms(7-i,:);
        nt=gca;
        nt.YLim=[0;1];
        set(gca, 'XTick', [1 2 3 4 5 6 7 8])
        xticklabels({'TT','TA','LD','LO','BI','DM','IC','MU'})
    end
end
%sgtitle('')
tl.Padding = 'compact'; tl.TileSpacing = 'compact';
cb = colorbar; cb.Layout.Tile = 'east';
cb.Label.String = '# trial'; cb.Label.FontSize = 11;
ech=linspace(0,1,7)-0.08333333; ech(:,1)=[];
cb.Ticks=ech; cb.TickLabels=6:-1:1;
end
end
%%                          Courbes VAF

load('CVafi.mat')
figure
Nbpart=19;
tl=tiledlayout(1,2);
cm = colormap(turbo(Nbpart-2)); % turbo, jet ou hsv
for i=1:2
    nexttile
    ind=1;
    for x=1:Nbpart
        if x==5 || x==12
            continue
        else
        plot(CVafi{12+i,ind},'Color',cm(ind,:),'LineWidth',1.7)
        hold on
        ind=ind+1;
        end
    end
    xlim([1;8])
    ylim([80;100])
    yline(90,'--'); xline(2,'--'); xline(3,'--');
    hold off               
    ylabel('VAF (%)') 
    xlabel('Number of synergies')
    if i==1
        text(4.2,85,'Without exoskeleton (NAT)')
    else
        text(4.4,85,'With exoskeleton (EXO)')
    end
end
tl.Padding = 'compact'; tl.TileSpacing = 'compact';
cb = colorbar; cb.Layout.Tile = 'east';
cb.Label.String = '# participant';
ech=linspace(0,1,18)-0.0294117647; ech(:,1)=[];
cb.Ticks=ech; cb.TickLabels=[1 2 3 4 6 7 8 9 10 11 13 14 15 16 17 18 19];
fontsize("increase")
fontsize("increase")
%%                   Enveloppes EMG et timestamps
p=19;           % N° Participant
tp=1;           % Type : 1-int  2-exoint
for m=1:2
figure
load('Tall.mat')
load('Meall.mat')
if tp==1
    type="int";
else
    type="exoint";
end
addpath("./Participants/P"+p+"/"+type);
load('Donnees.mat')
[C,D] = butter(6,2/500,'low');
AZf=filtfilt(C,D,AZ);
AXf=filtfilt(C,D,AX);
wYf=filtfilt(C,D,wY);
tl=tiledlayout(5,1);
if m==1
t1=nexttile; hold on; xline(Tall{tp,p},'--','LineWidth',1);
plot(AXf,'color','#F05576','LineWidth',1.1); plot(AZf,'color','#7F55F0','LineWidth',1.1);
ylabel('Acc. (m/s²)'); yyaxis right; plot(wYf,'color','#F1D41E','LineWidth',1.5); yline(0);
ax = gca; ax.YColor = 'k'; ylabel('Vit. ang. (rad/s)');
legend(t1,{'','','','','','','','','','','','','AX','AZ','wY',''})
t2=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(1,:),'FaceColor','r','EdgeColor',"#891c05","FaceAlpha",.4);
legend(t2,{'','','','','','','','','','','','','Trapèze Med.'})
t3=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(2,:),'FaceColor','g','EdgeColor',"#77AC30","FaceAlpha",.4);
legend(t3,{'','','','','','','','','','','','','Trapèze Asc.'})
t4=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(3,:),'FaceColor','b','EdgeColor',"#110589","FaceAlpha",.4);
legend(t4,{'','','','','','','','','','','','','Grd. Dorsal'})
t5=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(4,:),'FaceColor',"#e79a31",'EdgeColor',"#c57508","FaceAlpha",.4);
legend(t5,{'','','','','','','','','','','','','Longissimus'})
else
t6=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(5,:),'FaceColor','m','EdgeColor',"#7E2F8E","FaceAlpha",.4);
legend(t6,{'','','','','','','','','','','','','Biceps'})
t7=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(6,:),'FaceColor','y','EdgeColor',"#858905","FaceAlpha",.4);
legend(t7,{'','','','','','','','','','','','','Deltoïde'})
t8=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(7,:),'FaceColor',"c",'EdgeColor',"#058183","FaceAlpha",.4)
legend(t8,{'','','','','','','','','','','','','Iliocostalis'})
t9=nexttile;
hold on; xline(Tall{tp,p},'--','LineWidth',1)
area(Meall{tp,p}(8,:),'FaceColor','#c883e3','EdgeColor',"#9841bb","FaceAlpha",.4)
legend(t9,{'','','','','','','','','','','','','Multifide'})
end
xlabel('Temps (ms)');
tl.Padding = 'compact'; tl.TileSpacing = 'compact';
end
hold off
%%              Synergies et coef NAT EXO combinés
p=19;
figure
load('Walla.mat')
load('Calla.mat')

Wn=Walla{13,p};
We=Walla{14,p};
stdn=Walla{15,p};
stde=Walla{16,p};
[~,k]=size(Wn);
tl=tiledlayout(2,k);

if p==3 || p==6
    nbm=7;
else
    nbm=8;
end
W=zeros(nbm,2);
for i=1:k
    nexttile; hold on
    W=[Wn(:,i) We(:,i)];
    b=bar(W,1);
    b(1).FaceColor = "#63CA63";
    b(2).FaceColor = "#8966C6";
    set(gca, 'XTick', [1 2 3 4 5 6 7 8])
    xticklabels({'TT','TA','LD','LO','BI','DM','IC','MU'})
  
    stdW=[stdn(:,i) stde(:,i)];
    x=zeros(2,nbm);
    for j=1:2
        x(j,:)=b(j).XEndPoints;
    end
    e=errorbar(x',W,stdW,'.','Color','k');
    ylim([0;1]);
    hold off
    if i==1
        ylabel("Mean Synergies")
    end
end
legend('NAT','EXO','Standard deviation')

Cn=Calla{13,p};
Ce=Calla{14,p};
stdn=Calla{15,p};
stde=Calla{16,p};
Cpsn=Cn+stdn;
Cmsn=Cn-stdn;
Cpse=Ce+stde;
Cmse=Ce-stde;
[k,n]=size(Cn);

for i=1:k
    nexttile; hold on
    plot(Cn(i,:), 'Color','#2E802E', 'LineWidth',2)
    f=fill([1:1:n n:-1:1],[Cpsn(i,:) fliplr((Cmsn(i,:)))],'c');
    f.FaceColor="#79C679";
    f.FaceAlpha=0.3;
    plot(Ce(i,:), 'Color','#672E80', 'LineWidth',2)
    fe=fill([1:1:n n:-1:1],[Cpse(i,:) fliplr((Cmse(i,:)))],'c');
    fe.FaceColor="#C2AFE4";
    fe.FaceAlpha=0.3;
    ylim([0;max(max(Cpsn,[],'all'),max(Cpse,[],'all'))]);
    if i==1
        ylabel("Mean Activation Profiles")
    end
    if i==2
        xlabel("Movement completion (%)")
    end
end
legend('NAT','NAT std. deviation','EXO','EXO std. deviation')
tl.Padding = 'compact'; tl.TileSpacing = 'compact';

%%                           Dendrogramme
ind=1;
for i=1:19
    if i==3 || i==5 || i==6 || i==12
        continue
    else
    for tp=1:2
        if tp==1
            trial="N";
        else
            trial="E";
        end
        [~,k]=size(Wall{1,i});
        for j=1:k
            labels(ind,1)=(num2str(i)+"."+num2str(j)+"."+trial);
            ind=ind+1;
        end
    end
    end
end
ind=1;
for p=1:19
    if p==3 || p==5 || p==6 || p==12
        continue
    else
    [~,k]=size(Walla{12+tp,p});
    for tp=1:2
        for s=1:k
            X(ind,:)=(Walla{12+tp,p}(:,s))';
            ind=ind+1;
        end
    end
    end
end
figure
Y = pdist(X,'cosine');
SY = squareform(Y);
Z = linkage(Y,'average');
c = cophenet(Z,Y);
I = inconsistent(Z);
N=length(Z);
d=dendrogram(Z,0,'Labels',labels,'ColorThreshold',.15);
ylabel ("Cosine distance")
xlabel ("Synergies")
yline(0.15,'--')
eva=evalclusters(X,'linkage','silhouette','KList',1:5,'Distance','cosine');
linesColor = cell2mat(get(d,'Color'));
colorList = unique(linesColor, 'rows');
for x=1:N
    if linesColor(x,:)==colorList(2,:)
        d(x).Color=[.15 .85 .4];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(3,:)
        d(x).Color=[1 0.35 0.35];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(4,:)
        d(x).Color=[1	0.9	0];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(5,:)
        d(x).Color=[.8 .35 .87];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(6,:)
        d(x).Color=[0.18 0.88 1];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(7,:)
        d(x).Color=[.72 0.63 0.61];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(8,:)
        d(x).Color=[1 0.6 0.4];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(9,:)
        d(x).Color=[0.62 0.55 1];
        d(x).LineWidth=1.5;
    elseif linesColor(x,:)==colorList(10,:)
        d(x).Color=[1 0.22 0.88];
        d(x).LineWidth=1.5;
    end
end
%%                          Clusters
%ligne=participant, colonne=id synergie, chiffre=cluster
tabclust=[3 2 7 3 2 7;...
          3 9 0 4 9 0;...
          0 0 0 0 0 0;...     
          5 2 1 5 2 1;...
          0 0 0 0 0 0;...     % 5
          0 0 0 0 0 0;...
          8 5 1 9 3 1;...
          1 3 0 1 3 0;...
          3 2 0 3 2 0;...
          6 4 0 7 4 0;...     % 10
          1 5 0 1 5 0;...
          0 0 0 0 0 0;...
          3 6 4 5 6 4;...
          1 2 9 1 2 9;...
          1 2 0 6 4 0;...     % 15
          6 2 1 6 2 1;...
          9 7 3 4 7 5;...
          3 2 0 3 2 0;...
          4 8 0 4 8 0];

nbclust=13;
nbpart=15;
figure
nbcol=nbpart*2+1;
tl=tiledlayout(nbclust+1,nbcol);
cms=[0.62 0.55 1;1 0.6 0.4;0.18 0.88 1;1	0.9	0;...
    1 0.22 0.88;.15 .85 .4;1 0.35 0.35;.8 .35 .87];
nexttile
    ax=gca;
    ax.Visible=0;
for i=1:19
    if i==5 || i==12 || i==3 || i==6
        continue
    else
    nexttile([1,2])
    text(0.37,0.3,"P"+num2str(i))
    ax=gca;
    ax.Visible=0;
    end
end
for j=2:nbclust+1
    for i=1:nbcol
        nexttile((j-1)*nbcol+i)
        if i==1 && j<=nbclust-4
            text(0.1,0.47,"C"+num2str(j-1))
        end
        ax=gca;
        ax.Visible=0;
    end
end
ind=1;
for i=1:19
    if i==3 || i==5 || i==6 || i==12
        continue
    else
        for j=1:6
            if tabclust(i,j)==0
                continue
            else
                if j>3
                    tp=2;
                else
                    tp=1;
                end
                nexttile(tabclust(i,j)*nbcol+ind*2+tp-1)
                if tp==1
                    W=Walla{12+tp,i}(:,j);
                    stdW=Walla{14+tp,i}(:,j);
                else
                    W=Walla{12+tp,i}(:,j-3);
                    stdW=Walla{14+tp,i}(:,j-3);
                end
                W=W/norm(W);
                b=bar(W,.77);
                if tabclust(i,j)==9
                    b.FaceColor="#b7a09b";
                else
                    b.FaceColor=cms(tabclust(i,j),:);
                end
                ylim([0;1])
                set(gca, 'XTick', [])
                set(gca, 'YTick', [])
                if tp==1
                    set(gca,'Color',"#afe3af")
                else
                    set(gca,'Color',"#d6c7f0")
                end
            end
        end
        ind=ind+1;
    end
end
tl.Padding ='none'; tl.TileSpacing = 'none';
%%                      Silhouette par synergie
load('tabk.mat')
silh=zeros(19,3);
for p=1:19
    if p==3 || p==6
        nbm=7;
    else
        nbm=8;
    end
    k=tabk(p,1);
    for j=1:k
        X1=zeros(6,nbm);
        X2=zeros(6,nbm);
        for i=1:6
            temp=reshape(Wall{i,p},[],1);
            X1(i,:)=temp(1+(j-1)*nbm:j*nbm,1)';
            temp=reshape(Wall{i+6,p},[],1);
            X2(i,:)=temp(1+(j-1)*nbm:j*nbm,1)';
        end
     if p==3
        X1(5,:)=[]; X1(5,:)=[];
    elseif p==4
        X1(4,:)=[]; X2(2,:)=[];
    elseif p==6
        X1(2,:)=[]; X2(3,:)=[];
    elseif p==7
        X1(1,:)=[]; X1(1,:)=[];
    elseif p==13
        X1(5,:)=[];
    elseif p==14
        X2(6,:)=[];
    elseif p==15
        X1(6,:)=[];
    elseif p==18
        X2(3,:)=[];
    elseif p==17
        X1(5,:)=[];
    elseif p==19
        X1(2,:)=[];
     end
        temp=Silhou(X1,X2);
        silh(p,j)=temp(1,1);
    end
end
%%      Silhouette score par activation

load('tabk.mat')
silh=zeros(19,3);
for p=1:19
    k=tabk(p,1);
    X1=zeros(6,100);
    X2=zeros(6,100);

    for ki=1:3
        if ki==3 && k==2
            continue
        end
        for i=1:6
            xq=linspace(1,length(Call{i,p}(ki,:)),100);
            Ctemp=(interp1((Call{i,p}(ki,:))',xq))';
            X1(i,:)=reshape(Ctemp',[],1)';
            xq=linspace(1,length(Call{i+6,p}(ki,:)),100);
            Ctemp=(interp1((Call{i+6,p}(ki,:))',xq))';
            X2(i,:)=reshape(Ctemp',[],1)';
        end

    if p==3
        X1(5,:)=[]; X1(5,:)=[];
    elseif p==4
        X1(4,:)=[]; X2(2,:)=[];
    elseif p==6
        X1(2,:)=[]; X2(3,:)=[];
    elseif p==7
        X1(1,:)=[]; X1(1,:)=[];
    elseif p==13
        X1(5,:)=[];
    elseif p==14
        X2(6,:)=[];
    elseif p==15
        X1(6,:)=[];
    elseif p==18
        X2(3,:)=[];
    elseif p==17
        X1(5,:)=[];
    elseif p==19
        X1(2,:)=[];
    end
        temp=Silhou(X1,X2);
        silh(p,ki)=temp(1,1);
    end
end
%%          T Tests for each muscle weight (6*2) of each synergy

ttab=zeros(19,24);
for p=1:19
    k=tabk(p,1);
    if p==3 || p==6
        nbm=7;
    else
        nbm=8;
    end
    for j=1:k
        for i=1:nbm
            clear ntab etab
            ind=1;
            for l=1:6
                if isempty(Walla{l,p})
                    continue
                else
                ntab(ind,1)=Walla{l,p}(i,j);
                ind=ind+1;
                end
            end
            ind=1;
            for l=1:6
                if isempty(Walla{l+6,p})
                    continue
                else
                etab(ind,1)=Walla{l+6,p}(i,j);
                ind=ind+1;
                end
            end
            [pvalue,~]=ranksum(ntab,etab); % Wilcoxon rank sum / Mann–Whitney test
            ttab(p,i+(j-1)*nbm)=pvalue;
        end
    end
end
%%                             Amplitude C

tabkc=zeros(19,6);
for p=1:19
    Cint=Calla{13,p};
    Cexo=Calla{14,p};
    [k,~]=size(Cint);
    for i=1:k
        [~,lag]=max(xcorr(Cint(i,:),Cexo(i,:)));
        lag=lag-100;
        if lag<=0
            tempi=Cint(i,1:100+lag);
            tempe=Cexo(i,-lag+1:100);
        else
            tempi=Cint(i,lag+1:100);
            tempe=Cexo(i,1:100-lag);
        end
        tabkc(p,i)=(cumtrapz(tempe)/cumtrapz(tempi))-1;
        tabkc(p,i+3)=(max(xcorr(tempi,tempe))/max(xcorr(tempi,tempi)))-1;
    end
end