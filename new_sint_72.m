close all

n_loop=10;
perc_Out=0.03;

x_ruolo_assegnato=zeros(4,n_loop);
x_freq=zeros(4,n_loop);
x_precisione=zeros(4,n_loop);
ruolo=zeros(4,1);

% T =[1 0 0 0
%     0 0.8 0 0.2
%     0 0 1 0
%     0 0 0 1];

n_clique=2;
% n_clique=randsample(10,1,true)+2;
n_CS=3;
% n_CS=randsample(10,1,true)+2;

n_ES_max=5;
n_C_max=5;

n_B=5;  % Usi questo o le due righe sotto casuali
% n_B_max=n_clique*n_CS;
% n_B=randsample(n_B_max,1,true);

n_B_link_max=1; % In realtà divisi fra ES e C
% n_B_link_max=randsample(5,1,true);

n_C=randsample(n_C_max,n_clique,true)+2;
n_ES=randsample(n_ES_max,n_CS,true)+3;

n_C_tot=sum(n_C);
n_ES_tot=sum(n_ES);
n=n_CS+n_ES_tot+n_B+n_C_tot;

Gt=zeros(n,4,n_loop);
Ft=zeros(4,21,n_loop);
Vt=zeros(n,21,n_loop);
Gn=zeros(n,4,n_loop-1); % sarebbe Gt-1 calcolata dalla NMF
Tn=zeros(4,4,n_loop-1);
I=zeros(n,n_loop);


flag=zeros(n,1);

for k_loop=1:n_loop
% AL è adiacency list. Viene da A trasformata in logical
A=zeros(n);

% Center of Star
start=n_CS;
for i=1:n_CS
    for j=1:n_ES(i)
        A(i,j+start)=1;
    end
    start=start+n_ES(i);
end

% Clique
start=n-n_C_tot;
% start1=n-n_C_tot-n_B;
for k=1:n_clique
for i=1+start:start+n_C(k)
    for j=1:n_C(k)
        A(i,j+start)=1;
    end

%     for j=1:n_B
%         prob=rand;
%         if prob<0.1
%              A(i,j+start1)=1;
%         end
%     end
% Erano i B ma li metto tutti dai B
end
start=start+n_C(k);
end

% % Bridge
% start=n_CS;
% start1=n-n_C_tot;
% for i=n_CS+n_ES_tot+1:n-n_C_tot
%     prob=0.8;
%     conta=0;
%     for j=1:n_ES_tot
%         for k=1:n_CS
%             if j==n_ES(k)+1
%                 prob=0.8;
%             end
%         end
%         x=rand;
%         if x<prob
%              A(i,j+start)=1;
%              prob=0.01;
%              conta=conta+1;
%         end
%         if conta==n_B_link_max
%             break
%         end
%     end
%     conta=0;
%     prob=0.8;
%     for j=1:n_C_tot
%         if conta==n_B_link_max
%             break
%         end
%         for k=1:n_clique
%             if j==n_C(k)+1
%                 prob=0.8;
%             end
%         end
%         x=rand;
%         if x<prob
%              A(i,j+start1)=1;
%              prob=0.01;
%              conta=conta+1;
%         end
%     end
% 
% end

% Bridge

% Idea: per ogni bridge scelgo casualmente un certo numero di Stelle e un
% certo numero di Clique e all'interno di ognuna scelgo casualmente 1
% membro che collego.

for i=n_CS+n_ES_tot+1:n-n_C_tot

%     Questi sono gli indici di n_ES che indicano le stelle scelte.
    B_link_CS=randsample(n_CS,n_B_link_max);
%     Questi sono gli indici di n_C che indicano le clique scelte.
    B_link_C=randsample(n_clique,n_B_link_max);
     
    for k=1:n_B_link_max
        start=n_CS+sum(n_ES(1:B_link_CS(k)))-n_ES(B_link_CS(k));
        j=randsample(n_ES(B_link_CS(k)),1)+start;
        A(i,j)=1;
    end

    for k=1:n_B_link_max
        start=n-n_C_tot+sum(n_C(1:B_link_C(k)))-n_C(B_link_C(k));
        j=randsample(n_C(B_link_C(k)),1)+start;
        A(i,j)=1;
    end
end

% % Edge of stars
start=n_CS;
for k=1:n_CS
for i=1+start:start+n_ES(k)
    for j=1:n_ES(k)
        prob=rand;
        if prob<0.01
            A(i,j+start)=1;
        end
    end
end
start=start+n_ES(k);
end

%     Qui cambio il grafo con gli outliers

if k_loop>1
%     for i=1:n
%     scelta=0;
%     prob=rand;
%     if prob<perc_Out
%         A(i,:)=0;
%         A(:,i)=0;
%         scelta=randsample(4,1,true);
%         flag(i)=scelta;
% 
%         if scelta==1
%             k=randsample(n_CS,1,true);
%             start=n_CS+sum(n_ES(1:k))-n_ES(k);
%             for j=1:n_ES(k)
%                 A(i,j+start)=1;
%             end
% 
%         elseif scelta==2
%             k=randsample(n_CS,1,true);
%             A(i,k)=1;
% 
%         elseif scelta==3
%             B_link_CS=randsample(n_CS,n_B_link_max);
%             B_link_C=randsample(n_clique,n_B_link_max);
%      
%             for k=1:n_B_link_max
%                 start=n_CS+sum(n_ES(1:B_link_CS(k)))-n_ES(B_link_CS(k));
%                 j=randsample(n_ES(B_link_CS(k)),1)+start;
%                 A(i,j)=1;
%             end
% 
%             for k=1:n_B_link_max
%                 start=n-n_C_tot+sum(n_C(1:B_link_C(k)))-n_C(B_link_C(k));
%                 j=randsample(n_C(B_link_C(k)),1)+start;
%                 A(i,j)=1;
%             end
% 
%         elseif scelta==4
%             k=randsample(n_clique,1,true);
%             for j=1:n_C(k)
%                 A(i,j+(n-sum(n_C(k:n_clique))))=1;
%             end
%         end
%     end
%     end


    for i=1:n_ES_tot
        prob=rand;
        if prob<=perc_Out || flag(i+n_CS)~=0
            k=randsample(n_clique,1,true);
            flag(i+n_CS)=k;
            A(i+n_CS,:)=0;
            A(:,i+n_CS)=0;
            for j=1:n_C(k)
                A(i+n_CS,j+(n-sum(n_C(k:n_clique))))=1;
            end
            for j=1:i
                if flag(j+n_CS)==flag(i+n_CS)
                    A(i+n_CS,j+n_CS)=1;
                end
            end
        end
    end
end
% 
% if k_loop==10
%     for i=1:n_ES_tot
%         k=randsample(n_clique,1,true);
%         flag(i+n_CS)=k;
%         A(i+n_CS,:)=0;
%         A(:,i+n_CS)=0;
%         for j=1:n_C(k)
%             A(i+n_CS,j+(n-sum(n_C(k:n_clique))))=1;
%         end
%         for j=1:i
%             if flag(j+n_CS)==flag(i+n_CS)
%                 A(i+n_CS,j+n_CS)=1;
%             end
%         end
%     end
% end



% La rendo simmetrica e rimuovo 1 dalla diagonale (no auto link)
for i=1:n
    for j=1:n
        if A(i,j)~=A(j,i)
            A(i,j)=1;
            A(j,i)=1;
        end
    end
    A(i,i)=0;
end


AL=logical(A);
G=graph(AL);


[Gt(:,:,k_loop),Ft(:,:,k_loop),Vt(:,:,k_loop),max_G,I(:,k_loop),L,r] = DBMM(G,n);

% r


% if k_loop==1
%     G1=G_m;
%     F1=F_m;
%     V1=V_m;
% else
%     G2=G_m;
%     F2=F_m;
%     V2=V_m;
% end


% I(1:n_CS)
% I(n_CS+1:n_CS+n_ES_tot)
% I(n_CS+n_ES_tot+1:n-n_C_tot)
% I(n-n_C_tot+1:n)

if k_loop==1 || k_loop==n_loop

    x=linspace(1,n,n);
    % figure
    % plot(x,I);
    figure
    scatter(x,I(:,k_loop),'.');
    hold on
    xline([n_CS+1 n_CS+n_ES_tot+1 n-n_C_tot+1])


    figure
    h=plot(G);
    % Ora coloro: CS:y,     ES:r,    B:b,    C:g
    highlight(h,1:n_CS,'NodeColor','yellow');
    highlight(h,n_CS+1:n_CS+n_ES_tot,'NodeColor','red');
    highlight(h,n_CS+n_ES_tot+1:n-n_C_tot,'NodeColor','blue');
    highlight(h,n-n_C_tot+1:n,'NodeColor','green');
    % Qui coloro gli Outliers: per il momento non tengo conto di cosa sono ora
    % o di cosa erano, ma l'informazione c'è, quindi posso usarla se serve
    for i=1:n
        if flag(i)~=0
            highlight(h,i,'NodeColor','magenta');
            highlight(h,i);
        end
    end
end

% Ora cerco di calcolare quanto l'algoritmo sia preciso. Per farlo prima
% trovo quale ruolo ha assegnato ad ognuno dei quattro pattern,
% semplicemente tramite la moda. Poi calcolo la percentuale di volte che ha
% scelto la moda per quel pattern.

[x_ruolo_assegnato(1,k_loop),x_freq(1,k_loop)]=mode(I(1:n_CS,k_loop));
[x_ruolo_assegnato(2,k_loop),x_freq(2,k_loop)]=mode(I(n_CS+1:n_CS+n_ES_tot,k_loop));
[x_ruolo_assegnato(3,k_loop),x_freq(3,k_loop)]=mode(I(n_CS+n_ES_tot+1:n-n_C_tot,k_loop));
[x_ruolo_assegnato(4,k_loop),x_freq(4,k_loop)]=mode(I(n-n_C_tot+1:n,k_loop));
x_precisione(1,k_loop)=x_freq(1,k_loop)/n_CS;
x_precisione(2,k_loop)=x_freq(2,k_loop)/n_ES_tot;
x_precisione(3,k_loop)=x_freq(3,k_loop)/n_B;
x_precisione(4,k_loop)=x_freq(4,k_loop)/n_C_tot;

% for i=1:4
%     for j=1:4
%         if x_ruolo_assegnato(i,k_loop)==x_ruolo_assegnato(j,k_loop)&&i~=j
%             if x_precisione(i,k_loop)>=x_precisione(j,k_loop)
%                 
%             end
%         end
%     end
% end


% fine del loop
end 

% per Ax = B usi x = A\B
% T_alg=G1\G2
% T_confronto=T_ml/max(T_ml,[],'all')
% G2_prova=G1*T;
% [G_n,T_n]=nnmf(G2,4,'Algorithm','mult');

% [G1n,Tn,F] = NMF(Gt(:,:,2),r);

% Tn

% x_ruolo_assegnato
% x_precisione

for i=1:n_loop-1
    [Tn1(:,:,i),F_sum,k_sum] = NMF(Gt(:,:,i+1),r,Gt(:,:,i));
%     [Gn(:,:,i),Tn(:,:,i)] = nnmf(Gt(:,:,i+1),r,"algorithm","mult");
end
% Tn_media=sum(Tn,3)/(n_loop-1)
% Tn_media1=sum(Tn1,3)/(n_loop-1)

sommaT=sum(Tn1(:,:,1));
for k=1:r
    for j=1:r
        
        Tn1(k,j,1)=Tn1(k,j,1)/sommaT(j);
    end
end
Tn1(:,:,1)
Tn1(:,:,2)

for i=1:n_loop-1
    G_stacked1=[Gt(:,:,i)];
    G_stacked2=[Gt(:,:,i+1)];
end
[Tn_stacked,F_stacked,k_stacked] = NMF(G_stacked2,r,G_stacked1);
% [dump,Tn_stacked] = nnmf(G_stacked,r,"algorithm","mult");

Tn_stacked

sommaT=sum(Tn_stacked);
for k=1:r
    for j=1:r
        
        Tn_stacked(k,j)=Tn_stacked(k,j)/sommaT(j);
    end
end

Tn_stacked

% for k_loop=1:n_loop-1
%     ruolo=x_ruolo_assegnato(:,k_loop);
%     for i=1:4G_stacked2
%         if x_ruolo_assegnato(i,k_loop)~=x_ruolo_assegnato(i,k_loop+1)
%             ruolo_temp(i,k_loop)=17
%         end
%     end
% end






% T=zeros(4,4,n_loop-1);
% for k=1:n_loop-1
%     for i=1:n
%         T(I(i,k),I(i,k+1),k)=T(I(i,k),I(i,k+1),k)+1;
%     end
%     T_sum=sum(T(:,:,k),2);
%     for i=1:4
%         T(i,:,k)=T(i,:,k)/T_sum(i);
%     end
% end
% T_media=sum(T,3)/(n_loop-1)




G_vera=zeros(n,4,10);
G_vera_dump=zeros(n,4,10);
Tn_vera=zeros(4,4,n_loop-1);
Tn1_vera=zeros(4,4,n_loop-1);

for k=1:10
    for i=1:n_CS
        G_vera(i,:,k)=[1 0 0 0];
    end
    for i=n_CS+1:n_CS+n_ES_tot
        G_vera(i,:,k)=[0 1 0 0];
    end
    for i=n_CS+n_ES_tot+1:n_CS+n_ES_tot+n_B
        G_vera(i,:,k)=[0 0 1 0];
    end
    for i=n-n_C_tot+1:n
        G_vera(i,:,k)=[0 0 0 1];
    end
end

for k=1:10
    for i=n_CS+1:n_CS+floor((n_ES_tot/10)*k)
        G_vera(i,:,k)=[0 0 0 1];
    end
end

for k=1:10
    for i=n-n_C_tot+1:n-n_C_tot+floor((n_C_tot/10)*k)
        G_vera(i,:,k)=[0 1 0 0];
    end
end


for i=1:n_loop-1
    [G_vera_dump(:,:,i),Tn_vera(:,:,i)] = nnmf(G_vera(:,:,i+1),r,"algorithm","mult");
    [Tn1_vera(:,:,i)] = NMF(G_vera(:,:,i+1),r,G_vera(:,:,i));
end
Tn_media_vera=sum(Tn_vera,3)/(n_loop-1)
Tn1_media_vera=sum(Tn1_vera,3)/(n_loop-1)



beep











