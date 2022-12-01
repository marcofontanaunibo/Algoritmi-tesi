close all

% K=numero di comunità se sono uguali fra i due snapshot, sennò usi K1 e K2
K=8;
N=1250;
Ntot=N*K;
% Sia NOut il numero di outliers, una percentuale di Ntot
NOut=Ntot/10;

mu1=[1.5 1.5;3 2;3 1;4 2;4 1;1 4;2 4;3 4];
mu2=[1 1;1 2;2 1;3.5 1.5;2 2;1 4;2 4;3 4];
sigma1 = cat(3,[0.2 0.2],[0.05 0.05],[0.05 0.05],[0.05 0.05], ...
    [0.05 0.05],[0.05 0.05],[0.05 0.05],[0.05 0.05]); %cat fa in modo di avere un vettore di matrici
sigma2 = cat(3,[0.05 0.05],[0.05 0.05],[0.05 0.05],[0.2 0.2], ...
    [0.05 0.05],[0.05 0.05],[0.01 0.01],[0.1 0.05]);



nloop=10 ; % quante volte lo fai girare

time=zeros(nloop,1);
% AUC=zeros(nloop,1);
AUC_ROC=zeros(nloop,1);
AUC_ROC_controllo=zeros(nloop,1);
% n_Out_beccati=zeros(nloop,1);
% n_Out_beccati_controllo=zeros(nloop,1);
% alfa_ideale=zeros(nloop,1);
% alfa_indice_ideale=zeros(nloop,1);
% u_vector=zeros(nloop,1);

for k_loop=1:nloop
%     tic

gm1 = gmdistribution(mu1,sigma1);
gm2 = gmdistribution(mu2,sigma2);


% genero i punti:
X1=[];
X2=[];
for i=1:K
    X1 = [X1; mvnrnd(mu1(i,:),sigma1(1,:,i),N)];
    X2 = [X2; mvnrnd(mu2(i,:),sigma2(1,:,i),N)];
end


P = posterior(gm1,X1);  % prob condizionata che pt appartenga a cluster,
                      % cioè insomma P (Q).

Q = posterior(gm2,X2); 


% Ora calcolo S senza la presenza di outliers: è per il controllo con i
% dati sintetici.

S_noOut = mldivide(P,Q);


% inserisco forzatamente gli outliers:

iOut = randsample(Ntot,NOut); % gli indici degli outliers.
not_iOut=[1:Ntot];
not_iOut=not_iOut';
not_iOut(iOut)=[]; % gli indici non outliers

for i=1:NOut
    [q_max,j_max]=max(Q(iOut(i),:));
    [q_min,j_min]=min(Q(iOut(i),:));
    Q(iOut(i),j_max)=q_min;
    Q(iOut(i),j_min)=q_max*(1/2);
    j=randsample(K,1);
    Q(iOut(i),j)=Q(iOut(i),j)+q_max/2;
end
% Ho inserito quindi NOut outliers scambiando le entrate di Q.


% Ora calcolo A del metodo a due passi per il confronto

A_controllo=zeros(Ntot,K);
for o=1:Ntot
    for j=1:K
        A_controllo(o,j)=(Q(o,j)-P(o,:)*S_noOut(:,j))^2;
    end
end
A_sum=sum(A_controllo,2);
[M_controllo,I_controllo] = sort(A_sum,'descend');


[A,S,M,I,u,t] = ECO(P,Q);

% u_vector(k_loop)=u;


% Ora confronto gli outliers:

% Precisione=zeros(NOut,1);
% c=0;
% d=0;
% for i=1:NOut
%     for j=1:i
%         if I(j)-I_controllo(i)==0
%             c=c+1;
%         end
%         if (I(i)-I_controllo(j)==0 && j~=i)
%             c=c+1;
%         end
%     end
%     d=d+c;
%     c=0;
%     Precisione(i)=d/i;
% end

% % grafico di precisione su outlier ranking
% x=linspace(1,NOut,NOut);
% plot(x,Precisione);


% calcolo "AUC" come nel 70
% AUC(k_loop)=sum(Precisione)/NOut;


% Calcolo AUC-ROC facendo variare alfa come valore di A (M):
TPR=ones(Ntot*3,1);
FPR=ones(Ntot*3,1);
alfa=max(M);
k=1;
flag=0;

while alfa>0

    i=length(find(M>=alfa));

    TP = length(intersect(I(1:i),iOut));

    FP = length(intersect(I(1:i),not_iOut));

    FN = length(intersect(I(i:Ntot),iOut));

    TN = length(intersect(I(i:Ntot),not_iOut));

    TPR(k)=TP/(TP+FN);
    FPR(k)=FP/(FP+TN);

%     if flag==0 && TPR(k)==1
%         alfa_ideale(k_loop)=alfa;
%         flag=1;
%     end
    
    alfa=alfa-0.01;
    k=k+1;
end

% figure
% plot(FPR,TPR);

AUC_ROC(k_loop)=sum(TPR(1:k))/k;

% for i=1:Ntot
%     if M(i)<alfa_ideale(k_loop)
%         alfa_indice_ideale(k_loop)=i;
%         break
%     end
% end

% n_Out_beccati(k_loop)=length(intersect(I(1:NOut),iOut));
% n_Out_beccati_controllo(k_loop)=length(intersect(I_controllo(1:NOut),iOut));



% Calcolo AUC-ROC del confrontovfacendo variare alfa come valore di A (M):
TPR=ones(Ntot*3,1);
FPR=ones(Ntot*3,1);
alfa=max(M_controllo);
k=1;
flag=0;

while alfa>0

    i=length(find(M_controllo>=alfa));

    TP = length(intersect(I_controllo(1:i),iOut));

    FP = length(intersect(I_controllo(1:i),not_iOut));

    FN = length(intersect(I_controllo(i:Ntot),iOut));

    TN = length(intersect(I_controllo(i:Ntot),not_iOut));

    TPR(k)=TP/(TP+FN);
    FPR(k)=FP/(FP+TN);
    
    alfa=alfa-0.01;
    k=k+1;
end

% figure
% plot(FPR,TPR);

AUC_ROC_controllo(k_loop)=sum(TPR(1:k))/k;






% Ora i grafici

% figure
% %hold on
% scatter(X1(:,1),X1(:,2),10,'.');  % mostro i punti
% %scatter(X(I(1),1),X(I(1),2),50,'o','r','filled'); % i maggiori outliers
% %scatter(X(I(2),1),X(I(2),2),50,'o','y','filled');
% %scatter(X(I(3),1),X(I(3),2),50,'o','m','filled');
% %scatter(X(I(4),1),X(I(4),2),50,'o','k','filled');
% figure
% scatter(X2(:,1),X2(:,2),10,'.');  % mostro i punti
% 
% 
% % mostro i clusters
% 
% figure
% scatter(X1(:,1),X1(:,2),5,'.') % Scatter plot with points of size 10
% hold on
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm1,[x0 y0]),x,y);
% fcontour(gmPDF,[0 6 0 6]) % 0 6 0 6 è quanto grande è il grafico
% c1 = colorbar;
% ylabel(c1,'Probability Density Function');
% 
% figure
% scatter(X2(:,1),X2(:,2),5,'.') % Scatter plot with points of size 10
% hold on
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm2,[x0 y0]),x,y);
% fcontour(gmPDF,[0 6 0 6])
% c1 = colorbar;
% ylabel(c1,'Probability Density Function');



time(k_loop)=t;
k_loop
end
% u_media=mean(u_vector)
% x_media=mean(AUC)
% x_varianza=var(AUC);
x_AUC_ROC=mean(AUC_ROC)
x_AUC_ROC_controllo=mean(AUC_ROC_controllo)
x_AUC_ROC_var=var(AUC_ROC)
x_AUC_ROC_controllo_var=var(AUC_ROC_controllo)
% x_n_Out_beccati=mean(n_Out_beccati)
% x_n_Out_beccati_controllo=mean(n_Out_beccati_controllo)
% x_alfa_ideale=mean(alfa_ideale)
% x_alfa_indice_ideale=mean(alfa_indice_ideale)
x_time_tot=sum(time);
x_time_tot*5


% Mostro il variare dei valori di a (aggregati e ordinati):
% x=linspace(1,Ntot,Ntot);
% % figure
% % plot(x,M_controllo);
% % % hold on
% figure
% plot(x,M);
% figure
% plot(1:(NOut*3),M(1:(NOut*3)));

beep
