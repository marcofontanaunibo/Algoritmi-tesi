function [G_matrix,F,V,max_G,I,L,r] = DBMM(G,n)

% input: grafo. n è il numero di nodi

tic

% Trovare caratteristiche:
% Si usano il degree di ogni nodo e
% Egonet Features: si contano i link all'interno dell'Egonet 
% e quelli che collegano l'Egonet all'esterno.
% NB: Egonet è il nodo e i suoi vicini a passo 1.

% Quindi serve:
% deg del nodo
% EFI: Egonet Feature interno
% EFO: Egonet Feature esterno (out)
% Sono entrambi vettori di dim n: un egonet per nodo
% Aggiungo poi iterativamente (mi fermo al terzo grado):
% media e somma di ogni caratteristica

EFI=zeros(n,1);
EFO=zeros(n,1);
deg=degree(G);
media_deg_vicino=zeros(n,1);
somma_deg_vicino=zeros(n,1);
media_EFI_vicino=zeros(n,1);
somma_EFI_vicino=zeros(n,1);
media_EFO_vicino=zeros(n,1);
somma_EFO_vicino=zeros(n,1);

for i=1:n
%     EFI
    vicino=neighbors(G,i);
    egonet=[vicino;i];
    H = subgraph(G, egonet);
    EFI(i,1)=numedges(H);
    
%     EFO
    tutti=[];
    for k=1:size(vicino)
        vicini2=neighbors(G,vicino(k)); % vicini dei vicini di i
        tutti=[tutti;vicini2];
    end
    for j=1:size(egonet)
        tutti = tutti(tutti~=egonet(j));
    end
    [rig,~]=size(tutti);
    EFO(i,1)=rig;

%     Media del degree dei vicini
    x=degree(G,vicino);
    media_deg_vicino(i,1)=mean(x);
    somma_deg_vicino(i,1)=sum(x);
end
for i=1:n
    vicino=neighbors(G,i);
    x=EFI(vicino);
    y=EFO(vicino);
    media_EFI_vicino(i,1)=mean(x);
    somma_EFI_vicino(i,1)=sum(x);
    media_EFO_vicino(i,1)=mean(y);
    somma_EFO_vicino(i,1)=sum(y);
end


% La matrice scoperta è V ed è n*f, dove
% f è il numero di feature per ogni nodo.

f=21;
V=zeros(n,f);
V(:,1)=EFI;
V(:,2)=EFO;
V(:,3)=deg;
V(:,4)=media_deg_vicino;
V(:,5)=somma_deg_vicino;
V(:,6)= media_EFI_vicino;
V(:,7)= somma_EFI_vicino;
V(:,8)= media_EFO_vicino;
V(:,9)= somma_EFO_vicino;

grado2=zeros(n,12);
for j=1:6
    for i=1:n
    vicino=neighbors(G,i);
    grado2(i,j)=mean(V(vicino,j+3));
    grado2(i,j+6)=sum(V(vicino,j+3));
    end
end
V(:,10:21)=grado2;

% Ora troviamo i nodi:
% dobbiamo trovare G e F tc G*F approssima V, con
% G una matrice n*r e F r*f, dove le righe di G rappresentano
% l'appartenenza di un nodo ad ogni ruolo,
% r= numero di ruoli, e
% ogni colonna di F rappresenta come l'appartenenza a uno specifico ruolo
% contribuisca a stimare i valori delle features.

% G_matrix=zeros(n,r);
% F=zeros(r,f);

% Model selection: vogliamo non avere r preselezionato:
% usiamo il MDL per trovare l'r ideale:
% Minimizziamo L=M+E, con:
% M=r(n+f);
% E=somma...,



% Modo 1: trovare r ideale
% 
% L=zeros(f,1);
% r=0;
% for r=1:f
% 
%     [G_matrix,F] = NMF(V,r);
% % [G_matrix,F] = nnmf(V,r,"algorithm","mult");
% 
%     M=r*(n+f); % log2(n)*
% 
%     E=0;
%     GF=G_matrix*F;
%     for i=1:n
%         for j=1:f
%             if GF(i,j)==0
%                 GF(i,j)=0.1;
%             end
%             if V(i,j)==0
%                 V(i,j)=0.1;
%             end
%             E=E+V(i,j)*log(V(i,j)/GF(i,j))-V(i,j)+GF(i,j);
%         end
%     end
% 
%     L(r)=M+E;
%     if r>1
%         if L(r)-L(r-1)>=0
%             r=r-1;
%             [G_matrix,F] = NMF(V,r);
% % [G_matrix,F] = nnmf(V,r,"algorithm","mult");
%             break;
%         end
%     end
% end

% fine Modo 1



% Modo 2: r fissato

r=4;
[G_matrix,F] = nnmf(V,r,"algorithm","mult");
L=0;

% Fine modo 2



[max_G,I]=max(G_matrix,[],2);

time=toc;
end