function [A,S,M,I,u] = ECO(P,Q)

% N è il numero degli oggetti,
% X1 e X2 sono i due snapshot dei temporal dataset,
% K1 e K2 sono il numero di comunità in X1 e X2
% P e Q sono le relative matrici di appartenenza probabilistiche, ovvero
% ogni elemento di P rappresenta come un particolare oggetto appartenga
% a una particolare comunità.
% P è N*K1 e ogni elem di P sta fra 0 e 1, con somma su righe =1.
% P e Q sono l'input.
% u=1

tic
% tiene traccia del tempo impiegato dall'algoritmo

[Ntot,K1]=size(P);
[~,K2]=size(Q);

u=1;

% S è la matrice di corrispondenza, S è K1*K2, sempre con somma su righe=1,
% che rappresenta il cambiamento delle comunità da K1 a K2,
% va calcolata

S=zeros(K1,K2);
S_0=zeros(K1,K2);

% A è la matrice N*K2 che rappresenta l'outlierness score per ogni oggetto,
% va calcolata

A=zeros(Ntot,K2);
A_0=zeros(Ntot,K2);

% Ora si fa una doppia minimizzazione tenendo conto di S e A

% u è la somma di outlierness score, viene presa = 1.
% Poi u si aggiorna dopo.

% Inizializzo A, S e f, funzione obiettivo.
% Salvo i valori iniziali per non doverli ricalcolare ogni volta.

for i=1:K1
    for j=1:K2
        S_0(i,j)=1/K2;
    end
end

for i=1:Ntot
    for j=1:K2
        A_0(i,j)=1/(Ntot*K2);
    end
end

fo_0=0;
for o=1:Ntot
    for j=1:K2
        fo_0=fo_0+(log(1/A_0(o,j))*((Q(o,j)-P(o,:)*S_0(:,j))^2));
    end
end

% questa è la fz obiettivo che riscrivo di volta in volta
fo=zeros(2,1);


% Qui parte l'algoritmo  

eps=1e-06;      % tolleranza

for k_alg=1:10
    
    k_while=0;  % per tener conto dei cicli all'interno del while

    S=S_0;
    A=A_0;

    % fz obiettivo che riscrivo
    fo(1)=fo_0;
    indice=2;
    

        while true
            k_while=k_while+1;
            % per la convergenza devo fermarmi
            % quando la differenza di valore della
            % funzione obiettivo è inferiore alla tolleranza. 

            % aggiorna A

            den_A=0;
            for o=1:Ntot
                for j=1:K2
                    A(o,j)=((Q(o,j)-P(o,:)*S(:,j))^2)*u;
                    den_A=((Q(o,j)-P(o,:)*S(:,j))^2)+den_A;
                end
            end
            
            A=A/den_A;

            % aggiorna S

            % questo è per trovare B(i)

            B=zeros(K1,1);
            num=zeros(K1,K2);
            den=zeros(K1,K2);
            prodBj=1;
            sommACj=0;
            sommCj=0;

            for i=1:K1
                for j=1:K2
                    for o=1:Ntot
                        num(i,j)=num(i,j)+2*log(1/A(o,j))*P(o,i)*(Q(o,j)-P(o,:)*S(:,j)+P(o,i)*S(i,j));
                        den(i,j)=den(i,j)+2*log(1/A(o,j))*P(o,i)*P(o,i);
                    end
                    prodBj=prodBj*den(i,j);
                end
                for j=1:K2
                    sommACj=sommACj+num(i,j)*prodBj/den(i,j);
                    sommCj=sommCj+prodBj/den(i,j);
                end
                
                B(i)=(sommACj-prodBj)/sommCj;
                
                prodBj=1;
                sommACj=0;
                sommCj=0;
                
            end

            % ora finalmente trovo sij
            for i=1:K1
                for j=1:K2
                    S(i,j)=(num(i,j)-B(i))/den(i,j);
                end
            end
            
            % calcolo fo(i)
            % dove se i è 1 poi lo metto a 2, se i è 2 lo metto a 1
            % indice=2 di partenza.

            for o=1:Ntot
                for j=1:K2
                    fo(indice)=fo(indice)+(log(1/A(o,j))*((Q(o,j)-P(o,:)*S(:,j))^2));
                end
            end

            % condizione di convergenza
            if abs(fo(2)-fo(1))<eps
                break
            end

            if indice==1
                indice=2;
            else indice=1;
            end
            
            % la ferma a una certa per sicurezza
            if k_while==15
                break
            end

        end
        
    % aggiorno u
    
    u=0;
    for o=1:Ntot
        for j=1:K2
            u=u+((Q(o,j)-P(o,:)*S(:,j))^2);
        end
    end
    q = max(Q, [], 'all');
    u=u/(q^2);

%     Ora algoritmo ricomincia da capo, unica differenza
%     è il il valore di u
end


% ordina A per trovare outliers per ogni comunità,
% cioè i più grandi valori di A

A_sum=sum(A,2);
[M,I] = sort(A_sum,'descend');


time=toc;
end