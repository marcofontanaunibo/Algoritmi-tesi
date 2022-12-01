function [H,F,k] = NMF(V,r,W)

% W è Gt-1  cioè n*r
% H è T     cioè r*m
% V è Gt    cioè n*m

[n,m]=size(V);
% W=ones(n,r);
% W_save=W;
H=eye(r,m);

n_loop=100;
tol=1e-9;

F=zeros(n_loop,1);

for k=1:n_loop

%     W=W_save;
    WH=W*H;
    for a=1:r

        for i=1:n
            somma=0;
            %             WH=W*H;
            for u=1:m
                somma=somma+V(i,u)*H(a,u)/WH(i,u);
            end
            W(i,a)=W(i,a)*somma;
            sommaW=sum(W);
            W(i,a)=W(i,a)/sommaW(a);
        end



        for u=1:m
            somma=0;
            %             WH=W*H;
            for i=1:n
                somma=somma+W(i,a)*V(i,u)/WH(i,u);
            end
            H(a,u)=H(a,u)*somma;

        end
    end
    sommaT=sum(H);
    for k=1:r
        for j=1:r

            H(k,j)=H(k,j)/sommaT(j);
        end
    end

%     WH=W*H;
%     logWH=log(WH);
    
    for i=1:n
        for u=1:m
            F(k,1)=F(k,1)+V(i,u)*log(V(i,u)/WH(i,u))-V(i,u)+WH(i,u);
        end
    end

    if k>2
        if abs(F(k)-F(k-1))<=tol
            return
        end
        if F(k)-F(k-1)>=0 && F(k-1)-F(k-2)<=0
            return
        end
    end

end
end




