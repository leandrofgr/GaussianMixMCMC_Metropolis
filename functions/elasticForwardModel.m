function [G] = elasticForwardModel(uVp, uVs, wavelets, theta, type)

A = [];

    for i=1:length(theta)

        a = ones(length(uVs),1)*0.5*(1 + tand(theta(i))^2);
        b = -4*((uVs./uVp).^2)*sind(theta(i))^2;
        if nargin > 4
            % Fatti
            c = - 0.5*(tand(theta(i)).^2 + b);  
        else
            % Aki & Richards
            c = 0.5*(1 + b);    
        end  
        
         a(1) = [];
         b(1) = [];
         c(1) = [];
          
        A = [A ; diag(a) diag(b) diag(c)];
        
        n = size(wavelets,1);
        nl2 = int32((n-1)/2); 
        Saux = convmtx(wavelets(:,i),length(uVp)-1);
        Saux = Saux(nl2+1:end - nl2,:);
        
        sSaux = size(Saux);
        S(sSaux(1)*(i-1)+1:sSaux(1)*i,sSaux(2)*(i-1)+1:sSaux(2)*i) = Saux;
        
    end
    
    aux=1;
    aux1=0;
    for j = 1:3
        for i=1:length(uVp)-1
            D(aux,aux+aux1) = -1;
            D(aux,aux+aux1+1) = 1;
            aux=aux+1;
        end
        aux1 = aux1 + 1;
    end

    G = S*A*D;
    
end