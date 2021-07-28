% 
% clc
% clear all
  function [optimalpfigur] = constellationTX(M,entropy)
% M =64;
% entropy = [2:6];
L = sqrt(M);
m = log2(M);
l = (L/2)^2;
for i = 1:1:L
    A(1,i) = (2*i-1-L);
    B=A';
end
r=0;
ma=size(A,1);
mb=size(B,1);
[a,b]=ndgrid(1:ma,1:mb);
product = [A(a,:),B(b,:)];
c = product(:,L+1).^2;
product(:,L+1)=[];
product = product.^2;

amplitude = sqrt(c+product);
power2 = amplitude.^2;
% amp= amplitude(1:L/2,1:L/2);
quarter = amplitude(1:end/2,1:end/2);


row=1;
epsilon = 0.0009;
result = zeros (1,3+M);

%% Maxwell-Boltzman distribution for 10000 different amount of lambda
for lambda = 0:0.0001:1
    denom = sum (sum (exp(-power2*lambda)));
    %     for i = 1:L
    %         for j= 1:L
    %             prob(i,j) = exp((-power2(i,j)*lambda))/denom;
    %             h(i,j) = prob(i,j)*(log2(1/prob(i,j)));
    %             Energy(i,j) = (prob(i,j))*(amplitude(i,j)^2);
    %
    %         end
    %
    %
    %     end
    prob(1:L,1:L) = exp((-power2(1:L,1:L)*lambda))/denom;
    h(1:L,1:L) = prob(1:L,1:L).*(log2(1./prob(1:L,1:L)));
    Energy(1:L,1:L) = (prob(1:L,1:L)).*(amplitude(1:L,1:L).^2);
    
    H = sum(h,'all');
    avgE = (sum(Energy,'all'))/M;
    result(row,1) = H;
    result(row,2) = avgE;
    
    
    probs = reshape(prob,[1,M]);
    
    
    %     for ii = 3:1:M+2
    %         result(row,ii) = probs (1,ii-2);
    %     end
    result(row,3:M+2)=probs (1,:);
    
    row = row+1;
    
    
end
row = 1;

result = flip(result);

%% find the best probabilities based on average energy for each entropy
optimalp = zeros(numel(entropy),M+3);
optimalpfigur = zeros(numel(entropy),M);

for hh = entropy
    H_target = hh;
    H_down = H_target - epsilon;
    H_up = H_target + epsilon;
    accept = zeros(1,l+3);
    g = 1;
    
    
    
    for ii = 1:10001
        if H_target == 2
            %         for jj = 1:M+3
            %
            %             accept (g,jj) = result (1,jj);
            %             accept(g,M+3) = H_target;
            %         end
            accept (g,1:M+3) = result (1,1:M+3);
            accept(g,M+3) = H_target;
            g = g+1;
            break
        end
        if H_target == m
            break
        else
            if  result(ii,1) >= H_down && result(ii,1) <= H_up
                %                 for jj = 1:M+3
                %
                %                     accept (g,jj) = result (ii,jj);
                %                     accept(g,M+3) = H_target;
                %                 end
                accept (g,1:M+3) = result (ii,1:M+3);
                accept(g,M+3) = H_target;
                g=g+1;
                
            end
        end
    end
    
    if H_target == m
        %         for jj = 1:M+3
        %             accept (g,jj) = result (end,jj);
        %             accept(g,M+3) = H_target;
        %         end
        accept (g,1:M+3) = result (end,1:M+3);
        accept(g,M+3) = H_target;
    end
    
    energy = accept(:,2);
    min_Energy = min(energy(energy>0));
    for c = 1:size(accept,1)
        if  accept(c,2) == min_Energy
            r=r+1;
            optimalp (r,:)  = accept(c,:); %the results are here
            optimalpfigur (r,:) = optimalp (r,3:M+2);
        end
    end
    
    
    
    
    row = row+1;
    
    
end
  end