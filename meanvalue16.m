% clc
% clear all

function [config16] = meanvalue16(SNR,entropy,p,q)
%% function for defining mean value for different moments in AMC
% p = [4,4,2,8,6];
% q = [0,3,1,0,2];
% SNR = 30;
% entropy = [2:0.1:4];
M =16;
L = sqrt(M);

%% Maxwell-Boltzman distribution for 1000 different amount of lambda
load('result16');

%% find the best probabilities based on average energy for each entropy
if isequal(entropy,[2:0.1:4]) == 1
    load('optimalfigure16[2-0.1-4].mat')
end
if isequal(entropy,[2:0.5:4])== 1
    load('optimalfigure16[2-0.5-4].mat')
end
if isequal(entropy,[2:4]) ==1
    load('optimalfigure16[2-4].mat')
end


%% buil data stream in constellation



comp = sort(optimalpfigur,2,'descend');

N = 300;
s = 1;
% load 'traindata';
rng();
datatrain = rand(1,300000);
tosh1 = zeros(numel(p),N);
mupq = zeros(numel(entropy),numel(p));
mianpq = zeros(1,1);
for hhhh = 1:size(optimalpfigur,1)
    for kk = 1:(numel(datatrain)/N)
        
        data = datatrain((1+(kk-1)*N) : N*kk);
        % number of layers in constellation
        %         layernum = 1;
        %         for i = 1: M-1
        %             if comp(hhhh,i) ~= comp(hhhh,i+1)
        %                 layernum = layernum+1;
        %             end
        %         end
        layers = unique(comp(hhhh,:));
        layers=flip(layers);
        layernum = numel(layers);
        
        
        % number of symbols in each constellation
        
        layersymbol = zeros(1,layernum);
        layerprob = zeros(1,layernum);
        ii = 1;
        j = 1;
        x = 1;
        
        for k = 1:layernum
            layersymbol(k) = nnz(comp(hhhh,:) == layers(k));
            layerprob(k)= layersymbol(k)*layers(k);
        end
        
        %         for i = 1:M-1
        %             if j > layernum
        %                 break
        %             end
        %             w = 1;
        %             for ii = 1:M
        %                 if comp(hhhh,x) == comp(hhhh,ii+x)
        %                     w = w+1;
        %                     if x+ii == M
        %                         layersymbol(j)= w;
        %                         layerprob (j) = w*comp(hhhh,x);
        %                         j = j+1;
        %                         break
        %                     end
        %                 else
        %                     layersymbol(j)= w;
        %                     layerprob (j) = w*comp(hhhh,x);
        %                     j = j+1;
        %
        %                     break
        %                 end
        %             end
        %             x = ii+x;
        %
        %
        %         end
        
        
        % assign symbols
        data1 = zeros(1,numel(layerprob));
        for i = 1:numel(data)
            k = 1;
            b=0;
            x = 0;
            
            for j = 1: numel(layerprob)
                
                if data(i)>b && data(i)<layerprob(j)+b
                    data1(j) = data1(j)+1 ;
                    
                end
                b = layerprob(j)+b;
            end
        end
        
        symbolnum = data1./layersymbol;
        tedadsymbol = zeros (1,M);
        n = 1;
        z=1;
        cc = mod(data1 , layersymbol);
        for i = 1:numel(symbolnum)
            
            if cc(i) == 0
                
                %                 for j = n: layersymbol(i)+n-1
                %                     tedadsymbol(j) = symbolnum(i);
                %                     n = n+1;
                %                 end
                tedadsymbol(n:layersymbol(i)+n-1) = symbolnum(i);
                n = layersymbol(i)+n;
            else
                
                k = floor(symbolnum(i));
                L = layersymbol(i);
                
                tedadsymbol(1,n:(n+L-1)) = k;
                n = n+L;
                v = data1(i)-(k*(L));
                tedadsymbol (1,n-v:n-1) = tedadsymbol (1,n-v:n-1) +1;
                %
                
                %         n = n+1;
                
            end
            
        end
        
        
        data2 = zeros(1,N);
        n = 1;
        
        symbols = [5,13,15,7,1,4,12,9,11,14,6,3,0,8,10,2];
        for i = 1:M
            data2(1,n:n+tedadsymbol(i)-1) = symbols(i);
            n = n+tedadsymbol(i);
        end
        
        
        
        datai = qammod(data2,M);
        datatrain16 = awgn(datai,SNR);
        
        for j = 1:numel(p)
            
            tosh1(j,:) = (datatrain16 .^(p(j)-q(j))) .* (conj(datatrain16) .^q(j));
            mianpq(kk,j) = mean(tosh1(j,:));
            
            
        end
        
    end
    mupq(hhhh,:) = mean(mianpq);
    %     for j = 1:N
    %         test(hhhh,j) = (abs(datatrain16(1,j)))^2;
    %     end
    %     power(hhhh,1)=sum(test(hhhh,:))/N;
end

config16 = mupq;
end

