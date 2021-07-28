%
% clc
% clear all
function datai = dataSimulation(M,start,N,comp,data)
%% function for simulate data

 rng();
 data = rand(1,N);
% load 'simuldata128';
% if N == 128
%     data = simuldata128(start:N+start-1);
% end

%% number of layers in constellation

layers = unique(comp);
layers=flip(layers);
layernum = numel(layers);

%% number of symbols in each constellation

layersymbol = zeros(1,layernum);
layerprob = zeros(1,layernum);
ii = 1;
j = 1;
x = 1;
for k = 1:layernum
    layersymbol(k) = nnz(comp == layers(k));
    layerprob(k)= layersymbol(k)*layers(k);
end


%% assign symbols
data1 = zeros(1,numel(layerprob));
for i = 1:numel(data)
    k = 1;
    p=0;
    x = 0;
    
    for j = 1: numel(layerprob)
        
        if data(i)>p && data(i)<layerprob(j)+p
            data1(j) = data1(j)+1 ;
        end
        p = layerprob(j)+p;
    end
end

symbolnum = data1./layersymbol;
tedadsymbol = zeros (1,M);
n = 1;
z =1;
for i = 1:numel(symbolnum)
    
    if mod(data1(i),layersymbol(i)) == 0
        tedadsymbol(z:layersymbol(i)+z-1) = symbolnum(i);
        z = layersymbol(i)+z;
    else
        
        k = floor(symbolnum(i));
        L = layersymbol(i);
        

        tedadsymbol (1,z:(z+L-1)) = k;
        z = z+L;
        v = data1(i)-(k*(L));

        tedadsymbol(1,z-v : z-1) = tedadsymbol (1,z-v : z-1) +1;

        
    end
    
end

data2 = zeros(1,N);
n = 1;

if M == 16
    symbols = [5,13,15,7,1,4,12,9,11,14,6,3,0,8,10,2];
    for i = 1:M
        data2(1,n: n+tedadsymbol(i)-1) = symbols(i);
        n = tedadsymbol(i)+n;
    end
end

if M == 64
    symbols = [18,50,22,54,19,51,58,62,55,23,30,26 ,59,63,31,27 ,17,10,14,21,53,46,42,49 ,25,11,15,29,61,47,43,57,16,9,2,6,13,20,52,45,38,34,41,48 ,24,3,7,28,60,39,35,56 ,8,1,5,12,44,37,33,40,0,4,36,32];
    for i = 1:M
        data2(1,n: n+tedadsymbol(i)-1) = symbols(i);
        n = tedadsymbol(i)+n;
    end
end
datai = qammod(data2,M);
end