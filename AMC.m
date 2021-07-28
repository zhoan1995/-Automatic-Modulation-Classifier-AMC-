clc
clear all

%% define the constellations that the source can simulate

t = 30000;
N = 128;
SNR = [0:5:30];
p = [4,4,6,6,8,2,2,2];
q = [0,2,1,3,2,1,0,2];
entropy16 = [2:0.1:4];
entropy64 = [2:0.1:6];
optimalpfigur16 = constellationTX(16,entropy16);
comp16 = sort(optimalpfigur16,2,'descend');
optimalpfigur64 = constellationTX(64,entropy64);
comp64 = sort(optimalpfigur64,2,'descend');
omega = [1:numel([entropy16 entropy64])];
for snr = 1:numel(SNR)
    snr
    for M = [16 64]
        if M == 16
            mean16 = meanvalue16(SNR(snr),entropy16,p,q);
        else
            mean64 = meanvalue64(SNR(snr),entropy64,p,q);
        end
    end
    
    
    
    meanall = [mean16 ; mean64];
    a = 1;
    
    %     for e = 1:numel(omega)
    for i = 1:numel(entropy16)
        M = 16;
%         datai = dataSimulation(M,(i-1)*N+1,N,comp16(i,:));
        for trial = 1:t
            datai = dataSimulation(M,(N*(trial-1))+1,N,comp16(i,:));
            Y = awgn(datai,SNR(snr),'measured');
            % load 'Y';
            for ii = 1:numel(p)
                %                 tosh1 = (Y .^(p(ii)-q(ii))).*(conj(Y).^q(ii));
                tosh = (Y .^(p(ii)-q(ii)))*((Y).^q(ii))';
                %                 sumetosh = sum(sum(tosh1));
                %                 empmoment16(1,ii) = (1/N)*sumetosh;
                empmoment16(1,ii) = (1/N)*tosh;
                
            end
            %% covariance
            L=numel(p);
            cuv16 = zeros (L,L);
            
            for tt = 1: numel(omega)
                
                chap16(tt,:) = empmoment16 - meanall(tt,:);
                rast16(tt,:) = conj(chap16(tt,:));
                
                for m = 1:L
                    for n = 1:L
                        cuv16 (m,n) = mean(chap16(tt,m)*rast16(tt,n));
                    end
                end
                
                %liklihood16
                C = empmoment16(1,:) - meanall(tt,:);
                F = ctranspose(C);
                U = (1./cuv16);
                f(tt,1) = (1/((pi^L)*det(cuv16)))*exp(-C*U*F);
                
            end
            [argmaxf, argmax] = max(f);
            
            montecarlo(trial)= argmax;
            
        end
        for O = 1:numel(omega)
            %             arg = 0;
            %             for i = 1:trial
            %                 if montecarlo(i) == O
            %                     arg = arg+1;
            %                 end
            %             end
            %             validprob(a,O) = arg/trial;
            validprob(a,O) = nnz(montecarlo == O)/trial;
            
        end
        a=a+1;
        
        
    end
    
    for i = 1:numel(entropy64)
        M = 64;
%         datai = dataSimulation(M,((i-1)*N)+1,N,comp64(i,:));

        for trial = 1:t
        datai = dataSimulation(M,(N*(trial-1))+1,N,comp64(i,:));
        Y = awgn(datai,SNR(snr),'measured');
            %          load 'Y';
            for ii = 1:numel(p)
                %                 tosh = (Y .^(p(ii)-q(ii))).*(conj(Y).^q(ii));
                tosh = (Y .^(p(ii)-q(ii)))*((Y).^q(ii))';
                %                 sumetosh = sum(sum(tosh));
                %                 empmoment64(1,ii) = (1/N)*sumetosh;
                empmoment64(1,ii) = (1/N)*tosh;
            end
            %% covariance
            L=numel(p);
            cuv64 = zeros (L,L);
            
            for tt = 1: numel(omega)
                
                chap64(tt,:) = empmoment64(1,:) - meanall(tt,:);
                rast64(tt,:) = conj(chap64(tt,:));
                
                
                for m = 1:L
                    for n = 1:L
                        cuv64 (m,n) =mean(chap64(tt,m)*rast64(tt,n));
                    end
                end
                
                %liklihood16
                C = empmoment64(1,:) - meanall(tt,:);
                F = ctranspose(C);
                U = (1./cuv64);
                f(tt,1) = (1/((pi^L)*det(cuv64)))*exp(-C*U*F);
                
            end
            [argmaxf, argmax] = max(f);
            
            montecarlo(trial)= argmax;
        end
        for O = 1:numel(omega)
%                         arg = 0;
%                         for i = 1:trial
%                             if montecarlo(i) == O
%                                 arg = arg+1;
%                             end
%                         end
%                         validprob(a,O) = arg/trial;
             validprob(a,O) = nnz(montecarlo == O)/trial;
        end
        a=a+1;
        
    end
    
    %     end
    
    probc = trace(validprob);
    
    
    accuracy(snr) = probc/numel(omega);
    
    save(['accuracy_snr_' num2str(SNR(snr)) '.mat'],'validprob')
    
end
figure()
plot(SNR,accuracy,'b-x')




