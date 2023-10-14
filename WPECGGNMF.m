function [drb,B,V] = WPECGGNMF(mix, nb, L, derta, fftSize, shiftSize, ite, itNMF, beta, p, refMic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeNow = date;
if  strcmp(timeNow(8:end),'2023')

    [X, window] = STFT( mix, fftSize, shiftSize, 'hamming' );

    [F,T,M] = size(X);
    X_MFT = permute(X,[3,1,2]);
    X_MTF = permute(X,[3,2,1]);
    % X_TMF = permute(X,[2,3,1]);
    %初始化普参数
    B = max( rand( F, nb), eps );
    V = max( rand( nb, T ), eps );
    R(:,:) = B(:,:)*V(:,:);

    %构造信号
    xbar(:,:,1:derta+L) = zeros(F,M*L,derta+L);   % 构成预测滤波器的输入信号  F L T
    for t=derta+L+1:T
        for ail = 1:L
            xbar(:,(ail-1)*M+1:ail*M,t) = squeeze(X(:,t-derta-ail,:));
        end
    end
    fprintf('Iteration:    ');


    % 2. 迭代开始
    for iter = 1:ite
        fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', iter,ite);
        % (a) 去混响
        if(iter==1)
            for f = 1:F
                xEarly(f,:) = X_MFT(1,f,:); % 初值
            end
            P = max(abs(xEarly).^beta,1e-8);
            garma = abs(xEarly);
        else
            for f=1:F
                xEarly(f,:) = squeeze(X_MFT(1,f,:)).' - (G(:,f)'*(squeeze(xbar(f,:,:))));
            end
            P = max(abs(xEarly).^beta,1e-8);
            garma = abs(xEarly);
        end

        % 计算声源方差
        %%%%% Update B %%%%%
        for i = 1:itNMF
            B(:,:) = B(:,:) .* ( (P(:,:).*(R(:,:).^(-1-beta/p)))*beta*V(:,:).' ./ ( 2*(R(:,:).^(-1))*V(:,:).' ) ).^(p/(beta+p));
            B(:,:) = max(B(:,:),1e-8);
            R(:,:) = B(:,:)*V(:,:);
            %%%%% Update V %%%%%
            V(:,:) = V(:,:) .* ( beta*B(:,:).'*(P(:,:).*(R(:,:).^(-1-beta/p))) ./ ( 2*B(:,:).'*(R(:,:).^(-1)) ) ).^(p/(beta+p));
            V(:,:) = max(V(:,:),1e-8);
            R(:,:) = B(:,:)*V(:,:);
        end

        for f=1:F
            RR = (squeeze(R(f,:)).^(beta/p)) .*  ((squeeze(garma(f,:)).^(2-beta)));
            temp = squeeze(xbar(f,:,:));
            R_l = temp*diag(1./max(RR,1e-8))*temp';
            r_l = temp*diag(1./max(RR,1e-8))*X_MTF(refMic,:,f)';
            G(:,f) = pinv(R_l)*r_l;
        end
    end
    % Inverse STFT for each source
    drb = ISTFT( xEarly, shiftSize, window, size(mix,1) );
else
    error('Date Error occurred');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%