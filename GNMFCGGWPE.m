function [drb,B,V] = GNMFCGGWPE(mix, nb, L, delta, fftSize, shiftSize, ite, itNMF, beta, p, refMic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, window] = STFT( mix, fftSize, shiftSize, 'hamming' );

[F,T,M] = size(X);
X_MFT = permute(X,[3,1,2]);
X_MTF = permute(X,[3,2,1]);

B = max( rand( F, nb), eps );
V = max( rand( nb, T ), eps );
R(:,:) = B(:,:)*V(:,:);

myeps = 1e-8;

xbar(:,:,1:delta+L) = zeros(F,M*L,delta+L);
for t=delta+L+1:T
    for ail = 1:L
        xbar(:,(ail-1)*M+1:ail*M,t) = squeeze(X(:,t-delta-ail,:));
    end
end

fprintf('Iteration:      ');
for iter = 1:ite
    fprintf('\b\b\b\b\b%2d/%2d', iter,ite);
    if iter==1
        for f = 1:F
            desiredSignal(f,:) = X_MFT(refMic,f,:);
        end
        P = max(abs(desiredSignal).^beta,myeps);
        absDesiredSignal = abs(desiredSignal);
    else
        for f=1:F
            desiredSignal(f,:) = squeeze(X_MFT(refMic,f,:)).' - (G(:,f)'*(squeeze(xbar(f,:,:))));  % equation (3)
        end
        P = max(abs(desiredSignal).^beta,myeps);
        absDesiredSignal = abs(desiredSignal);
    end

    for i = 1:itNMF
        B(:,:) = B(:,:) .* ( (P(:,:).*(R(:,:).^(-1-beta/p)))*beta*V(:,:).' ./ ( 2*(R(:,:).^(-1))*V(:,:).' ) ).^(p/(beta+p)); % equation (36)
        B(:,:) = max(B(:,:),myeps);
        R(:,:) = B(:,:)*V(:,:);
        V(:,:) = V(:,:) .* ( beta*B(:,:).'*(P(:,:).*(R(:,:).^(-1-beta/p))) ./ ( 2*B(:,:).'*(R(:,:).^(-1)) ) ).^(p/(beta+p)); % equation (37)
        V(:,:) = max(V(:,:),myeps);
        R(:,:) = B(:,:)*V(:,:);
    end

    for f=1:F
        RR = (squeeze(R(f,:)).^(beta/p)) .*  ((squeeze(absDesiredSignal(f,:)).^(2-beta))); % equation (30)
        temp = squeeze(xbar(f,:,:));
        AutoCorrelationMatrix = temp*diag(1./max(RR,myeps))*temp'; % equation (28)
        RelationVector = temp*diag(1./max(RR,myeps))*X_MTF(refMic,:,f)'; % equation (29)
        G(:,f) = pinv(AutoCorrelationMatrix)*RelationVector; % equation (27)
    end
end
drb = ISTFT( desiredSignal, shiftSize, window, size(mix,1) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%