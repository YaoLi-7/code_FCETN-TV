function RecGroup=AW_LRTC(input,mask,lambda,eta,maxiter)
% output = X
imagesize=size(input);N=ndims(input); 
mask = logical(mask);
% RecGroup = zeros(imagesize);

    for n = 1 : N
        M{n} = zeros(imagesize);
        Lambda{n} = zeros(imagesize);
    end    
    RecGroup = input;     
    beta = 1;
%     omega=1/N*ones(1,N);
    omega=[1,1,0.1,10]/12.1;
    epsilon=1e-4;
    
    for it = 1: maxiter
        beta = beta * 1.05;
       %%%%%% to solve the M problem (auxiliary variable)      
        tau = omega/beta;
        for n = 1 : N
            temp2 = Unfold(M{n}, imagesize, n);
            w{n} = svd(temp2);
            weig{n}= w{n};
            w{n} = w{n} + epsilon;
            w{n} = 1./w{n};
            Z = Unfold(RecGroup + Lambda{n}/beta, imagesize, n);
            M{n} = Pro2TraceNorm_logdet(Z, tau(n)*w{n});
            M{n} = Fold(M{n}, imagesize, n);
        end
        
        
        %%%%%% to solve the X_tensor problem
        RecGroupOld = RecGroup;
        M_sum = 0; Lambda_sum = 0;
        for n = 1 : N
            M_sum = M_sum + M{n};
            Lambda_sum = Lambda_sum + Lambda{n};
        end
        RecGroup = (beta*M_sum - Lambda_sum)/(N*beta);
        RecGroup(mask) = RecGroupOld(mask);
               
        %%%%%% check the convergence
        res(it) = norm(RecGroup(~mask) - RecGroupOld(~mask)) / norm(RecGroupOld(:));
        fprintf('AW_LRTC: Iter %.0f, RSE %.6f \n',it,res(it));
%         if res(it) < 1e-8 && it>400
%             break;
%         end
        
%          if (mod(it,20)==0 || iter ==1)
%         difference=norm(T(:)-X(:))^2;
%         fprintf('AW_LRTC: Iter %.0f, RSE %.3f, diff %.3f, k %.1d, %.1d, %.1d, omega %.1d, %.1d, %.1d\n'...
%             ,iter,RSE(1),difference,k(1),k(2),k(N),omega(1),omega(2),omega(N));
%         end
        
        % omega update
        if it>50
        k=k_finder(weig,eta);
        omega=omega_updater(omega,k);
        end
        %%%%%% update the Lagrange multiplier
        for n = 1 : N
            Lambda{n} = Lambda{n} +  beta*(RecGroup - M{n}); 
        end
    end
end

    
    
% N=ndims(T); S=size(T);
% % initialization
% X=zeros(S); X(W==1)=T(W==1);X(W==0)=mean(T(:));
% omega=1/N*ones(1,N);
% M=cell(1,N); Y=M;Sigma=Y;
% for n=1:N
%     Y{n}=zeros(S);
%     M{n}=Y{n};
%     Sigma{n}=Y{n};
% end
% iter=0;
%  while iter~=maxiter
%     iter=iter+1;
%     
% %  ADMM update
% Msum=zeros(S);
% Ysum=zeros(S);
% for n=1:N
%     Df=mytenmat(X+(1/K)*Y{n},n);
%     [Shrink,~,Sig]=Pro2TraceNorm_hsi(Df,omega(n)/K);
%    % Shrink=Pro2TraceNorm_hsi(Df,omega(n)/K);
%     Sigma{n}=Sig;
%     M{n}=tenmat2tensor(Shrink,S,n);
%     Msum=Msum+M{n};
%     Ysum=Ysum+Y{n};
% end
% lastX=X;
% update_X=(Msum-Ysum/K)/N;
% X(W==0)=update_X(W==0);
% for n=1:N
%     Y{n}=Y{n}-K*(M{n}-X);
% end
% K=K*1.01;
% 
% % omega update
%  k=k_finder(Sigma,Q);
%  omega=omega_updater(omega,k);
% 
% % evaluations
%  if (mod(iter,20)==0 || iter ==1)
%         RSE=RSE_fun(T,X,W);
%         difference=norm(T(:)-X(:))^2;
%         fprintf('AW_LRTC: Iter %.0f, RSE %.3f, diff %.3f, k %.1d, %.1d, %.1d, omega %.1d, %.1d, %.1d\n'...
%             ,iter,RSE(1),difference,k(1),k(2),k(N),omega(1),omega(2),omega(N));
%  end
% end