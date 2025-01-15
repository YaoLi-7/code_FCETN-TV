function k=k_finder(Sigma,T)
N=numel(Sigma);
k=zeros(1,N);
for n=1:N
for i=1:numel(Sigma{n})
    rho=sum(Sigma{n}(1:i))/sum(Sigma{n});
    if rho>=T
        k(n)=i/numel(Sigma{n});
         break;
    end
end
end
end