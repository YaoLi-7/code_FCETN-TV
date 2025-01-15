function Xn=mytenmat(X,n)
S=size(X);N=numel(S);
if n==1
    Xn=reshape(X,S(1),prod(S(2:end)));
elseif n==N
    Xn=permute(reshape(X,prod(S(1:end-1)),S(N)),[2,1]);
else
    arr=[n 1:n-1 n+1:N];
    Xn=reshape(permute(X,arr),S(n),prod(S)/S(n));
end
end