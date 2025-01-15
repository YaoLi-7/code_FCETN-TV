function omega=omega_updater(omega,k)
N=numel(k);
ksum=sum(k);
denominator=0;
for n=1:N
    factor=exp(100*k(n)/ksum);
    denominator=denominator+factor;
end

for n=1:N
    omega(n)=exp(100*k(n)/ksum)/denominator;
end
end