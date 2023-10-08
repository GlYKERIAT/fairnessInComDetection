k=2;# blue red
h=4;# poses dikaies omades thes

a=0.5;
b=0.009;
c=0.005;
d=0.001;

n=200;
block_sizes=(n/(k*h))*ones(1,k*h)  
sensitive=zeros(n,1);

#block_sizes=[[2,2],[2,2],[2,2]]
#block_sizes=[3*(n/10)*[0.5,0.3],3*(n/10)*[0.5,0.3],2*(n/10)*[0.5,0.3],1*(n/10)*[0.5,0.3],1*(n/10)*[0.5,0.3]]

 
    labels=zeros(n,1);
    for yyy=1:k
        for zzz=1:h
            block=(yyy-1)*h+zzz;
            sensitive(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=zzz;
            labels((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=yyy;
        end
    end
    
 function A=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes)
%generates
%
%INPUT:
%n ... number of elements
%a,b,c,d ... parameters / probabilities
%k ... number of clusters
%h ... number of groups
%block_sizes ... vector of length k*h with sum(block_sizes)=n; the 
%                first h entries correspond to the first cluster for 
%                the h groups, and so on ... 
%
%OUTPUT
%A ... adjacency matrix of size n x n


if (sum(block_sizes)~=n)||(length(block_sizes)~=(k*h))
    error('wrong input')
end

adja=random('Binomial',1,d,n,n);

for ell=1:k
    for mmm=1:k
        for ggg=1:h
            for fff=1:h
                if ell==mmm
                    if ggg==fff
                        adja((sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))),...
                            (sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))))=random('Binomial',1,a,block_sizes((ell-1)*h+(ggg)),block_sizes((ell-1)*h+(ggg)));
                    else
                        adja((sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))),...
                            (sum(block_sizes(1:((ell-1)*h+(fff-1))))+1):(sum(block_sizes(1:((ell-1)*h+(fff))))))=random('Binomial',1,c,block_sizes((ell-1)*h+(ggg)),block_sizes((ell-1)*h+(fff)));
                    end
                    
                else
                    if ggg==fff
                        adja((sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))),...
                                (sum(block_sizes(1:((mmm-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((mmm-1)*h+(ggg))))))=random('Binomial',1,b,block_sizes((ell-1)*h+(ggg)),block_sizes((mmm-1)*h+(ggg)));
                    end
                end
            end
        end
    end
end
            



A=triu(adja,1);
A=A+A';

end




adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);
adja;
sensitive
labels
save('C:\Users\Glykeria Toulina\Desktop/data4.txt', 'adja', '-ascii');
save('C:\Users\Glykeria Toulina\Desktop/labels4.txt', 'labels', '-ascii');
save('C:\Users\Glykeria Toulina\Desktop/sensitive4.txt', 'sensitive', '-ascii');