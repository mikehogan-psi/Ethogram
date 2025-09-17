function[X] = generate_X(b,R,T,template,P)
if size(b)
   [Ltrial,Nshape,Ntrial] = size(b);
else
   Nshape = 0;
   [Ltrial,Nd,Ntrial] = size(T);
   P = [];
end
[Np,Nd] = size(template);
X = NaN*zeros(Np,Nd,Ltrial,Ntrial);
for t = 1:Ntrial
    for l = 1:Ltrial
        X(:,:,l,t) = template;
        if size(b)            
            for s = 1:Nshape
                X(:,:,l,t) = X(:,:,l,t) + b(l,s,t)*P(:,:,s);
            end
        end
        X(:,:,l,t) = X(:,:,l,t)*R(:,:,l,t) + repmat(T(l,:,t),Np,1);
    end
end
