function[result1,result2,result3,result4,result5] = MY_PGdemo1(X,W,Ln,Ls,n_parameter,s_parameter,u_parameter,v_parameter)
options.u_parameter = u_parameter;
options.v_parameter = v_parameter;
options.k = 200;
options.maxiter = 1000;
[go,pro] = size(X);

[U,V] = My_reset_all(X,W,Ln,Ls,options);
Ln = Ln/n_parameter;
Ls = Ls/s_parameter;
for i=1:pro
    Ln(i,:)=Ln(i,:)/sum(Ln(i,:));
end

for i=1:go
    Ls(i,:)=Ls(i,:)/sum(Ls(i,:));
end
ppi = Ln;
gogo = Ls;
result1 = U*V;
result4 = gogo*(result1+X)+(result1+X)*ppi+result1;
result5 = gogo*(result1+X)+result1+(result1+X)*ppi+(result1+X)*ppi*ppi;
result2 = result1+(result1+X)*ppi;
result3 = gogo*(result1+X)+result1;



end
