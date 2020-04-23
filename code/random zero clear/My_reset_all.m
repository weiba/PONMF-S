function [U, V] = My_reset_all(X,W, Ln,Ls, options)
u_parameter = options.u_parameter;
v_parameter = options.v_parameter;
maxiter = options.maxiter;
k = options.k;
Dn = diag(sum(Ln,2));
Ds = diag(sum(Ls,2));

[m,n]=size(X);
U=abs(rand(m,k));
V=abs(rand(k,n));

niter=0;
obj1 = 100;
while niter<maxiter
    U_up = (W.*X)*V'+ Ls*U;
    U_down = (W.*(U*V))*V'+ Ds*U+u_parameter*U;
    U = U.*(U_up./(max(U_down,1e-10)));
    
    V_up = U'*(W.*X)+V*Ln;
    V_down = U'*(W.*(U*V))+V*Dn+v_parameter*V;
    V = V.*(V_up./max(V_down,1e-10));

    c3=W.*(X-U*V);
 	obj2 = sum(sum(c3.^2))+trace(U'*(Ds-Ls)*U)+trace(V*(Dn-Ln)*V')+u_parameter*sum(sum(U.^2))+v_parameter*sum(sum(V.^2));
	judge = abs(obj1-obj2)/obj2;
	niter = niter+1;
	if judge<10e-6
		break;
	end
	obj1 = obj2;
    	
end
n = niter
end
