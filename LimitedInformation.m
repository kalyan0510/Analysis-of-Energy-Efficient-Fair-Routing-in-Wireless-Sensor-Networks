M1 = [1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0; 0 1 0 0  0 1 0 0  0 1 0 0  0 0 0 0; 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
 M2 = [1 1 1 1 0 0 0 0  0 0 0 0  0 0 0 0 ;0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 1 1 1  0 0 0 0];
 M3 = [0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0;0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0;0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0];

 sz=3;
 const=12;
 n = (sz+1)*(sz+1)+4*(sz);
 R = [4;3;6];
 alf = [1/3;1/3;1/3];
 Con = 1e6;
 P = [0;4;2;3;1;0;3;2;4;1;0;6;0;0;0;0];
 E = [300;300;300];
 e = ones(3,1);
 fmin = ones(3,1)*(1e-10);

 A = [ (M1-M2) eye(3) zeros(3) zeros(3) zeros(3);
      (M1-M2) zeros(3) -eye(3) zeros(3) zeros(3);
      M2-M1-diag(alf)*M3 zeros(3) zeros(3)  eye(3) zeros(3);
      M3 zeros(3) zeros(3) zeros(3)  -eye(3)];
 b = [0;0;0;R;0;0;0;fmin];

 f = @(X)[(e'*M1) 0 0 0 0 0 0 0 0 0 0 0 0]*X;
 c = @(X) A*X-b;
 fx = @(y) [(e'*M1) 0 0 0 0 0 0 0 0 0 0 0 0]';
 cx = A';
 fxx = zeros(28,28);
 cxx = zeros(28,28);

 lam = ones(const,1);
 X = A\b;
 
 mu = 15;
 alpha = 0.5;
 X = X .+ 1e-9;
 z = mu ./ X';
 sX = diag(X);	
 Z = diag(z);
 e = ones(n,1);
 size(z)
 W = fxx + 1*cxx;

 for i = 1:25,
   


    AA = [W   cx        -eye(n);
        cx'  zeros(const,const)    zeros(const,n);
        Z    zeros(n,const) sX];
    
    bb = -[fx(X) + cx * lam - z';
        c(X);
        sX*Z*e-mu*e];
    
    dd = inv(AA) * bb;
    
    X = X + alpha * dd(1:28,1);
    lam = lam + alpha * dd(29:40,1);
    z = z + alpha * dd(41:68,1)';

    disp([f(X)])
    sX = diag(X);	
 	Z = diag(z);
     mu = mu / 10;
 end

 disp('we can see the value of f(x) decreases through the iterations')

 Ans = [X(1:4,1)' ; X(5:8,1)' ; X(9:12,1)' ; X(13:16,1)' ]