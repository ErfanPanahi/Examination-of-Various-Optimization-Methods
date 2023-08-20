%% Problem 1

clc
clear
x1 = -10:0.01:10;
x2 = -10:0.01:10;
[X1,X2] = meshgrid(x1,x2);
f = F(X1,X2);
figure(1)
mesh(X1,X2,db(f));
title('db($f(x)$) in terms of $x_1$ and $x_2$)','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('db{$f(x)$}','Interpreter','latex');

%% Problem 2

figure(2)
subplot(1,2,1);
contour(X1,X2,f,"ShowText","on");
title('Contours of $f(x)$','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
subplot(1,2,2);
contour(X1,X2,db(f),"ShowText","on");
title('Contours of db($f(x)$)','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');

%% Problem 5

% mu = 0.1
mu1 = 0.1;
i1 = 0;
X_hat1 = [6;6];
X_i = X_hat1;
f1 = [];
while (F(X_hat1(1),X_hat1(2)) <= F(X_i(1),X_i(2)))
    X_i = X_hat1;
    f1 = [f1 F(X_hat1(1),X_hat1(2))];
    i1 = i1+1;
    X_hat1 = X_hat1 - mu1 * G(X_hat1(1),X_hat1(2));
end
fprintf("\nSteepest Descend Method \nMinimum value of f(x) is %f (x1 = %f ,x2 = %f) \n Needed iterantions (SD Method - mu = 0.1): %d \n",f1(end),X_i(1),X_i(2),i1);

% mu = 0.01
mu2 = 0.01;
i2 = 0;
X_hat2 = [6;6];
X_i = X_hat2;
f2 = [];
while (F(X_hat2(1),X_hat2(2)) <= F(X_i(1),X_i(2)))
    X_i = X_hat2;
    f2 = [f2 F(X_hat2(1),X_hat2(2))];
    i2 = i2+1;
    X_hat2 = X_hat2 - mu2 * G(X_hat2(1),X_hat2(2));
end
fprintf("\nSteepest Descend Method \nMinimum value of f(x) is %f (x1 = %f ,x2 = %f) \n Needed iterantions (SD Method - mu = 0.01): %d \n",f2(end),X_i(1),X_i(2),i2);

figure(5)
plot(1:i1,f1,'b',1:i2,f2,'r');
title('SD method - $\mu = 0.1$ and $\mu = 0.01$','Interpreter','latex');
legend('mu = 0.1','mu = 0.01');
xlabel('iteration','Interpreter','latex');
ylabel('$f(x)$','Interpreter','latex');

%% Problem 6

H = [2 1;
     1 2];
iN = 0;
X_hatN = [6;6];
X_i = X_hatN;
while (F(X_hatN(1),X_hatN(2)) <= F(X_i(1),X_i(2)))
    X_i = X_hatN;
    X_hatN = X_hatN - inv(H) * G(X_hatN(1),X_hatN(2));
    iN = iN+1;
end
fprintf("\nNewton Method \nMinimum value of f(x) is %f (x1 = %f ,x2 = %f) \n Needed iterantions (Newton Method): %d \n",F(X_i(1),X_i(2)),X_i(1),X_i(2),iN);

%% Problem 7

iAM = 0;
Xvec = [];
X_hatAM = [6;6];
X_i = X_hatAM;
while (F(X_hatAM(1),X_hatAM(2)) <= F(X_i(1),X_i(2)))
    % x2 is fixed
    X_i = X_hatAM;
    Xvec(:,end+1) = X_hatAM;
    X_hatAM(1) = (4-X_hatAM(2))/2;
    iAM = iAM+1;
    if(F(X_hatAM(1),X_hatAM(2)) > F(X_i(1),X_i(2)))
        break
    end
    % x1 is fixed
    X_i = X_hatAM;
    Xvec(:,end+1) = X_hatAM;
    X_hatAM(2) = (6-X_hatAM(1))/2;
    iAM = iAM+1;
end
fprintf("\nAlternation Minimization Method \nMinimum value of f(x) is %f (x1 = %f ,x2 = %f) \n Needed alterantions (AM Method): %d \n",F(X_i(1),X_i(2)),X_i(1),X_i(2),iAM);
figure(7)
hold on
contour(X1,X2,db(f));
plot(Xvec(1,:),Xvec(2,:));
title('Alternation Minimization Method','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
hold off

%% Problem 8

iGP = 0;
X_hatGP = [6;6];
X_i = X_hatGP;
while (F(X_hatGP(1),X_hatGP(2)) <= F(X_i(1),X_i(2)))
    X_i = X_hatGP;
    X_hatGP = X_hatGP - mu1 * G(X_hatGP(1),X_hatGP(2));
    X_hatGP = X_hatGP/(sqrt(X_hatGP(1)^2+X_hatGP(2)^2));
    iGP = iGP+1;
end
fprintf("\nGradient Projection Method \nMinimum value of f(x) is %f (x1 = %f ,x2 = %f) \n Needed iterantions (GP Method): %d \n",F(X_i(1),X_i(2)),X_i(1),X_i(2),iGP);

%% Functions 

function [f] = F (x1,x2)
    f = x1.^2 + x2.^2 - 4*x1 - 6*x2 + 13 + x1.*x2;  
end

function [g] = G (x1,x2)
    g = [2*x1-4+x2;
         2*x2-6+x1];
end
