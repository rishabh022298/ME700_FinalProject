% For theoretical calculation
nu = 0.3;
k = (3 - nu)/(1+nu);
mu = 80769.2;
h = 1;

%% Please import JInt.txt file first before running the script. Drag and drop works in MATLAB

Jint1 = table2array(JInt);

u_th = linspace(0,Jint1(length(Jint1(:,1)),1),100);
J_th = (2*mu*(k+1)/((k-1)*h)).*(u_th.^2);            % Theoretical Values

plot(u_th, J_th, '--')
hold on
plot(Jint1(:,1),Jint1(:,2),'r')

title("$J$-Integral values for different $q(r)$", Interpreter="latex")
xlabel("$u_0$", Interpreter="latex")
ylabel("$J$-Integral", Interpreter="latex")
legend("Theoretical","Numerical",Interpreter="latex",Location="northwest")
