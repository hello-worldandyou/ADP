%% ��С���� one shot 
clc;clear
Q = [1,0;0,2];
R = 1;
alpha = 0.9;
B = [0;1];
deltat = 0.001;

% ����ʽ��������ʼȨ�� 
W0 = zeros(5,1);
W(:,1)=W0;
% ����״̬�������ݼ�S1 S2 S3....  res
for j=1:10  
x1_=-0.05*j+(0.05*j-(-0.05*j))*rand(1,50); 
x2_=-0.05*j+(0.05*j-(-0.05*j))*rand(1,50); 
[m, n] = meshgrid(x1_, x2_');
[res(:,1),res(:,2)] = deal(reshape(m,[],1), reshape(n,[],1));
% ������Ҿ���
rowrank = randperm(size(res, 1)); % size���a��������randperm���Ҹ��е�˳��
res = res(rowrank,:);              % ����rowrank�������и��У�ע��rowrank��λ��

% ����ʽ������
%phi = [x(1);x(2);x(1)^2;x(2)^2;x(1)*x(2)];
%dphi'/dx = [1��0, 2*x(1), 0, x(2);                             %2X5ά��
%            0, 1, 0, 2*x(2), x(1)];
phi_X = [res(:,1),res(:,2),res(:,1).*res(:,1),res(:,2).*res(:,2),res(:,1).*res(:,2)];
PhiPhi = inv(phi_X'*phi_X);  % 5*5
for i = 1:2500
    %ƫF/X=[ƫfl/x1,ƫf1/x2��ƫf2/x1;ƫf2/x2;]       ֮ǰд�Ĳ�ͬ% A = [0, alpha*(-2*x(1))*x(2)-1; 1, alpha*(1-x(1)^2)]           %�ſ˱Ⱦ���
    A = [0, 1; alpha*(-2*res(i,1))*res(i,2)-1, alpha*(1-res(i,1)^2)];   
    dphiT = [1,0, 2*res(i,1), 0, res(i,2);0, 1, 0, 2*res(i,2), res(i,1)];
    u(i) = -inv(R)*B'*inv(A')*(dphiT*W(:,j)-Q*res(i,:)');
    % ��ɢ����ѧ����
    Xk1 = [res(i,1);res(i,2)]+([res(i,2);alpha*(1-res(i,1)^2)*res(i,2)-res(i,1)] + B*u(i))*deltat;     %  delta t  discrete
    phi_Xk1 = [Xk1(1),Xk1(2),Xk1(1)*Xk1(1),Xk1(2)*Xk1(2),Xk1(1)*Xk1(2)];
    RHS(i,1) = 0.5*(res(i,:)*Q*res(i,:)'+u(i)'*R*u(i))+phi_Xk1*W(:,j); 
end 
   W(:,j+1) = PhiPhi*phi_X'*RHS; 
end   


