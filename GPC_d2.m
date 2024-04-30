%2019-04-05

%function���Ľ���GPCֱ���㷨�����ڶ����Ծ��ӦǰNʱ�̹���ֵ��GPCֱ���㷨
%parameter: lambda������ʱȡ1.25
%           forget_ele���������ӣ�����ʱȡ0.95
%           soft_ele���ữ���ӣ�����ʱȡ0
function [y,ud_space,Theta_His] = GPC_d2(pN,PA,PB,na,nb,N,Nu,lambda,forget_ele,soft_ele,Cn)
%############################## �û��༭�� ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N);
M=0.01;    %�ɽ����������
C=2;       %���򻯷�������
%#######################################################################

%--------------------------------------------------1 ��ʾ������֪��
u=zeros(1, Last_t+1-First_t);    %������
ud_space=zeros(1, Last_t+1-First_t);   %���������Ĵ洢
y=zeros(1,Last_t+1-First_t);    %�����
w=rand(1,1+Last_t-First_t)*2-1;   %�Ŷ������ȷ�Χ��-1,1�����������
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %��ʾ��(k)/?
%------------------------------1.2 �������P��Q
close_ele=0.5;     %�����ⶨ��Ծ��ӦǰN��g0~gN-1��һ��������
                   %ʵ��Ӧ���У�g���ʹ�ò���ֵ
                   %��Ȼ������֤�ã�g�ʵ������ȡֵ��ϵͳҲ�������õĹ���
Itmd_val=close_ele*ones(N,1);   %��g0~gN-1
if(N>1)
    for i=2:N
        Itmd_val(i) = (1-Itmd_val(i-1))*close_ele+Itmd_val(i-1);
    end
end 
%Itmd_val=[0.40;1.06;1.45;1.48;1.30;1.13];
G=zeros(N,N);    %����G
for row=1:N
    i=0;
    while row+i<=N
        G(row+i,1+i)=Itmd_val(row);
        i=i+1;
    end
end
G = [G(:,1:Nu)];
Q = G'*G + lambda*eye(Nu);
Q = inv(Q);
P = Q * G';
P = [P(1,:)];
Q = [Q(1,:)];

%-------------------------------------------------2 �ռ俪��
Yr = zeros(N,1);       %�����ữ����Ӱ���yr����
Yd = zeros(na+N,1);    %�洢[��y(k),...,��y(k-(na+N-1))]
Ud = zeros(nb+N,1);    %�洢[��u(k-1),...,��u(k-N-nb)]
ud = 0;                %���ڴ洢��u(k)
Theta = ones(na+nb,1)*0;    %�ȵĳ�ֵ
Theta_His = ones(na+nb,1+Last_t-First_t)*0;
Y = zeros(na+N,1);  %���򻯹��������õ����洢[y(k),...,y(k-(na+N-1))]����ľ���ֵ����2
U = zeros(nb+N,1);  %�洢[u(k-1),...,u(k-nb-N)]����ľ���ֵ����2

%--------------------------------------------------3 ���ƹ���
tic;
for k = 1-First_t : 1+Last_t-First_t   %ʱ��ָ��t, t_p=1-First_t��Ӧʵ��ʱ��t=0
    %------------------------------3.1 ���t_pʱ�����
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %��ȡy(k)
    %if(k<71-First_t)
    %   y(k) = 0.7*y(k-1)+0.9*u(k-1)-0.4*u(k-2)+c*w(k);
    %else   %t=50ʱ���ض���ı�
    %   y(k) = -0.3*y(k-1)-0.7*y(k-2)+0.9*u(k-1)+0.7*u(k-2)+c*w(k);
    %end
    %------------------------------3.2 ��������
    %----------1 ����Yd��Yabs_mtpb_2��Yr
    Yd = [y(k)-y(k-1);[Yd(1:na+N-1,1)]];
    Y = [y(k);[Y(1:na+N-1,1)]];
    Itmd_val = yr(k)*(1-soft_ele);
    Yr(1) = y(k)*soft_ele + Itmd_val;
    for j=2:N
        Yr(j) = Yr(j-1)*soft_ele + Itmd_val;
    end
    %----------2 �����
    phi = 0;
    for i=1:Nu
        phi = phi + Ud(N+1-i)*Q(i);
    end
    phi = phi * lambda;        %��*Q*��u
    for i=N:-1:1
        phi = phi + P(i)*Y(N+1-i);
    end                        %P*y
    phi = phi - sum(P)*Y(N+1); %-sum(P)*y(k-N)
    %----------3 ��ȡX(k-N),X(k)
    X_mns_N = [[Yd(N+1:N+na,1)];[Ud(N+1:N+nb,1)]];  %na+nbά����
        %X(t-N)=[��y(k-N);...;��y(k-N-na+1);��u(t-1-N);...;��u(t-nb-N)]
    X = [[Yd(1:na,1)];[Ud(1:nb,1)]];
    %----------4 ��ȡ��������
    regu_ele = max(abs(Y))*2;    %ȫ�� regularization element ��������
    Itmd_val = max(abs(U))*2;    %Itmd_valȫ��intermediate value ��Ϊ�м���
    regu_ele = max([regu_ele,Itmd_val,C]);
    %----------5 ����
    phi = phi / regu_ele;
    X_mns_N = X_mns_N / regu_ele;
    %----------6 ���(k),  f(M,��(k))
    epsilon = phi - X_mns_N'*Theta - Ud(N)/regu_ele;
    if(epsilon>M)
        epsilon = epsilon - M;
    elseif(epsilon<-M)
        epsilon = epsilon + M;
    else
        epsilon = 0;
    end
    %----------7 ���¦�(k)
    Theta = Theta + forget_ele*epsilon*X_mns_N/(1+X_mns_N'*X_mns_N);
    Theta_His(:,k) = Theta;
    %----------8 �����u(k)
    ud = P*Yr - sum(P)*y(k) - X'*Theta;
    u(k) = u(k-1) + ud;  %���u(t)
    if u(k)<u_min           %��u(t)һЩ��������
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
    %----------9 ����Ud��U 
    Ud = [Ud(1:nb+N-1,1)];
    Ud = [ud;Ud];
    U = [U(1:nb+N-1,1)];
    U = [u(k);U];
end
run_time=toc;

figure(6);
uicontrol('Style','text','Position',[10 0 150 20],'String',['����ʱ��: ',num2str(run_time),'s']);
        %����ʱ����ʾ
subplot(211);
plot(t,yr,'b--',t,y,'r-');
xlim([First_t Last_t]);
%axis([First_t,Last_t,-10,60]);
legend('�趨ֵyr','ʵ��ֵy','Location','Best');
grid on;
title('���ڱ��ض���ǰN����Ծ��Ӧ��GPCֱ���㷨');
ylabel('���y(k)');
subplot(212);
plot(t,ud_space);
xlim([First_t Last_t]);
grid on;
ylabel('��������u_d(k)');
xlabel('k');

%--------------------------------------------------5 �༭����ֵ
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
Theta_His = [Theta_His(:,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\GPC_d2__data.txt','w');%��txt�ļ�
fprintf(DataTXT,'%s\t%s\t%s\t\t%s\r\n','ʱ��k','���y','��������ud','������Ʋ�����(1~?)');
data_len = length(y(:));
theta_len = length(Theta_His(:,1));
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\t\t',y(j),ud_space(j));
    for i = 1:theta_len
        fprintf(DataTXT,'%.4f\t',Theta_His(i,j));
    end
    fprintf(DataTXT,'\r\n');
end
fclose(DataTXT);
