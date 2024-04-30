  %2019-03-14
  
%function��������GPC
%parameter: 
function [y,ud_space,Theta_His] = GPC_d3(pN,PA,PB,na,nb,N,Nu,lambda,soft_ele,Cn)
%############################## �û��༭�� ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N);
%#######################################################################

%--------------------------------------------------1 ��ʾ��֪��
u=zeros(1, Last_t+1-First_t);    %������
ud_space=zeros(1, Last_t+1-First_t);   %���������Ĵ洢
y=zeros(1,Last_t+1-First_t);    %�����
w=rand(1,1+Last_t-First_t)*2-1;   %�Ŷ������ȷ�Χ��-1,1�����������
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %��ʾ��(k)/?

%��A��B��ȡG
%A=[1 -0.7 0];B=[0.9 -0.4];   %��һ�׶α��ض������Ϊ��
rho=1;

%----------------------------��ȡ����ʽP��Q
M=N+Nu+(na+1)+nb;   %�������
Theta = ones(M,1)*0.3;   %��
Theta_His = 0.3*ones(M,1+Last_t-First_t);  %�洢�ȵĹ켣������ѧϰ�о�
X = zeros(M,1);    %X'=[y(t),...,y(t-(N-1)),�ˡ�u(t+Nu-1-N),...,�ˡ�u(t-N)
                       %-y(t-N),...,-y(t-na-N),-��u(t-N-1),...,-��u(t-N-nb)];
Ud = zeros(N+nb,1);  %[��u(t-1),...,��u(t-N-nb)]T
ud = 0;
P_theta = ones(M,M);
P_theta = P_theta+eye(M);  %�������������P
X_taddN = zeros(M,1);   %X(t+N), �������u(t)ʱ���������˷�

tic;
%31-First_t
for k = 1-First_t : 1+Last_t-First_t   %ʱ��ָ��k, k=1-First_t��Ӧʵ��ʱ��t=0
%for t = 0         :  Last_t
    %----------------------------------------3.1 ���kʱ�����y(k)
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %��ȡy(k)
    %if(k==1-First_t+50 || k==1-First_t+130)  %����˹�����
        %y(k)=y(k)+10;
    %end
    %----------------------------------------3.2 ��������
    %------------------------------3.2.1 ����Ud,X
    Ud = [Ud(1:N+nb-1,1)];
    Ud = [ud;Ud];
    X = [X(1:M-1,1)];
    X = [y(k);X];
    X(N+1) = lambda*Ud(N+1-Nu);
    X(N+Nu+1) = -y(k-N);
    X(N+Nu+na+2) = -Ud(N+1);
    %-----------------------------3.2.2 ���¦�
    epsilon = Ud(N)-X'*Theta;
    Mid = P_theta*X;         %P(t-1)*X(t)
    Mid = Mid/(X'*Mid+rho);  %P(t-1)*X(t)/(rho+X(t)'*P(t-1)*X(t))
    Theta = Theta+Mid*epsilon;  
    P_theta = (P_theta-Mid*X'*P_theta)/rho;
    Theta_His(:,k) = Theta;  %��¼�ȹ켣
    %----------------------------------------3.3 ʵʱ���¿�����
    %�����u(t),u(t)
    for i=M:-1:N+Nu+2
        X_taddN(i) = X_taddN(i-1);
    end
    X_taddN(N+Nu+1) = -y(k);  %��ʷ�����
    X_taddN(N+Nu+na+2) = -Ud(1);  %��ʷ������
    %----- �趨ֵ��
    mid = yr(k)*(1-soft_ele);
    X_taddN(N)=y(k)*soft_ele + mid;
    for i=N-1:-1:1
        X_taddN(i) = X_taddN(i+1)*soft_ele + mid;
    end
    
    ud = X_taddN'*Theta;
    u(k) = u(k-1) + ud;
    if u(k)<u_min           %��u(t)һЩ��������
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
end
run_time=toc;

figure(5);
uicontrol('Style','text','Position',[10 0 150 20],'String',['����ʱ��: ',num2str(run_time),'s']);
        %����ʱ����ʾ

subplot(211);
plot(t,yr,'b--',t,y,'r-');
grid on;
title('�򵥵� ����Ԥ������Ӧ���� ֱ���㷨����');
%xlim([First_t Last_t]);
axis([First_t,Last_t,-10,60]);
ylabel('���y(k)');
legend('�趨ֵyr','ʵ��ֵy','Location','Best');

subplot(212);
plot(t,ud_space);
grid on;
xlim([First_t Last_t]);
%axis([First_t,Last_t,-40,40]);
ylabel('��������u_d(k)');
xlabel('k')

%--------------------------------------------------5 �༭����ֵ
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
Theta_His = [Theta_His(:,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\GPC_d3__data.txt','w');%��txt�ļ�
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
