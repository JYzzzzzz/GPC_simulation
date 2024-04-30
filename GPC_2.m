  %2019-03-14
  
%function��GPC����Ӧ�㷨
%parameter: 
%return: Theta_His: ���Ʋ����Ĺ켣
%note:����Ӧ�㷨����ȡСһ��Ч���ã���С�ȶ��Բ�
%     A��B��ֵ��ϵͳ��Ӱ��
function [y,ud_space,Theta_His] = GPC_2(pN,PA,PB,na,nb,N1,Nu,lambda,rho,soft_ele,Cn)
%############################## �û��༭�� ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N1);
%#######################################################################

%--------------------------------------------------2 ��ʾ������֪��
u=zeros(1, Last_t+1-First_t);    %������
ud_space=zeros(1, Last_t+1-First_t);   %���������Ĵ洢
y=zeros(1,Last_t+1-First_t);    %�����
w=rand(1,1+Last_t-First_t)*2-1;   %�Ŷ������ȷ�Χ��-1,1�����������
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %��ʾ��(k)/?

A = ones(1,na+1);   %A(1)��Զ����1
B = ones(1,nb+1);
n_theta = na+nb+1;                 %�����ȵĽ���
Theta = [-3;0;0.5;0];        %�ȵĳ�ֵ
  %�ȸ����ʾ��
Theta_His = ones(n_theta,1+Last_t-First_t)*1;
P_theta = ones(n_theta,n_theta)*2;
P_theta = P_theta+eye(n_theta)*1;  %�������������P
X = zeros(n_theta,1);          %΢��������[1~na,na+1~n_theta]

ud = 0;                        %���ڴ洢��u(t)
yd = 0;                        %���ڴ洢��y(t)
Ud = zeros(1,nb); Ud_p=1;  %�洢[��u(t-1),...,��u(t-nb)]
Yr = zeros(N1,1);          %�洢�ữ����趨ֵ��

%--------------------------------------------------3 ���ƹ���
tic;
for k = 1-First_t : 1+Last_t-First_t   %ʱ��ָ��t, t_p=1-First_t��Ӧʵ��ʱ��t=0
%for t = 0         :  Last_t
    %------------------------------3.1 ���t_pʱ�����
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %��ȡy(k)
    %if(k==1-First_t+50 || k==1-First_t+130)  %����˹�����
    %    y(k)=y(k)+10;
    %end
    %------------------------------3.2 ��������
    for i=n_theta:-1:2
         X(i) = X(i-1);  %����X
    end
    X(na+1) = ud;  %��u(t-1)����X
    X(1) = -yd;    %-��y(t-1)����X
    yd = y(k)-y(k-1); %update ��y(t)
    %-------------------- update ��(t),P(t-1)
    Mid = P_theta*X;         %P(t-1)*X(t)
    Mid = Mid/(X'*Mid+rho);  %P(t-1)*X(t)/(rho+X(t)'*P(t-1)*X(t))
    Theta = Theta+Mid*(yd-X'*Theta);  
    P_theta = (P_theta-Mid*X'*P_theta)/rho;
    Theta_His(:,k) = Theta;
    %------------------------------3.3 ʵʱ���¿�����
    for i=1:na
        A(i+1) = Theta(i);
    end
    for i=1:nb+1
        B(i) = Theta(i+na);
    end
    [P,Alpha,Beta]=GPC_getCtrlRule(A,na,B,nb,N1,Nu,lambda);
    %------------------------------3.4 ���¿�������u(t)
    ud = 0;    %��u(t)
    
    mid = yr(k)*(1-soft_ele);
    Yr(1) = y(k)*soft_ele + mid;
    for j=2:N1
        Yr(j) = Yr(j-1)*soft_ele +mid;
    end
    ud = ud + P*Yr;    %����P�Ǻ�����
    
    for j=1:na+1
        ud = ud - Alpha(j)*y(k+1-j);
    end
    
    for j=1:nb
        ud = ud - Beta(j)*Ud(mod(Ud_p-2+j,nb)+1);  %Ud_pָ���u(t-1),Ud_p+1ָ���u(t-2)
    end    %�������u(t)
    
    %----------2 ��u(k)
    u(k) = u(k-1) + ud;
    if u(k)<u_min           %��u(t)һЩ��������
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
    %---------- ���¡�u(k-1)�洢����
    Ud_p = Ud_p-1;
    if Ud_p < 1
        Ud_p = nb;
    end
    Ud(Ud_p) = ud;     %����5�С����¡�u(t-1)
end
run_time=toc;

%--------------------------------------------------4 ��ͼ
figure(4);
uicontrol('Style','text','Position',[10 0 150 20],'String',['����ʱ��: ',num2str(run_time),'s']);
        %����ʱ����ʾ
subplot(211);
plot(t,yr,'b--',t,y,'r-');
legend('�趨ֵyr','ʵ��ֵy','Location','Best');
axis([First_t Last_t -10 50]);
%xlim([First_t Last_t]);   %ֻ���ƺ����귶Χ
grid on;
title('����Ԥ������Ӧ�����㷨����');
ylabel('���y(k)');
subplot(212);
plot(t,ud_space);
%plot(t,Theta_His(1,:));
%axis([First_t Last_t -10 100]);
xlim([First_t Last_t]);
grid on;
ylabel('��������u_d(k)');
%ylabel('a_1(k)');
xlabel('k')

%--------------------------------------------------5 �༭����ֵ
%------------------------------ ������ʷ�켣
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
Theta_His = [Theta_His(:,1-First_t:1+Last_t-First_t)];

%------------------------------ ��ʷ�켣�����txt
DataTXT = fopen('.\GPC_simulate_data\GPC_2__data.txt','w'); %��txt�ļ�
fprintf(DataTXT,'%s\t%s\t%s\t\t%s\r\n','ʱ��k','���y','��������ud','������Ʋ�����(1~?)');
    %�������
data_len = length(y(:));  %���ݳ��ȣ�����ʱ�䣩��ȡ
theta_len = length(Theta_His(:,1));  %����������ȡ
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\t\t',y(j),ud_space(j));
    for i = 1:theta_len
        fprintf(DataTXT,'%.4f\t',Theta_His(i,j));
    end
    fprintf(DataTXT,'\r\n');
end

%------------------------------ ��ȡ��ӦPID���������ʱ�̣�
[E,F,G,H] = GPC_getEFGH(A,na,B,nb,N1,Nu);
Kp = 0;
Ki = 0;
Kd = 0;
for j=1:N1
    Kp = Kp - P(j)*(F(j,2)+2*F(j,3));
    Ki = Ki + P(j);
    Kd = Kd + P(j)*F(j,3);
end
fprintf('Kp = %.4f\n',Kp)
fprintf('Ki = %.4f\n',Ki)
fprintf('Kd = %.4f\n\n',Kd)

fclose(DataTXT);
