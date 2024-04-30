  %2019-03-08
  
%function��������GPC
%parameter: pN: the number of plants
%           PA: a group of parameters of the group of "y"(the output)
%               Plant1:a1,a2,a3,... Plant2:a1,a2,a3,... Plant3:a1,a2,a3,...
%           PB: a group of parameters of the group of "u"(the input)
%               Plant1:b0,b1,b2,... Plant2:b0,b1,b2,... Plant3:b0,b1,b2,...
%           Cn: the range of noise.����Cn��
%return: y: ���
%        u: ������
function [y,ud_space] = GPC_1(pN,PA,PB,na,nb,N1,Nu,lambda,soft_ele,Cn)
%############################## �û��༭�� ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N1);
%#######################################################################

%--------------------------------------------------1 ��ʾ������֪��
u=zeros(1, Last_t+1-First_t);    %������
ud_space=zeros(1, Last_t+1-First_t);   %���������Ĵ洢
y=zeros(1,Last_t+1-First_t);    %�����
w=rand(1,1+Last_t-First_t)*2-1;   %�Ŷ������ȷ�Χ��-1,1�����������
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %��ʾ��(k)/?

A=[1,PA(1,:)];
B=PB(1,:);

%--------------------------------------------------2 ���߼��������
[P,Alpha,Beta] = GPC_getCtrlRule(A,na,B,nb,N1,Nu,lambda);

%--------------------------------------------------3 ���ƹ���
Ud = zeros(1,nb); Ud_p=1;  %�洢��u(t-1),...,��u(t-nb)
Yr = zeros(N1,1);          %�洢�ữ����趨ֵ��
tic;
for k = 1-First_t : 1+Last_t-First_t   %ʱ��ָ��t, t_p=1-First_t��Ӧʵ��ʱ��t=0
%for t = 0         :  Last_t
    %----------------------------------------3.1 ���kʱ�����
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %��ȡy(k)
    %----------------------------------------3.2 ���¿�����
    ud = 0;    %���ڴ洢��u(t)
    %----- �趨ֵ��
    mid = yr(k)*(1-soft_ele);
    Yr(1) = y(k)*soft_ele + mid;
    for j=2:N1
        Yr(j) = Yr(j-1)*soft_ele + mid;
    end
    ud = ud + P*Yr;    %����P�Ǻ�����
    %for j=1:N1
        %if(t_p+j>Last_t+1-First_t)
        %    ud = ud + P(j)*yr(Last_t+1-First_t);
        %else
        %    ud = ud + P(j)*yr(t_p+j);
        %end
        %ud = ud + P(j)*yr(k);
    %end
    %----- ��ʷ�����
    for j=1:na+1
        ud = ud - Alpha(j)*y(k+1-j);
    end
    %----- ��ʷ������
    for j=1:nb
        ud = ud - Beta(j)*Ud(mod(Ud_p-2+j,nb)+1);  %Ud_pָ���u(t-1),Ud_p+1ָ���u(t-2)
    end    %�������u(t)
    %---------- ���u(t)
    u(k) = u(k-1) + ud;
    if u(k)<u_min           %��u(t)һЩ��������
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
    %---------- ���¡�u(t-1)�洢����
    Ud_p = Ud_p-1;
    if Ud_p < 1
        Ud_p = nb;
    end
    Ud(Ud_p) = ud;     %����5�С����¡�u(t-1)
end
run_time=toc;

figure(3);
uicontrol('Style','text','Position',[10 0 150 20],'String',['����ʱ��: ',num2str(run_time),'s']);
        %����ʱ����ʾ
%set(gca,'position',[0.1 0.1 0.9 0.9]);
subplot(211);
plot(t,yr,'b--',t,y,'r-');
xlim([First_t Last_t]);
legend('�趨ֵyr','ʵ��ֵy','Location','Best');
grid on;
title('����Ԥ������㷨���棨����GPC��ģ����֪��');
ylabel('���y(k)');
subplot(212);
plot(t,ud_space);
grid on;
xlim([First_t Last_t]);
ylabel('��������u_d(k)');
xlabel('k')

%--------------------------------------------------5 �༭����ֵ
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\GPC_1__data.txt','w');%��txt�ļ�
fprintf(DataTXT,'%s\t%s\t%s\r\n','ʱ��k','���y','��������ud');
data_len = length(y(:));
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\r\n',y(j),ud_space(j));
end
fclose(DataTXT);

%�о��ʼǣ�
%1������lambda������ʱ��䳤������ƽ��
%2��NuԽС��������Խ��ϵͳ����Խ����Ӧʱ������ʱ
