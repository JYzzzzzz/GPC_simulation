
    %2019-03-02 

%function: ����GPC�����ʵ��м䲽�衪������Ԥ��ģ�Ͳ���EFGH��
          %���ڻ���GPC��ʼʱ�Ŀ����ʼ��㣬�Լ��������Ӧ�㷨
          %ÿ��tʱ�̲������ƺ���еĲ��裨���������ϲ�ѯԭ������
%parameter: A�����顣�洢���ض���Ķ���ʽ��ϵ���������Ǿ�ȷֵҲ�����ǹ���ֵ
           %na����ֵ��Ϊ��A�Ľ���-1��
           %B������A
           %nb����ֵ��Ϊ��B�Ľ���-1����B�Ľ����涨>=2
           %N1�����Ԥ��ʱ�򣬹涨��СԤ��ʱ��N0Ϊ1
           %Nu������ʱ��1<=Nu<=N1
%return: �����Ϊ���顣Ϊ�����ʵ�ϵ��
function [E,F,G,H] = GPC_getEFGH(A,na,B,nb,N1,Nu)
%--------------------------------------------------2 ���������
E=zeros(1,N1);
F=zeros(N1,na+1);
G=zeros(N1,Nu);
H=zeros(N1,nb);
%----------------------------------------2.1 ����E��F,G,H����
%------------------------------2.1.1 ��ֵ
E(1) = 1;                  %e0
Ad = conv(A,[1 -1]);
for i=1:na+1
   F(1+((i-1)*N1)) = -Ad(i+1); %f1,0  f1,1  f1,2
end
G(1)=E(1)*B(1);
for i=1:nb
   H(1+((i-1)*N1)) = B(i+1); 
end
%------------------------------2.1.2 ���ƹ���
if (N1>1)
    for j=2:N1         %���ƹ���
        E(j) = F(j-1+((1-1)*N1));
        for i=1:na
            F(j+((i-1)*N1)) = F(j-1+((i+1-1)*N1))-Ad(i+1)*E(j);
        end
        F(j+((na+1-1)*N1)) = -Ad(na+1+1)*E(j);
        G(j) = E(j)*B(1)+H(j-1+((1-1)*N1));
        if(nb>1)
            for i=1:nb-1
                H(j+((i-1)*N1)) = E(j)*B(i+1)+H(j-1+((i+1-1)*N1));
            end
        end
        H(j+((nb-1)*N1)) = E(j)*B(nb+1);
    end
end
for j=1:N1  %����G
    i=1;
    while j+i<=N1 && 1+i<=Nu
        G(j+i+((1+i-1)*N1)) = G(j+((1-1)*N1));
        i=i+1;
    end
end