%2019-04-17

%function������GPC����������ʱ�������༭�������ݵĳ�ʼ�����滻����
%parameter: retType: ����ֵ����ʽ��
%                    dat: EditData�е����ݣ���
%                    idx: EditData��Ӧ���±�
function retData = GPCmain_edits_find(EditName,EditData,str,retType)
EditLen = length(EditName(:));
retData = 0;
for i=1:EditLen  %˳�����
    if(strcmp(EditName(i),str)==1) %�ҵ���Ӧ��һ��
        if(strcmp(retType,'dat')==1) %����EditData�е�����
            retData = EditData(i);
        elseif(strcmp(retType,'idx')==1) %����EditData��Ӧ���±�
            retData = i;
        end
        return; %��ǰ��������
    end
end