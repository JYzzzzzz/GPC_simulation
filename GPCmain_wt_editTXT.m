%2019-04-20 JYZ

%function: ������༭���ڵĲ���д��txt
function GPCmain_wt_editTXT(EditName,EditData)
EditTXT = fopen('.\GPCmain_reset.txt','w');%��txt�ļ�
EditLen = length(EditData(:));   %��������
for i=1:EditLen
   fprintf(EditTXT,'%s\t%.4f\r\n',EditName{i},EditData(i));
end
fclose(EditTXT);