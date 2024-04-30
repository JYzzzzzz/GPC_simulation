%2019-04-20 JYZ

%function: 将界面编辑框内的参数写入txt
function GPCmain_wt_editTXT(EditName,EditData)
EditTXT = fopen('.\GPCmain_reset.txt','w');%打开txt文件
EditLen = length(EditData(:));   %数据组数
for i=1:EditLen
   fprintf(EditTXT,'%s\t%.4f\r\n',EditName{i},EditData(i));
end
fclose(EditTXT);