% % -------------------------------------- % %
% % ����Ӧѡ����ѧϰ�Ĳ�ֽ����㷨������   % %
% % -------------------------------------- % %
%% ����ͼ����
clear;
clc;
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30; %���д���
D=30;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "Algs";  % ���ڱ���Ա��㷨��־���ַ���

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
   
DODEMatrix=zeros(runNumber,Max_FES);
OMLDEMatrix=zeros(runNumber,Max_FES);
OBAMatrix=zeros(runNumber,Max_FES);
NBOLDEMatrix=zeros(runNumber,Max_FES);
GPODEMatrix=zeros(runNumber,Max_FES);
JaDEMatrix=zeros(runNumber,Max_FES);
ACDE_FMatrix=zeros(runNumber,Max_FES);

Algs_FES = zeros(7,Max_FES);
result_gbest=zeros(runNumber,7);
for k=1:30
    func_num=k;
    if k==2
       continue;
    end
    fprintf("\n----------------------------------\n");
    fprintf("��ʼ���ԶԱ��㷨��%dά-F%d���� >>>>\n",D,k);
    fprintf("----------------------------------\n");
    
for i=1:runNumber
    fprintf("-----��%d������-----\n",i);
    fprintf("OMLDE--->");
    [Pb,~,FEs_fitness]=OMLDE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    OMLDEMatrix(i,:)=FEs_fitness;
    result_gbest(i,1)=Pb;
    
    fprintf("DODE--->");
    [Pb,~,FEs_fitness]=DODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    DODEMatrix(i,:)=FEs_fitness;
    result_gbest(i,2)=Pb;
    
    fprintf("OBA--->");
    [Pb,~,FEs_fitness]=OBA(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    OBAMatrix(i,:)=FEs_fitness;
    result_gbest(i,3)=Pb;    
    
    fprintf("NBOLDE--->");
    [Pb,~,FEs_fitness]=NBOLDE(func_num,fhd,D,NP,CR,gen_max,Max_FES,border,func_num);
    NBOLDEMatrix(i,:)=FEs_fitness;
    result_gbest(i,4)=Pb;
    
    fprintf("GPODE--->\n");
    [Pb,~,FEs_fitness]=GPODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    GPODEMatrix(i,:)=FEs_fitness;
    result_gbest(i,5)=Pb;
   
    fprintf("JaDE--->");
    [Pb,~,FEs_fitness]=JaDE(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    JaDEMatrix(i,:)=FEs_fitness;
    result_gbest(i,6)=Pb;
    
    fprintf("ACDE_F\n");
    [Pb,~,FEs_fitness]=ACDE_F(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    ACDE_FMatrix(i,:)=FEs_fitness;
    result_gbest(i,7)=Pb;
    
    
end
Algs_FES(1,:)= mean(OMLDEMatrix,1);
Algs_FES(2,:)= mean(DODEMatrix,1);
Algs_FES(3,:)= mean(OBAMatrix,1);
Algs_FES(4,:)= mean(NBOLDEMatrix,1);
Algs_FES(5,:) = mean(GPODEMatrix,1);
Algs_FES(6,:) = mean(JaDEMatrix,1);
Algs_FES(7,:) = mean(ACDE_FMatrix,1);

if D==30
    path='E:\MATLAB\project\ASODE\save_data\boxplotData\30D\';
elseif D==50
    path='E:\MATLAB\project\ASODE\save_data\boxplotData\50D\';
end
% filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
% save([path,filename],'Algs_FES');
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
save([path,filename],'result_gbest');

fprintf("\nOMLDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,1)),std(result_gbest(:,1)));
fprintf("\nDODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,2)),std(result_gbest(:,2)));
fprintf("\nOBA:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,3)),std(result_gbest(:,3)));
fprintf("\nNBOLDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,4)),std(result_gbest(:,4)));
fprintf("\nGPODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,5)),std(result_gbest(:,5)));
fprintf("\nJaDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,6)),std(result_gbest(:,6)));
fprintf("\nACDE_F:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,7)),std(result_gbest(:,7)));

%����Ա��㷨�ľ�ֵ�ͷ��excel
% save_data(k,func_num,s,D,str);

end  % ���Ժ���ѭ������

%% �Ա��㷨�ۺϲ���
clear;
clc;
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30; %���д���
D=10;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "Algs";  % ���ڱ���Ա��㷨��־���ַ���

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
   
DODEMatrix=zeros(runNumber,Max_FES);
OMLDEMatrix=zeros(runNumber,Max_FES);
OBAMatrix=zeros(runNumber,Max_FES);
NBOLDEMatrix=zeros(runNumber,Max_FES);
GPODEMatrix=zeros(runNumber,Max_FES);
JaDEMatrix=zeros(runNumber,Max_FES);
ACDE_FMatrix=zeros(runNumber,Max_FES);

Algs_FES = zeros(7,Max_FES);
result_gbest=zeros(runNumber,7);
for k=1:30
    func_num=k;
    if k==2
       continue;
    end
    fprintf("\n----------------------------------\n");
    fprintf("��ʼ���ԶԱ��㷨��%dά-F%d���� >>>>\n",D,k);
    fprintf("----------------------------------\n");
    
for i=1:runNumber
    fprintf("-----��%d������-----\n",i);
    fprintf("OMLDE--->");
    [Pb,~,FEs_fitness]=OMLDE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    OMLDEMatrix(i,:)=FEs_fitness;
    result_gbest(i,1)=Pb;
    
    fprintf("DODE--->");
    [Pb,~,FEs_fitness]=DODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    DODEMatrix(i,:)=FEs_fitness;
    result_gbest(i,2)=Pb;
    
    fprintf("OBA--->");
    [Pb,~,FEs_fitness]=OBA(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    OBAMatrix(i,:)=FEs_fitness;
    result_gbest(i,3)=Pb;    
    
    fprintf("NBOLDE--->");
    [Pb,~,FEs_fitness]=NBOLDE(func_num,fhd,D,NP,CR,gen_max,Max_FES,border,func_num);
    NBOLDEMatrix(i,:)=FEs_fitness;
    result_gbest(i,4)=Pb;
    
    fprintf("GPODE--->\n");
    [Pb,~,FEs_fitness]=GPODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    GPODEMatrix(i,:)=FEs_fitness;
    result_gbest(i,5)=Pb;
   
    fprintf("JaDE--->");
    [Pb,~,FEs_fitness]=JaDE(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    JaDEMatrix(i,:)=FEs_fitness;
    result_gbest(i,6)=Pb;
    
    fprintf("ACDE_F\n");
    [Pb,~,FEs_fitness]=ACDE_F(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    ACDE_FMatrix(i,:)=FEs_fitness;
    result_gbest(i,7)=Pb;
end
Algs_FES(1,:)= mean(OMLDEMatrix,1);
Algs_FES(2,:)= mean(DODEMatrix,1);
Algs_FES(3,:)= mean(OBAMatrix,1);
Algs_FES(4,:)= mean(NBOLDEMatrix,1);
Algs_FES(5,:) = mean(GPODEMatrix,1);
Algs_FES(6,:) = mean(JaDEMatrix,1);
Algs_FES(7,:) = mean(ACDE_FMatrix,1);


if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
end
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
save([path,filename],'Algs_FES');

fprintf("\nOMLDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,1)),std(result_gbest(:,1)));
fprintf("\nDODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,2)),std(result_gbest(:,2)));
fprintf("\nOBA:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,3)),std(result_gbest(:,3)));
fprintf("\nNBOLDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,4)),std(result_gbest(:,4)));
fprintf("\nGPODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,5)),std(result_gbest(:,5)));
fprintf("\nJaDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,6)),std(result_gbest(:,6)));
fprintf("\nACDE_F:\n��ֵΪ:%d\nstdΪ:%d\n",mean(result_gbest(:,7)),std(result_gbest(:,7)));

%����Ա��㷨�ľ�ֵ�ͷ��excel
save_data(k,func_num,s,D,str);

end  % ���Ժ���ѭ������


%% ��������ASODE
clear;
clc; 
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30; %���д���
D=50;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D;
gen_max = Max_FES / NP;  %���������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
%func_num=1;
fhd=str2func('cec17_func');
str = "ASODE"; % ���ڱ�־ASODE

ASODEMatrix=zeros(runNumber,Max_FES);

s=zeros(1,runNumber);
ASODE_FES = zeros(1,Max_FES);

for k=26:26
    func_num=k;
    if k==2
       continue;
    end
    fprintf("--------------------------\n");
    fprintf("��ʼ����%s��%dά-F%d���� >>>>\n",str,D,k);
    fprintf("--------------------------\n");
    
    %%% ��ȡ�����Ժ���Ӧ����������  ����ASODEʹ��%%%
    if D==10 
        RD_path='E:\MATLAB\project\ASODE\save_data\RANK_DIV_10\';
        RF_path='E:\MATLAB\project\ASODE\save_data\RANK_FIT_10\';
    elseif D==30
        RD_path='E:\MATLAB\project\ASODE\save_data\RANK_DIV_30\';
        RF_path='E:\MATLAB\project\ASODE\save_data\RANK_FIT_30\';
    elseif D==50
        RD_path='E:\MATLAB\project\ASODE\save_data\RANK_DIV_50\';
        RF_path='E:\MATLAB\project\ASODE\save_data\RANK_FIT_50\';
    end
    RD_filename=strcat('RD','_',int2str(D),'D_F',int2str(k),'.mat');
    load([RD_path,RD_filename]); 
    RF_filename=strcat('RF','_',int2str(D),'D_F',int2str(k),'.mat');
    load([RF_path,RF_filename]); 
    
for i=1:runNumber
    fprintf("------ASODE��%d������------\n",i);
    [Pb,~,FEs_fitness]=ASODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,RANK_DIV,RANK_FIT,func_num);
    ASODEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
end
ASODE_FES(1,:)= mean(ASODEMatrix,1);

if D==10
    path='E:\MATLAB\project\ASODE\save_data\ASODE_data_10\ASODE_10\' ;
elseif D==30
    path='E:\MATLAB\project\ASODE\save_data\ASODE_data_30\ASODE_30\';
    path_box = 'E:\MATLAB\project\ASODE\save_data\boxplotData\30D\ASODE\';
elseif D==50
    path='E:\MATLAB\project\ASODE\save_data\ASODE_data_50\ASODE_50\';
    path_box = 'E:\MATLAB\project\ASODE\save_data\boxplotData\50D\ASODE\';
end
filename=strcat('ASODE','_',int2str(D),'D_F',int2str(k),'.mat');
% save([path,filename],'ASODE_FES');
save([path_box,filename],'s');
% 
fprintf("\nASODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));


% %%% ����ASODE���ݵ�excel�ļ� %%%
% save_data(k,func_num,s,D,str); 
% 
% 
% % ���������㷨���ݣ���������ͼ
% if D==10
%     path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
% elseif D==30
%     path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
% elseif D==50
%     path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
% end
% filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
% load([path,filename]);
% 
% figure(k);
% xx = 1:Max_FES;
% plot(xx,log(Algs_FES(1,xx)),'r-*',...
%      xx,log(Algs_FES(2,xx)),'k-o',...
%      xx,log(Algs_FES(3,xx)), 'g-^',...
%      xx,log(Algs_FES(4,xx)), 'b-d',...
%      xx,log(Algs_FES(5,xx)), 'm-s',...
%      xx,log(ASODE_FES(1,xx)), 'c-p',...
%      'MarkerIndices',1:(Max_FES/10):Max_FES,'LineWidth',1);
% legend('OMLDE','DODE','OBA','NBOLDE','GODE','ASODE');
% xlabel('Number of Function Evaluations');
% ylabel('Average Function Values(Log)');

end

%% ��������OMLDE
clear;
clc;
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=10; %���д���
D=10;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "OMLDE";  % ���ڱ���Ա��㷨��־���ַ���

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
   
OMLDEMatrix=zeros(runNumber,Max_FES);

Local_Algs_FES = zeros(1,Max_FES);
s=zeros(1,runNumber);
for k=5:5
    func_num=k;
    if k==2
       continue;
    end
    fprintf("\n----------------------------------\n");
    fprintf("��ʼ����%s��%dά-F%d���� >>>>\n",str,D,k);
    fprintf("----------------------------------\n");
    
for i=1:runNumber
    fprintf("-----%s��%d������-----\n",str,i);
    [Pb,~,FEs_fitness]=OMLDE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    OMLDEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
     
end
Local_Algs_FES(1,:)= mean(OMLDEMatrix,1);

if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
end
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
load([path,filename]);
Algs_FES(1,:)=Local_Algs_FES(1,:);
save([path,filename],'Algs_FES');

fprintf("\nOMLDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));

%����OMLDE��Best,Worst,Median,Mean,Std��excel
save_data(k,func_num,s,D,str);
end

%% ��������DODE
clear;
clc;
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30; %���д���
D=10;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "DODE";  % ���ڱ���Ա��㷨��־���ַ���

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
   
OMLDEMatrix=zeros(runNumber,Max_FES);

Local_Algs_FES = zeros(1,Max_FES);
s=zeros(1,runNumber);
for k=5:5
    func_num=k;
    if k==2
       continue;
    end
    fprintf("\n----------------------------------\n");
    fprintf("��ʼ����%s��%dά-F%d���� >>>>\n",str,D,k);
    fprintf("----------------------------------\n");
    
for i=1:runNumber
    fprintf("-----%s��%d������-----\n",str,i);
    [Pb,~,FEs_fitness]=DODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    OMLDEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
     
end
Local_Algs_FES(1,:)= mean(OMLDEMatrix,1);

if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
end
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
load([path,filename]);
Algs_FES(2,:)=Local_Algs_FES(1,:);
save([path,filename],'Algs_FES');

fprintf("\nDODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));

%����DODE��Best,Worst,Median,Mean,Std��excel
save_data(k,func_num,s,D,str);
end

%% ��������OBA
clear;
clc;
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30; %���д���
D=10;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "OBA";  % ���ڱ���Ա��㷨��־���ַ���

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
   
OMLDEMatrix=zeros(runNumber,Max_FES);

Local_Algs_FES = zeros(1,Max_FES);
s=zeros(1,runNumber);
for k=1:1
    func_num=k;
    if k==2
       continue;
    end
    fprintf("\n----------------------------------\n");
    fprintf("��ʼ����%s��%dά-F%d���� >>>>\n",str,D,k);
    fprintf("----------------------------------\n");
    
for i=1:runNumber
    fprintf("-----%s��%d������-----\n",str,i);
    [Pb,~,FEs_fitness]=OBA(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    OMLDEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
     
end
Local_Algs_FES(1,:)= mean(OMLDEMatrix,1);

if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
end
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
load([path,filename]);
Algs_FES(3,:)=Local_Algs_FES(1,:);
save([path,filename],'Algs_FES');

fprintf("\nOBA:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));

%����OBA��Best,Worst,Median,Mean,Std��excel
save_data(k,func_num,s,D,str);
end

%% ��������NBOLDE
clear;
clc;
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30; %���д���
D=10;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "NBOLDE";  % ���ڱ���Ա��㷨��־���ַ���

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
   
OMLDEMatrix=zeros(runNumber,Max_FES);

Local_Algs_FES = zeros(1,Max_FES);
s=zeros(1,runNumber);
for k=1:1
    func_num=k;
    if k==2
       continue;
    end
    fprintf("\n----------------------------------\n");
    fprintf("��ʼ����%s��%dά-F%d���� >>>>\n",str,D,k);
    fprintf("----------------------------------\n");
    
for i=1:runNumber
    fprintf("-----%s��%d������-----\n",str,i);
    [Pb,~,FEs_fitness]=NBOLDE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    OMLDEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
     
end
Local_Algs_FES(1,:)= mean(OMLDEMatrix,1);

if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
end
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
load([path,filename]);
Algs_FES(4,:)=Local_Algs_FES(1,:);
save([path,filename],'Algs_FES');

fprintf("\nNBOLDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));

%����NBOLDE��Best,Worst,Median,Mean,Std��excel
save_data(k,func_num,s,D,str);
end

%% ��������GPODE
clear;
clc;
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30; %���д���
D=10;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "GPODE";  % ���ڱ���Ա��㷨��־���ַ���

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
   
GPODEMatrix=zeros(runNumber,Max_FES);

Local_Algs_FES = zeros(1,Max_FES);
s=zeros(1,runNumber);
for k=1:1
    func_num=k;
    if k==2
       continue;
    end
    fprintf("\n----------------------------------\n");
    fprintf("��ʼ����%s��%dά-F%d���� >>>>\n",str,D,k);
    fprintf("----------------------------------\n");
    
for i=1:runNumber
    fprintf("-----%s��%d������-----\n",str,i);
    [Pb,~,FEs_fitness]=GPODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    GPODEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
     
end
Local_Algs_FES(1,:)= mean(GPODEMatrix,1);

if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
end
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
load([path,filename]);
Algs_FES(5,:)=Local_Algs_FES(1,:);
save([path,filename],'Algs_FES');

fprintf("\nGPODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));

%����GPODE��Best,Worst,Median,Mean,Std��excel
save_data(k,func_num,s,D,str);
end

%% ��ӶԱ��㷨ACDE_F��JaDE��������
clear;
clc; 
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=20;%���д���
D=50;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D;
gen_max = Max_FES / NP;  %���������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
%func_num=1;
fhd=str2func('cec17_func');

global fbias
%����ֵƫ����%
fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];

JaDEMatrix=zeros(runNumber,Max_FES);
ACDE_FMatrix=zeros(runNumber,Max_FES);

s=zeros(2,runNumber);
Local_Algs_FES = zeros(2,Max_FES);

for k=1:30
    func_num=k;
    if k==2
       continue;
    end
    fprintf("----------------------------------\n");
    fprintf("��ʼ���ԶԱ��㷨ACDE_F��SaDE��%dά-F%d����\n",D,k);
    fprintf("----------------------------------\n");
for i=1:runNumber
    fprintf("-----��%d������-----\n",i);
    fprintf("JaDE--->");
    [Pb,~,FEs_fitness]=JaDE(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    JaDEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
    
    fprintf("ACDE_F\n");
    [Pb,~,FEs_fitness]=ACDE_F(func_num,fhd,D,NP,gen_max,Max_FES,border,func_num);
    ACDE_FMatrix(i,:)=FEs_fitness;
    s(2,i)=Pb;
end

Local_Algs_FES(1,:) = mean(JaDEMatrix,1);
Local_Algs_FES(2,:) = mean(ACDE_FMatrix,1);

if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\';
end
filename=strcat('Algs','_',int2str(D),'D_F',int2str(k),'.mat');
load([path,filename]);
Algs_FES(6,:)=Local_Algs_FES(1,:);
Algs_FES(7,:)=Local_Algs_FES(2,:);
save([path,filename],'Algs_FES');

fprintf("\nJaDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));
fprintf("\nACDE_F:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(2,:)),std(s(2,:)));
Algs_Save_datas(k,D,s);

end

%% OBL�Ա��㷨���ݱ������
clear;
clc; 
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=20;%���д���
D=50;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D; % �����������
gen_max = Max_FES/NP;  % �������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');

%��¼ÿ����������Ӧֵ
ODEMatrix=zeros(runNumber,Max_FES);
QODEMatrix=zeros(runNumber,Max_FES);
QRODEMatrix=zeros(runNumber,Max_FES);
EODEMatrix=zeros(runNumber,Max_FES);
REODEMatrix=zeros(runNumber,Max_FES);
FQRODEMatrix=zeros(runNumber,Max_FES);
GODEMatrix=zeros(runNumber,Max_FES);
RODEMatrix=zeros(runNumber,Max_FES);
CODEMatrix=zeros(runNumber,Max_FES);
PODEMatrix=zeros(runNumber,Max_FES);
SODEMatrix=zeros(runNumber,Max_FES);
OCDEMatrix=zeros(runNumber,Max_FES);
COODEMatrix=zeros(runNumber,Max_FES);


Algs_FES = zeros(13,Max_FES);
s=zeros(13,runNumber);
for k=22:30
    func_num=k;
    if k==2
       continue;
    end
    fprintf("-----------------------------------\n");
    fprintf("OBL���Գؿ�ʼ����%dά-F%d���� >> >>\n",D,k);
    fprintf("-----------------------------------\n");
for i=1:runNumber
    fprintf("-----��%d������-----\n",i);
    fprintf("ODE--->");
    [Pb,~,~,FEs_fitness]=ODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    ODEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
    
    fprintf("QODE--->");
    [Pb,~,~,FEs_fitness]=QODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    QODEMatrix(i,:)=FEs_fitness;
    s(2,i)=Pb;
    
    fprintf("QRODE--->");
    [Pb,~,~,FEs_fitness]=QRODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    QRODEMatrix(i,:)=FEs_fitness;
    s(3,i)=Pb;   
    
    fprintf("EODE--->");
    [Pb,~,~,FEs_fitness]=EODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    EODEMatrix(i,:)=FEs_fitness;
    s(4,i)=Pb;
    
    fprintf("REODE--->");
    [Pb,~,~,FEs_fitness]=REODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    REODEMatrix(i,:)=FEs_fitness;
    s(5,i)=Pb;
    
    fprintf("FQRODE--->");
    [Pb,~,~,FEs_fitness]=FQRODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    FQRODEMatrix(i,:)=FEs_fitness;
    s(6,i)=Pb;
    
    fprintf("GODE--->");
    [Pb,~,~,FEs_fitness]=GODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    GODEMatrix(i,:)=FEs_fitness;
    s(7,i)=Pb; 
    
    fprintf("RODE--->");
    [Pb,~,~,FEs_fitness]=RODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    RODEMatrix(i,:)=FEs_fitness;
    s(8,i)=Pb;
    
    fprintf("CODE--->");
    [Pb,~,~,FEs_fitness]=CODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    CODEMatrix(i,:)=FEs_fitness;
    s(9,i)=Pb;
    
    fprintf("PODE--->");
    [Pb,~,~,FEs_fitness]=PODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    PODEMatrix(i,:)=FEs_fitness;
    s(10,i)=Pb;
    
    fprintf("SODE--->");
    [Pb,~,~,FEs_fitness]=SODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    SODEMatrix(i,:)=FEs_fitness;
    s(11,i)=Pb;
    
    fprintf("OCDE--->");
    [Pb,~,~,FEs_fitness]=OCDE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    OCDEMatrix(i,:)=FEs_fitness;
    s(12,i)=Pb;
    
    fprintf("COODE\n");
    [Pb,~,~,FEs_fitness]=COODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,func_num);
    COODEMatrix(i,:)=FEs_fitness;
    s(13,i)=Pb;   

end

Algs_FES(1,:)= mean(ODEMatrix,1);Algs_FES(2,:)= mean(QODEMatrix,1);
Algs_FES(3,:)= mean(QRODEMatrix,1);Algs_FES(4,:)= mean(EODEMatrix,1);
Algs_FES(5,:)= mean(REODEMatrix,1);Algs_FES(6,:)= mean(FQRODEMatrix,1);
Algs_FES(7,:)= mean(GODEMatrix,1);Algs_FES(8,:)= mean(RODEMatrix,1);
Algs_FES(9,:)= mean(CODEMatrix,1);Algs_FES(10,:)= mean(PODEMatrix,1);
Algs_FES(11,:)= mean(SODEMatrix,1);Algs_FES(12,:)= mean(OCDEMatrix,1);
Algs_FES(13,:)= mean(COODEMatrix,1);

fprintf("ODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));
fprintf("\nQODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(2,:)),std(s(2,:)));
fprintf("\nQRODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(3,:)),std(s(3,:)));
fprintf("\nEODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(4,:)),std(s(4,:)));
fprintf("\nREODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(5,:)),std(s(5,:)));
fprintf("\nFQRODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(6,:)),std(s(6,:)));
fprintf("\nGODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(7,:)),std(s(7,:)));
fprintf("\nRODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(8,:)),std(s(8,:)));
fprintf("\nCODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(9,:)),std(s(9,:)));
fprintf("\nPODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(10,:)),std(s(10,:)));
fprintf("\nSODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(11,:)),std(s(11,:)));
fprintf("\nOCDE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(12,:)),std(s(12,:)));
fprintf("\nCOODE:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(13,:)),std(s(13,:)));
Algs_Save_datas(k,D,s);

%%% ��������Ժ���Ӧ�ȵ����� %%%
if D==10 
    FES_path='E:\MATLAB\project\DE\save_data\OBLS_data_10\';
elseif D==30
    FES_path='E:\MATLAB\project\DE\save_data\OBLS_data_30\';
elseif D==50
    FES_path='E:\MATLAB\project\DE\save_data\OBLS_data_50\';
end
FES_filename=strcat('OBLS','_',int2str(D),'D_F',int2str(k),'.mat');
save([FES_path,FES_filename],'Algs_FES');  %��¼����ֵ��������ͼ
fprintf("\nOBL���Գ����ݱ���ɹ�����\n\n");

end

%% ASODE���ҶԱȲ���
clear;
clc; 
close all; 
addpath(genpath(pwd));  %��һ������matlab����·���ϵ��ļ�
runNumber=30;%���д���
D=50;         %ά��
NP=100;       %NPΪ��Ⱥ��ģ
F=0.5;        %ͻ������
CR=0.9;       %�������
Max_FES = 10000 * D;
gen_max = Max_FES / NP;  %���������� 
border=100;     %���½����ֵ��һ��Ϊ�Գ������ռ䣩
fhd=str2func('cec17_func');
str = "ASODE"; % ���ڱ�־ASODE

ASODE_1Matrix=zeros(runNumber,Max_FES);
ASODE_2Matrix=zeros(runNumber,Max_FES);
ASODE_3Matrix=zeros(runNumber,Max_FES);
ASODE_4Matrix=zeros(runNumber,Max_FES);

s=zeros(4,runNumber);
msASODE_FES = zeros(4,Max_FES);

for k=20:30
    func_num=k;
    if k==2
       continue;
    end
    fprintf("--------------------------\n");
    fprintf("��ʼ����%s��%dά-F%d���� >>>>\n",str,D,k);
    fprintf("--------------------------\n");
    
    %%% ��ȡ�����Ժ���Ӧ����������  ����ASODEʹ��%%%
    if D==10 
        RD_path='E:\MATLAB\project\DE\save_data\RANK_DIV_10\';
        RF_path='E:\MATLAB\project\DE\save_data\RANK_FIT_10\';
    elseif D==30
        RD_path='E:\MATLAB\project\DE\save_data\RANK_DIV_30\';
        RF_path='E:\MATLAB\project\DE\save_data\RANK_FIT_30\';
    elseif D==50
        RD_path='E:\MATLAB\project\DE\save_data\RANK_DIV_50\';
        RF_path='E:\MATLAB\project\DE\save_data\RANK_FIT_50\';
    end
    RD_filename=strcat('RD','_',int2str(D),'D_F',int2str(k),'.mat');
    load([RD_path,RD_filename]); 
    RF_filename=strcat('RF','_',int2str(D),'D_F',int2str(k),'.mat');
    load([RF_path,RF_filename]);
    
for i=1:runNumber
    fprintf("------ASODE���ҶԱȵ�%d������------\n",i);
    fprintf("ABM-1--->");
    % ABM-1  DE/rand/1
    [Pb,~,FEs_fitness]=ASODE_1(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,RANK_DIV,RANK_FIT,func_num);
    ASODE_1Matrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
    
    fprintf("ABM-2--->");
    % ABM-2   ASODE/ABM-2
    [Pb,~,FEs_fitness]=ASODE_2(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,RANK_DIV,RANK_FIT,func_num);
    ASODE_2Matrix(i,:)=FEs_fitness;
    s(2,i)=Pb;
    
    fprintf("ICA-1--->");
    % ICA 3:5:2
    [Pb,~,FEs_fitness]=ASODE_3(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,RANK_DIV,RANK_FIT,func_num);
    ASODE_3Matrix(i,:)=FEs_fitness;
    s(3,i)=Pb;
    
    fprintf("ICA-2\n");
    % ICA 1:5:4
    [Pb,~,FEs_fitness]=ASODE_4(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,RANK_DIV,RANK_FIT,func_num);
    ASODE_4Matrix(i,:)=FEs_fitness;
    s(4,i)=Pb;
    
end
msASODE_FES(1,:)= mean(ASODE_1Matrix,1);
msASODE_FES(2,:)= mean(ASODE_2Matrix,1);
msASODE_FES(3,:)= mean(ASODE_3Matrix,1);
msASODE_FES(4,:)= mean(ASODE_4Matrix,1);

if D==30
    path='E:\MATLAB\project\DE\save_data\myself_comparsion_data\msASODE_30D\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\myself_comparsion_data\msASODE_50D\';
end
filename=strcat('msASODE','_',int2str(D),'D_F',int2str(k),'.mat');
save([path,filename],'msASODE_FES');
% 
fprintf("\nASODE_1:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(1,:)),std(s(1,:)));
fprintf("\nASODE_2:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(2,:)),std(s(2,:)));
fprintf("\nASODE_3:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(3,:)),std(s(3,:)));
fprintf("\nASODE_4:\n��ֵΪ:%d\nstdΪ:%d\n",mean(s(4,:)),std(s(4,:)));
% %%% ����ASODE���ݵ�excel�ļ� %%%
myself_com_savedata(k,D,s); 


% ����ASODE���ݣ���������ͼ
if D==10
    path='E:\MATLAB\project\DE\save_data\ASODE_data_10\ASODE_10\';
elseif D==30
    path='E:\MATLAB\project\DE\save_data\ASODE_data_30\ASODE_30\';
elseif D==50
    path='E:\MATLAB\project\DE\save_data\ASODE_data_50\ASODE_50\';
end
filename=strcat('ASODE','_',int2str(D),'D_F',int2str(k),'.mat');
load([path,filename]);

figure(k);
xx = 1:Max_FES;
plot(xx,log(ASODE_FES(1,xx)), 'c-p',...
     xx,log(msASODE_FES(1,xx)),'k-o',... 
     xx,log(msASODE_FES(2,xx)), 'b-d',...
     xx,log(msASODE_FES(3,xx)), 'g-^',...
     xx,log(msASODE_FES(4,xx)), 'm-s',...
     'MarkerIndices',1:(Max_FES/10):Max_FES,'LineWidth',1);
legend('ASODE','ASODE\_1','ASODE\_2','ASODE\_3','ASODE\_4');
xlabel('Number of Function Evaluations');
ylabel('Average Function Values(Log)');

end