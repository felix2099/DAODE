clear;
clc; 
close all; 
addpath(genpath(pwd));  
runNumber=30; 
D=50;         
NP=100;       
F=0.5;        
CR=0.9;       
Max_FES = 10000 * D;
gen_max = Max_FES / NP;   
border=100;     
%func_num=1;
fhd=str2func('cec17_func');
str = "DAODE"; 

DAODEMatrix=zeros(runNumber,Max_FES);

s=zeros(1,runNumber);
DAODE_FES = zeros(1,Max_FES);

for k=1:30
    func_num=k;
    if k==2
       continue;
    end
    fprintf("--------------------------\n");
    fprintf("The Initiation of Testing %s's %d-Dimensional -F%d Function >>>>\n",str,D,k);
    fprintf("--------------------------\n");
    
    %%% Read diversity and fitness ranking data%%%
    if D==10 
        RD_path='..\DAODE\save_data\RANK_DIV_10\';
        RF_path='..\DAODE\save_data\RANK_FIT_10\';
    elseif D==30
        RD_path='..\DAODE\save_data\RANK_DIV_30\';
        RF_path='..\DAODE\save_data\RANK_FIT_30\';
    elseif D==50
        RD_path='..\DAODE\save_data\RANK_DIV_50\';
        RF_path='..\DAODE\save_data\RANK_FIT_50\';
    end
    RD_filename=strcat('RD','_',int2str(D),'D_F',int2str(k),'.mat');
    load([RD_path,RD_filename]); 
    RF_filename=strcat('RF','_',int2str(D),'D_F',int2str(k),'.mat');
    load([RF_path,RF_filename]); 
    
for i=1:runNumber
    fprintf("------DAODE:The %dth run------\n",i);
    [Pb,~,FEs_fitness]=DAODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,RANK_DIV,RANK_FIT,func_num);
    DAODEMatrix(i,:)=FEs_fitness;
    s(1,i)=Pb;
end
DAODE_FES(1,:)= mean(DAODEMatrix,1);

fprintf("\nDAODE:\nmean is:%d\nstd is:%d\n",mean(s(1,:)),std(s(1,:)));

end
