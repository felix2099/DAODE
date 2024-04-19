
function [Pb,trace,FEs_fitness]=ASODE(func_num,fhd,D,NP,F,CR,gen_max,Max_FES,border,RANK_DIV,RANK_FIT,varargin)

fbias=[100,200,300,400,500,600,700,...
       800,900,1000,1100,1200,1300,...
       1400,1500,1600,1700,1800,1900,...
       2000,2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000];
eps = 1e-10;
trace=zeros(gen_max,2);   
bounds=border*ones(D,2);  %设置边界
bounds(:,1)=-1*bounds(:,1);
rng=(bounds(:,2)-bounds(:,1))';    
x=(ones(NP,1)*rng).*(rand(NP,D))+(ones(NP,1)*bounds(:,1)');%初始种群
r=zeros(gen_max,14);
Pb=inf; %cost(1);%存放最优值
Xb=x(1,:);%存放最优位置
count=1;
cost=zeros(1,NP);

for i=1:NP         
    cost(1,i)=feval(fhd,x(i,:)',varargin{:})-fbias(func_num);
    if(cost(1,i) <= Pb)
        Pb=cost(1,i);
        Xb=x(i,:);                
    end
end

%种群初始化采用多样性最大的策略
for i = 1:NP
    ox(i,:) = INI_OBL_POOL(D,NP,cost,RANK_DIV,count,x,i,Xb,Pb);
    if (feval(fhd,ox(i,:)',varargin{:})-fbias(func_num)) < (feval(fhd,x(i,:)',varargin{:})-fbias(func_num))
            x(i,:) = ox(i,:);
    end
end

trial=zeros(1,D);
for i=1:NP
    cost(1,i)=feval(fhd,x(i,:)',varargin{:})-fbias(func_num);
    if(cost(1,i)<=Pb)
        Pb=cost(1,i);
        Xb=x(i,:);                
   end
end

fitFEs_count = NP;
initial_FEs = 1;
new_FEs = fitFEs_count;
FEs_fitness(initial_FEs:new_FEs) = Pb;
old_FEs = new_FEs;

trace(1,1)=1;
trace(1,2)=Pb;

for count = 2 : gen_max
    
    if fitFEs_count > Max_FES
        break;
    end
    
%     superior_archive_gbest = zeros(0.2*NP,D); %此处设置存档，保存每代的优势个体 
%     middle_archive_gbest = zeros(0.5*NP,D);   %此处设置存档，保存每代的中间个体
%     inferior_archive_gbest = zeros(0.3*NP,D); %此处设置存档，保存每代的劣势个体
    
    %为排序结果分配内存，加快运行速度，每代运行后归零
    superior_fit = zeros(1,NP);
    middle_fit = zeros(1,NP);
    inferior_fit = zeros(1,NP);
    %在此对种群个体进行排序
    [pop_value,pop_index] = sort(cost); % pop_value 排序后适应值；pop_index 排序后适应值的索引值。
    [~,pop_rank] = sort(pop_index);     % pop_rank 得出种群适应值排名  排名对应于cost原始值

    for i = 1 : length(cost)
        if pop_rank(i) <= 20  % 优势个体排名为[1,20]  
           superior_fit(1,i) = pop_value(pop_rank(i));   % 此处得到的是适应值，而非个体的位置 
           superior_fit(superior_fit == 0) = [];       % 去除数组中0元素
           j = 0;
           for ii = 1 : NP
                if (cost(1,ii) <= max(superior_fit))  % 找到对应的个体
                    j = j + 1;
                    superior_archive_gbest(j,:) = x(ii,:); 
                end
           end           
           
        elseif (pop_rank(i) >= 21)  && (pop_rank(i) <= 70)  % 中间个体排名为[21,70]            
           middle_fit(1,i) = pop_value(pop_rank(i));    
           middle_fit(middle_fit == 0) = [];
           k = 0;
           for ii = 1 : NP
                if (cost(1,ii) >= min(middle_fit) & cost(1,ii) <= max(middle_fit))  % 找到对应的个体
                    k = k + 1;
                    middle_archive_gbest(k,:) = x(ii,:);                    
                end
           end
           
        else                                            % 劣势个体排名为[71,100]            
           inferior_fit(1,i) = pop_value(pop_rank(i));
           inferior_fit(inferior_fit == 0) = [];
           m = 0;
           for ii = 1 : NP
                if (cost(1,ii) >= min(inferior_fit))  % 找到对应的个体
                    m = m + 1;
                    inferior_archive_gbest(m,:) = x(ii,:);            
                end
           end      
        end           
    end
    
    for i=1:NP
     
        c=randi([1 size(middle_archive_gbest,1)]);
        a=randi([1 size(superior_archive_gbest,1)]);
        b=randi([1 size(inferior_archive_gbest,1)]);

        jrand=floor(rand*D+1);
        for k=1:D
            if(rand<CR||jrand==k)
%                 trial(k)=x(c,k)+F*(x(a,k)-x(b,k));
%                 trial(k)=middle_archive_gbest(c,k)+F*(superior_archive_gbest(a,k)-inferior_archive_gbest(b,k));
                trial(k)=superior_archive_gbest(a,k)+F*(middle_archive_gbest(c,k)-inferior_archive_gbest(b,k));
            else
                trial(k)=x(i,k);
            end
            if trial(k)<bounds(k,1)
                    trial(k)=bounds(k,1);
            end
            if trial(k)>bounds(k,2)
                   trial(k)=bounds(k,2);
            end
        end

        %读取多样性和适应度排名
%         p=count/gen_max; %sqrt(2*count/gen_max); 
        p = 0.5+0.5*cos(i/gen_max*(D/10)*pi); %三角函数
%         p = acos(2*i/gen_max-1)/pi;   %反三角函数
        [~,index] = sort((p * RANK_DIV(count,2:14) + (1-p) * RANK_FIT(count,2:14)));   %综合排名
        [~,r(count,2:14)]=sort(index);
       
        if ((feval(fhd,trial(:),varargin{:})-fbias(func_num)) <= max(superior_fit))   % 判断此个体是否属于优势个体
            T = 1;
            OBL_POOL(D,NP,pop_rank,count,x,r,i,Xb,T);   % 选择多样性较小和收敛较快的策略，排名略差的策略
             
        elseif ((feval(fhd,trial(:),varargin{:})-fbias(func_num)) <= max(middle_fit)) & ((feval(fhd,trial(:),varargin{:})-fbias(func_num)) >= min(middle_fit))  % 判断此个体是否属于中间个体
            T = 2;
            OBL_POOL(D,NP,pop_rank,count,x,r,i,Xb,T);   % 选择多样性和收敛速度中间的策略，高斯产生随机策略            
             
        elseif ((feval(fhd,trial(:),varargin{:})-fbias(func_num)) >= min(inferior_fit))   % 判断此个体是否属于劣势个体
            T = 3;
            OBL_POOL(D,NP,pop_rank,count,x,r,i,Xb,T);   % 选择排名较好的策略
        end              

        
        oxscore=feval(fhd,ox(i,:)',varargin{:})-fbias(func_num);                            
        trialscore=feval(fhd,trial(:),varargin{:})-fbias(func_num);  
        fitFEs_count = fitFEs_count + 1;    %记录评估次数
        xiscore=cost(i);
        
        xmin=xiscore;
        if xmin>trialscore
            x(i,1:D)=trial(1:D);
            cost(i)=trialscore;
            xmin=trialscore;
        end
        if xmin>oxscore
            x(i,1:D)=ox(1:D);
            cost(i)=oxscore;
        end
              
        if cost(i) <= Pb
            Pb=cost(i);
%             if Pb <= eps
%                 Pb = 0;
%                 cost(i) = 0;
%             end
            Xb(1:D)=x(i,1:D);
        end  
    new_FEs = fitFEs_count;
    FEs_fitness(old_FEs:new_FEs) = Pb;
    old_FEs = new_FEs;        
    end  
    
    trace(count,1)=count;
    trace(count,2)=Pb;
    
end  %迭代循环结束

end

