% % -------------------------------------- % %
% % Initial Population OBL Strategy Pool   % %
% % -------------------------------------- % %
function zh_ox = INI_OBL_POOL(D,NP,cost,RANK_DIV,count,x,i,Xb,Pb)

a=min(x);
b=max(x);

% % FQRODE % %
[~,index] = sort(cost);
[~,F_index] = sort(index);
for ii=1:NP
    if (cost(ii) == Pb)
        F_SORT = F_index(ii);
        break;
    end
end
K = F_SORT/NP;

% % ---------------ODE----------------------- % %
if RANK_DIV(count,2) == 1 
    ox(i,:) = a + b - x(i,:);
    zh_ox = ox(i,:);
    
% % ---------------QODE----------------------- % %
elseif RANK_DIV(count,3) == 1 
     for j = 1 : D
         oox(i,j) = a(j) + b(j) - x(i,j);
         M(j) = (a(j) + b(j)) / 2;
         if (x(i,j) < M(j))
             ox(i,j) = M(j) + rand * (oox(i,j) - M(j));
         else
             ox(i,j) = M(j) + rand * (M(j) - oox(i,j));
         end      
     end 
     zh_ox = ox(i,:);
     
% % --------------QRODE----------------------- % %
elseif RANK_DIV(count,4) == 1 
      for j = 1 : D
          M(j) = (a(j) + b(j)) / 2;
          if (x(i,j) < M(j))
              ox(i,j) = M(j) + rand * (x(i,j)-M(j));
           else
              ox(i,j) = M(j) + rand * (M(j)-x(i,j));
           end
      end  
      zh_ox = ox(i,:);
      
% % ---------------EODE---------------------- % %
elseif RANK_DIV(count,5) == 1 
      for j = 1 : D
          o(j) = (a(j) + b(j)) / 2;
          oox(i,j) = a(j) + b(j) - x(i,j);
          if x(i,j) < o(j)
              ox(i,j) = oox(i,j) + rand * (b(j) - oox(i,j));
          elseif x(i,j) > o(j)
              ox(i,j) = a(j) + rand * (oox(i,j) - a(j));
          else
              ox(i,j) = a(j) + b(j) - x(i,j);
          end
      end    
      zh_ox = ox(i,:);
      
% % ---------------REODE----------------------- % %
elseif RANK_DIV(count,6) == 1 
    for j = 1 : D
        o(j) = (a(j) + b(j))/2;
        if x(i,j) > o(j)
            ox(i,j) = a(j) + rand * (x(i,j)-a(j));
        elseif x(i,j) < o(j)
            ox(i,j) = x(i,j) + rand * (b(j)-x(i,j));
        else
            ox(i,j) = a(j) + b(j) - x(i,j);
        end
    end    
    zh_ox = ox(i,:);
    
% % ---------------FQRODE---------------------- % %
elseif RANK_DIV(count,7) == 1           
      for j = 1 : D
          o(j) = (a(j) + b(j))/2;
          oox(i,j) = a(j) + b(j) - x(i,j);
          if oox(i,j) <= o(j)
              ox(i,j) = oox(i,j) + (o(j)-oox(i,j)) * K;
          else
              ox(i,j) = o(j) + (oox(i,j) - o(j)) * (1 - K);
          end
      end     
      zh_ox = ox(i,:);
      
% % ---------------GODE---------------------- % %
elseif RANK_DIV(count,8) == 1 
      for j = 1 : D
         ox(i,j) = rand() * (a(j) + b(j)) - x(i,j);
         if ox(i,j) < a(j) 
             ox(i,j)=(b(j)-a(j))*rand() + a(j);
         end
         if ox(i,j) > b(j)
             ox(i,j) = (b(j) - a(j)) * rand() + a(j);
         end
      end
      zh_ox = ox(i,:);

% % ----------------RODE---------------------- % %
elseif RANK_DIV(count,9) == 1
    beta=180 * normrnd(1,0.25);
    for j = 1 : D
        o(j) = (a(j) + b(j)) / 2;
        u(i,j) = (x(i,j) - o(j));
        v(i,j) = sqrt((x(i,j) - a(j)) * (b(j) - x(i,j)));
        ox(i,j) = o(j) + (u(i,j) * cosd(beta)- v(i,j) * sind(beta));
    end    
    zh_ox = ox(i,:);
    
% % ----------------CODE--------------------- % %
elseif RANK_DIV(count,10) == 1 
      for k = 1 : D  
          M(k) = 0;
          for j = 1 : NP
              M(k) = M(k) + x(j,k); 
          end 
              M(k)=M(k) / NP;
      end
            
      for k = 1 : D
          ox(i,k) = 2 * M(k) - x(i,k);
          if (ox(i,k) > b(k))
              d=M(k);
              c=b(k);                          
              ox(i,k)=d + rand * (c-d);    
          end
          if (ox(i,k) < a(k))           
              d=a(k); 
              c=M(k);
              ox(i,k)=d + rand * (c-d);    
          end             
      end
      zh_ox = ox(i,:);
      
% % ----------------PODE---------------------- % %
elseif RANK_DIV(count,11) == 1 
    d = round(1 + (D - 1) * rand(1,1));
    for j = 1 : D
        if j == d
            ox(i,j) = x(i,j);
        else
            ox(i,j) = a(j) + b(j) - x(i,j);
        end
    end    
    zh_ox = ox(i,:);
    
% % ----------------SODE--------------------- % %
elseif RANK_DIV(count,12) == 1 
    for j = 1 : D
        o(j)=(a(j) + b(j))/2;
        oox(i,j) = a(j) + b(j) - x(i,j);
        if x(i,j) > o(j)
            ox(i,j) = a(j) + (oox(i,j) - a(j)) * rand(1,1);
        elseif x(i,j) == o(j)
            ox(i,j) = (a(j) + (b(j) - a(j)) * rand(1,1)) - x(i,j);
        else                   
            ox(i,j) = oox(i,j) + (b(j) - oox(i,j)) * rand(1,1);
        end
    end    
    zh_ox = ox(i,:);
    
% % ----------------OCDE--------------------- % %
elseif RANK_DIV(count,13) == 1 
    for j = 1 : D
       o(j) = (a(j) + b(j))/2;
       if x(i,j) <= o(j) & x(i,j) >= a(j)
            ox(i,j) = (x(i,j) + 2*b(j))/3;
       else
            ox(i,j) = (x(i,j) + 2*a(j))/3;
       end
    end    
    
% % ----------------COODE---------------------- % %
elseif RANK_DIV(count,14) == 1 
    ox(i,:) = 2 * Xb- x(i,:);         
    if ((ox(i,:) < a) | (ox(i,:) > b))
        ox(i,:) = a + rand * (b-a);
    end
    zh_ox = ox(i,:);
    
end

                           
end

