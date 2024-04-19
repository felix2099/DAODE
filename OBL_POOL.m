% % -------------------------------------- % %
% % OBL POOL                             % %
% % -------------------------------------- % %
function zh_ox = OBL_POOL(D,NP,pop_rank,count,x,r,i,Xb,T)

a=min(x);
b=max(x);
F_SORT = pop_rank(i);
K = F_SORT/NP;

if T == 3
    rank_random = randi([1 2]);
    % % ---------------ODE----------------------- % %
    if (r(count,2) <= 2 && r(count,2) >= 1 && r(count,2) == rank_random) 
        ox(i,:) = a + b - x(i,:);
        zh_ox = ox(i,:);
    % % ---------------QODE----------------------- % %
    elseif (r(count,3) <= 2 && r(count,3) >= 1 && r(count,3) == rank_random)  
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
    elseif (r(count,4) <= 2 && r(count,4) >= 1 && r(count,4) == rank_random)  
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
    elseif (r(count,5) <= 2 && r(count,5) >= 1 && r(count,5) == rank_random) 
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
    elseif (r(count,6) <= 2 && r(count,6) >= 1 && r(count,6) == rank_random)  
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
    elseif (r(count,7) <= 2 && r(count,7) >= 1 && r(count,7) == rank_random)   
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
    elseif (r(count,8) <= 2 && r(count,8) >= 1 && r(count,8) == rank_random)  
          for j = 1 : D
             ox(i,j) = rand() * (a(j) + b(j)) - x(i,j);
             if ox(i,j) < a(j) 
                 ox(i,j)=(b(j)-a(j))*rand + a(j);
             end
             if ox(i,j) > b(j)
                 ox(i,j) = (b(j) - a(j)) * rand + a(j);
             end
           end
           zh_ox = ox(i,:);

    % % ----------------RODE---------------------- % %
    elseif (r(count,9) <= 2 && r(count,9) >= 1 && r(count,9) == rank_random) 
        beta=180 * normrnd(1,0.25);
        for j = 1 : D
            o(j) = (a(j) + b(j)) / 2;
            u(i,j) = (x(i,j) - o(j));
            v(i,j) = sqrt((x(i,j) - a(j)) * (b(j) - x(i,j)));
            ox(i,j) = o(j) + (u(i,j) * cosd(beta)- v(i,j) * sind(beta));
        end    
        zh_ox = ox(i,:);

    % % ----------------CODE--------------------- % %
    elseif (r(count,10) <= 2 && r(count,10) >= 1 && r(count,10) == rank_random)  
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
    elseif (r(count,11) <= 2 && r(count,11) >= 1 && r(count,11) == rank_random)  
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
    elseif (r(count,12) <= 2 && r(count,12) >= 1 && r(count,12) == rank_random)  
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
    elseif (r(count,13) <= 2 && r(count,13) >= 1 && r(count,13) == rank_random)  
        for j = 1 : D
           o(j) = (a(j) + b(j))/2;
           if x(i,j) <= o(j) & x(i,j) >= a(j)
                ox(i,j) = (x(i,j) + 2*b(j))/3;
           else
                ox(i,j) = (x(i,j) + 2*a(j))/3;
           end
        end    
        zh_ox = ox(i,:);

    % % ----------------COODE---------------------- % %
    elseif (r(count,14) <= 2 && r(count,14) >= 1 && r(count,14) == rank_random)  
        ox(i,:) = 2 * Xb- x(i,:);         
        if ((ox(i,:) < a) | (ox(i,:) > b))
            ox(i,:) = a + rand * (b-a);
        end
        zh_ox = ox(i,:);
    end
    
elseif T == 2
    rank_random = round(normrnd(13/2,7/6,1,1));
%     rank_random = randi([3 10]);
    % % ---------------ODE----------------------- % %
    if (r(count,2) <= 10 && r(count,2) >= 3 && r(count,2) == rank_random) 
        ox(i,:) = a + b - x(i,:);
        zh_ox = ox(i,:);
    % % ---------------QODE----------------------- % %
    elseif (r(count,3) <= 10 && r(count,3) >= 3 && r(count,3) == rank_random)
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
    elseif (r(count,4) <= 10 && r(count,4) >= 3 && r(count,4) == rank_random) 
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
    elseif (r(count,5) <= 10 && r(count,5) >= 3 && r(count,5) == rank_random) 
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
    elseif (r(count,6) <= 10 && r(count,6) >= 3 && r(count,6) == rank_random) 
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
    elseif (r(count,7) <= 10 && r(count,7) >= 3 && r(count,7) == rank_random) 
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
    elseif (r(count,8) <= 10 && r(count,8) >= 3 && r(count,8) == rank_random) 
          for j = 1 : D
             ox(i,j) = rand() * (a(j) + b(j)) - x(i,j);
             if ox(i,j) < a(j) 
                 ox(i,j)=(b(j)-a(j))*rand + a(j);
             end
             if ox(i,j) > b(j)
                 ox(i,j) = (b(j) - a(j)) * rand + a(j);
             end
           end
           zh_ox = ox(i,:);

    % % ----------------RODE---------------------- % %
    elseif (r(count,9) <= 10 && r(count,9) >= 3 && r(count,9) == rank_random)
        beta=180 * normrnd(1,0.25);
        for j = 1 : D
            o(j) = (a(j) + b(j)) / 2;
            u(i,j) = (x(i,j) - o(j));
            v(i,j) = sqrt((x(i,j) - a(j)) * (b(j) - x(i,j)));
            ox(i,j) = o(j) + (u(i,j) * cosd(beta)- v(i,j) * sind(beta));
        end    
        zh_ox = ox(i,:);

    % % ----------------CODE--------------------- % %
    elseif (r(count,10) <= 10 && r(count,10) >= 3 && r(count,10) == rank_random) 
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
    elseif (r(count,11) <= 10 && r(count,11) >= 3 && r(count,11) == rank_random) 
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
    elseif (r(count,12) <= 10 && r(count,12) >= 3 && r(count,12) == rank_random) 
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
    elseif (r(count,13) <= 10 && r(count,13) >= 3 && r(count,13) == rank_random) 
        for j = 1 : D
           o(j) = (a(j) + b(j))/2;
           if x(i,j) <= o(j) & x(i,j) >= a(j)
                ox(i,j) = (x(i,j) + 2*b(j))/3;
           else
                ox(i,j) = (x(i,j) + 2*a(j))/3;
           end
        end    
        zh_ox = ox(i,:);

    % % ----------------COODE---------------------- % %
    elseif (r(count,14) <= 10 && r(count,14) >= 3 && r(count,14) == rank_random) 
        ox(i,:) = 2 * Xb- x(i,:);         
        if ((ox(i,:) < a) | (ox(i,:) > b))
            ox(i,:) = a + rand * (b-a);
        end
        zh_ox = ox(i,:);

    end
    
else
        rank_random = randi([11 13]);
    % % ---------------ODE----------------------- % %
    if (r(count,2) <= 13 && r(count,2) >= 11 && r(count,2) == rank_random)
        ox(i,:) = a + b - x(i,:);
        zh_ox = ox(i,:);
    % % ---------------QODE----------------------- % %
    elseif (r(count,3) <= 13 && r(count,3) >= 11 && r(count,3) == rank_random)
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
    elseif (r(count,4) <= 13 && r(count,4) >= 11 && r(count,4) == rank_random)
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
    elseif (r(count,5) <= 13 && r(count,5) >= 11 && r(count,5) == rank_random) 
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
    elseif (r(count,6) <= 13 && r(count,6) >= 11 && r(count,6) == rank_random)
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
    elseif (r(count,7) <= 13 && r(count,7) >= 11 && r(count,7) == rank_random) 
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
    elseif (r(count,8) <= 13 && r(count,8) >= 11 && r(count,8) == rank_random)
          for j = 1 : D
             ox(i,j) = rand() * (a(j) + b(j)) - x(i,j);
             if ox(i,j) < a(j) 
                 ox(i,j)=(b(j)-a(j))*rand + a(j);
             end
             if ox(i,j) > b(j)
                 ox(i,j) = (b(j) - a(j)) * rand + a(j);
             end
           end
           zh_ox = ox(i,:);

    % % ----------------RODE---------------------- % %
    elseif (r(count,9) <= 13 && r(count,9) >= 11 && r(count,9) == rank_random)
        beta=180 * normrnd(1,0.25);
        for j = 1 : D
            o(j) = (a(j) + b(j)) / 2;
            u(i,j) = (x(i,j) - o(j));
            v(i,j) = sqrt((x(i,j) - a(j)) * (b(j) - x(i,j)));
            ox(i,j) = o(j) + (u(i,j) * cosd(beta)- v(i,j) * sind(beta));
        end    
        zh_ox = ox(i,:);

    % % ----------------CODE--------------------- % %
    elseif (r(count,10) <= 13 && r(count,10) >= 11 && r(count,10) == rank_random) 
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
    elseif (r(count,11) <= 13 && r(count,11) >= 11 && r(count,11) == rank_random)
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
    elseif (r(count,12) <= 13 && r(count,12) >= 11 && r(count,12) == rank_random) 
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
    elseif (r(count,13) <= 13 && r(count,13) >= 11 && r(count,13) == rank_random) 
        for j = 1 : D
           o(j) = (a(j) + b(j))/2;
           if x(i,j) <= o(j) & x(i,j) >= a(j)
                ox(i,j) = (x(i,j) + 2*b(j))/3;
           else
                ox(i,j) = (x(i,j) + 2*a(j))/3;
           end
        end    
        zh_ox = ox(i,:);

    % % ----------------COODE---------------------- % %
    elseif (r(count,14) <= 13 && r(count,14) >= 11 && r(count,14) == rank_random) 
        ox(i,:) = 2 * Xb- x(i,:);         
        if ((ox(i,:) < a) | (ox(i,:) > b))
            ox(i,:) = a + rand * (b-a);
        end
        zh_ox = ox(i,:);

    end
    
end  

end

