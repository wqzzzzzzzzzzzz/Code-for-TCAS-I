%% Methods about the allocation of weights to the nodes in the network topology 
%% In this paper, we We employ the Uniform criterion (type=1) in this paper, where the weight of each node is the same
function [C] = Weight(num ,A,SIGMA_EM,WEIGHT_EM,type,Sig,orders)
%%
if type==1 % Uniform  
 for i=1:num
    temp=find(A(i,:));     
    for j=1:length(temp)
        C(i,temp(j))= 1/length(temp);
    end        
 end
end
%%
if type==2 %N
  for i=1:num
      Ad=0;
      Max_N=0;
      temp=find(A(i,:));
      if length(temp)> Max_N
        Max_N=length(temp);
      end      
  end
  
  for i=1:num
   temp=find(A(i,:));   
   for j=1:length(temp)
    C(i,temp(j))= 1/Max_N;
    if i==temp(j)
      C(i,temp(j))= 1-(length(temp)-1)/Max_N;   
    end   
   end
  end
end

%%
if type==3 %n_max
  for i=1:num
   Ad=0;
   temp=find(A(i,:));
   Max_N=length(temp);
    for j=1:length(temp)
       if length(find(A(temp(j),:))) >Max_N
           Max_N=length(find(A(temp(j),:)));
       end 
    end
    
    for j=1:length(temp)
        C(i,temp(j))=1/Max_N;
        if i==temp(j)
           C(i,temp(j))= 1-(length(temp)-1)/Max_N; 
        end
    end   
  end 
end
%%
if type==4 %metropolis           
 for i=1:num
        temp=find(A(i,:));
        Ad=0;
        for j=1:length(temp)
            C(i,temp(j))=1/max(length(find(A(temp(j),:))),length(temp));
            if i~=temp(j)
            Ad=Ad+C(i,temp(j));
            end
        end
        C(i,i)=1-Ad;
 end
end

%%
if type==5 %                
 for i=1:num
         Ad=0;
        temp=find(A(i,:));
        for j=1:length(temp)
            Ad=Ad+Sig(temp(j))^(-1);
        end     
        for j=1:length(temp)
            C(i,temp(j))= Sig(temp(j))^(-1)/Ad;
        end
        
 end
end
%%
if type==6 %Relative degree
 for i=1:num
        temp=find(A(i,:));
        Ad=0;
        for j=1:length(temp)
           Ad=Ad+length(find(A(temp(j),:)));
        end
        for j=1:length(temp)
             C(i,temp(j))= length(find(A(temp(j),:)))/Ad; 
        end
 end
end
%%
if type==7 %R-d-v                
 for i=1:num
         Ad=0;
        temp=find(A(i,:));
        for j=1:length(temp)
            Ad=Ad + length(find(A(temp(j),:)))*Sig(temp(j))^(-1);
        end     
        for j=1:length(temp)
            C(i,temp(j))= length(find(A(temp(j),:)))*Sig(temp(j))^(-1)/Ad;
        end        
 end
end
%%
%%
if type==8 %            
 for i=1:num
         Ad=0;
        temp=find(A(i,:));
        for j=1:length(temp)
            Ad=Ad + length(find(A(temp(j),:)))*Sig(temp(j));
        end     
        for j=1:length(temp)
            C(i,temp(j))= length(find(A(temp(j),:)))*Sig(temp(j))/Ad;
        end        
 end
end
end


