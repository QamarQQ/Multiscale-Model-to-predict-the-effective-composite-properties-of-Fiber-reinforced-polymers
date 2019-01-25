%INPUTS

%number of layers

N=5;

%width 
b= 10;

%Thickness
t=[.5  .5  .5 .5 .5];

%Applied shear load
Q=1;

%layerwise youngs modulus
E=[10 20 30 40 50];

%Calculation of boundaries of layers : z_1(i) or zk

c=sum(t)/2;
z_1(1)=c;
for i=1:N
    
    z_1(i+1)=z_1(i)-(t(i))
end

clear i;

%Calculation of normalization factors : n(i) 

E_a=sum(E)/N;
for i=1:N
    n(i)=E(i)/E_a;
end

clear i;

%Calculation of nuetral axes : z_g or z*
for i=1:N
    c(i)=((z_1(i+1)+z_1(i))*n(i)*t(i));
end

clear i;

a=sum(c)/2;
for i=1:N
    C1(i)=(t(i)*n(i));
end
w=sum(C1);
z_g = a/w;


%Calculation of Boundaries wrt Nuetral axes or z or z_cap

 for i=1:N+1;
 z(i) = z_1(i)-z_g;
end

clear i;

%Calculation of Moment of Inertia

for i=1:N
e(i)= (b/12)* n(i)*(t(i)*t(i)*t(i));
f(i)= (b/4)*n(i)*(t(i))*((z(i+1)+z(i))^2);
end

I = sum(e)+sum(f);

clear i;


% Layerwise/lamina wise shear stress distribution Calculation.

for i=1:N
    
    
    t_21=0;
    t_22=0;
    t_23=0;
    t_2_1=0;
    t_2_2=0;
    t_2_3=0;
    
    r=((z_1(i)-(z_1(i+1)))/100);
    x(:,:,i) = [ z_1(i) : -r :(z_1(i+1))]
 
    
    if z_1(i)*z_1(i+1)>=0 
      
      % t-1 calutaion starts  
        if z_1(i)<0
            h1=-1*z_1(i);
        else
        h1=z_1(i);
        end
    
        if z_1(i+1)<0
           h2=-1*z_1(i+1);
        else
           h2=z_1(i+1);
        end 
    
        if h1>h2
          H= z(i)
        else 
         H= z(i+1) 
        end
        t_1(:,:,i) = (b/2)*n(i)*((x(:,:, i )-z_g).^2-(H)^2)
    
     % t-1 calutaion ends
     
     % t-2 calutaion starts 
     
                                     if z_1(i)*z_1(i+1)>0 && z_1(i)>0
                                     display case1
                                     k=i-1;
                                     m=k+1;
                                     
                                     
                                         
                                     for m=k+1:N
                                     t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     
                                     
                                     t_21=t_2_1;
                                     t_2_1=0;
                                     t_2(i)=t_21
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                  
                                  elseif z_1(i)*z_1(i+1)>0 && z_1(i)<0
                                      display case2
                                      
             
                                     k=i;
                                     m=k+1;
                                     
                                     
                                     for m=k+1:N
                                     t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     t_22=t_2_2;
                                     t_2_2=0;
                                     t_2(i)=t_22
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                   
                                  elseif  z_1(i)*z_1(i+1)==0 && z_1(i)==0
                                      display case3
                                      k=i;
                                     m=k+1;
                                     
                                     
                                     for m=k+1:N
                                     t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     t_22=t_2_2;
                                     t_2_2=0;
                                     t_2(i)=t_22
                                      Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                      T_x(:,:,i) = -(Q_x(:,:,i)*Q)/I
                                   
                                  elseif z_1(i)*z_1(i+1)==0 && z_1(i+1)==0
                                     display case4
                                     k=i-1;
                                     m=k+1;
                                     
                                     
                                         
                                     for m=k+1:N
                                     t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                     end
                                     
                                     
                                     t_21=t_2_1;
                                     t_2_1=0;
                                     t_2(i)=t_21
                                     Q_x(:,:,i) = t_1(:,:,i) + t_2(i); 
    
                                     T_x(:,:,i) = -(Q_x(:,:,i)*Q)/I
     
                                     end
     
                                     
         % t-2 calutaion ends                            
     
    else
        
        % t-1 and t-2 calutaion starts
        for v= 1:length(x)
           
         if x(:,v,i)>=0
         display case6  
         k=i-1;
         m=k+1;
         t_1(:,v,i) = (b/2)*n(i)*((x(:,v, i )-z_g).^2-(z(m))^2)
                                        
                                     
                                     
                                         
                                        for m=k+1:N
                                        t_2_1=t_2_1+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                        end
                                     
                                     
                                        t_21=t_2_1
                                        t_2_1=0;
                                        t_2(i)=t_21
                                        Q_x(:,v,i) = t_1(:,v,i) + t_2(i); 
    
                                        T_x(:,v,i) = -(Q_x(:,v,i)*Q)/I
         
         else
         display case7
         k=i;
         m=k+1;
         t_1(:,v,i) = (b/2)*n(i)*((x(:,v, i )-z_g).^2-(z(m))^2)
                                        k=i;
                                        m=k+1;
                                     
                                     
                                        for m=k+1:N
                                        t_2_2=t_2_2+(n(m)*(b/2)*(z(m)^2-z(m+1)^2));
         
                                        end
                                        t_22=t_2_2
                                        t_2_2=0;
                                        t_2(i)=t_22
                                        Q_x(:,v,i) = t_1(:,v,i) + t_2(i); 
    
                                        T_x(:,v,i) = -(Q_x(:,v,i)*Q)/I
       
         end
         
        v=v+1; 
        
        end
        
        
       % t-1 and t-2  calutaion ends
    end  
     
                                     
                                   
end



M=max(T_x(:,:,:))

% PLOTTING %   

% figure (1)
% plot (  T_x(:,:,1), x(:,:, 1) ) ;
 
% figure (1)
% plot (  T_x(:,:,1), x(:,:, 1) , T_x(:,:,2), x(:,:, 2) );
 
 
 %figure (1)
 %plot (  T_x(:,:,1), x(:,:, 1) , T_x(:,:,2), x(:,:, 2) ,T_x(:,:,3) , x(:,:, 3) );
  

% figure (2)
% plot (  T_x(:,:,1), x(:,:, 1) , T_x(:,:,2), x(:,:, 2) ,T_x(:,:,3) ,x(:,:, 3) , T_x(:,:,4) ,x(:,:, 4) ) ;

%figure (2)
 %plot (  T_x(:,:,1), x(:,:, 1) , T_x(:,:,2), x(:,:, 2) ,T_x(:,:,3) ,x(:,:, 3) , T_x(:,:,4) ,x(:,:, 4) ) ;
 
figure (2)
plot (  T_x(:,:,1), x(:,:, 1) , T_x(:,:,2), x(:,:, 2) ,T_x(:,:,3) ,x(:,:, 3) , T_x(:,:,4) ,x(:,:, 4)  ,T_x(:,:,5) , x(:,:, 5) ) ;

% figure (1)
% plot (  T_x(:,:,1), x(:,:, 1) , T_x(:,:,2), x(:,:, 2) ,T_x(:,:,3) ,x(:,:, 3) , T_x(:,:,4) ,x(:,:, 4)  ,T_x(:,:,5) , x(:,:, 5) ,T_x(:,:,6) , x(:,:, 6) ) ;


