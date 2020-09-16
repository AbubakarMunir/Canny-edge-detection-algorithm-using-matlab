classdef CannyAlgo
    
    methods(Static)
        
        function [sum] = prodsum(mat1,mat2)
          s=size(mat1,1);
          mat=zeros(s,s);
          for i=1:s
              for j=1:s
                mat(i,j)=mat1(i,j)*mat2(i,j);
              end

          end
          sum=0;
          for i=1:s
              for j=1:s
                sum=sum+mat(i,j);
              end
          end

          

        end

  

        
        %Function to take out smaller matrices out from larger 0 padded
        %matrix for each iteration of convolution
        
        function [r] = takeoutMatrices(mat,kernel,a,b)
        s=size(kernel,1);
        r=zeros(s,s);
        c=1;
        d=1;
        for i=a:(a+s-1)
            for j=b:(b+s-1)
              r(c,d)=mat(i,j);
              d=d+1;
            end
            d=1;
            c=c+1;
        end
        end
        
        
        
        
        
        % Final Function for convolution and it uses above two functions
        function [res] = convulee(mat,kernel)
                    r=size(mat,1)+2;
                    c=size(mat,2)+2;
                    
                    res=zeros(r-2,c-2);
                    I=zeros(r,c);          %for mat with 0 padding 
                    R=zeros(r,c);          %resultant of convolution
                    k=1;
                    l=1;
                    for i=2:r-1
                        for j=2:c-1
                            I(i,j)=mat(k,l);
                            l=l+1;
                        end
                        k=k+1;
                        l=1;
                    end
                    
                    ks=size(kernel,1);
                    for i=1:(r-(ks-1))
                        for j=1:(c-(ks-1))
                            [s] = CannyAlgo.takeoutMatrices(I,kernel,i,j);
                            [prod] = CannyAlgo.prodsum(s,kernel);
                            R(i+1,j+1)=prod;
                            
                        end
                        
                    end
                   
                    k=1;
                    l=1;
                    for i=2:r-1
                        for j=2:c-1
                            res(k,l)=R(i,j);
                            l=l+1;
                        end
                        k=k+1;
                        l=1;
                    end
                    
        end
        
        function [R] = convolution(I,kernel)
            thirdD=size(I,3);
            %Translation for 2D greysclae image
            if(thirdD==1)        
                R=convulation.convulee(I,kernel);
            %Translation for 3D colored image
            else
                r=I(:,:,1);
                g=I(:,:,2);
                b=I(:,:,3);
                r=convulation.convulee(r,kernel);
                g=convulation.convulee(g,kernel);
                b=convulation.convulee(b,kernel);
                R(:,:,1)=r;
                R(:,:,2)=g;
                R(:,:,3)=b;
                
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Laplacian of gaussian mask applied%%%%%%%%%%%%%%%%%%%%%
        function R = LOG(I)
           
            
            log=[1/16  2/16 1/16; 2/16 4/16 2/16; 1/16  2/16 1/16];
            R=rgb2gray(I);     %convert image to greyscale image 
           
            R=CannyAlgo.convulee(R,log);        % for smoothening
            R=uint8(R);
        end
        
        function [mag,deg]=add(I)
            verticalMask=[-1 -2 -1;0 0 0;1 2 1];
            horizontalMask=[-1 0 1;-2 0 2;-1 0 1];
           
            s1=CannyAlgo.convulee(I,horizontalMask);
            
            s2=CannyAlgo.convulee(I,verticalMask);
           
            mag=zeros(size(s1,1),size(s1,2));
            for i=1:size(s1,1)
                for j=1:size(s1,2)
                    mag(i,j)=sqrt((s1(i,j)^2 + s2(i,j)^2));
                    mag(i,j)=round(mag(i,j));
                end
            end
            deg=zeros(size(s1,1),size(s1,2));
            for i=1:size(s1,1)
                for j=1:size(s1,2)
                    deg(i,j)=atan2(s2(i,j),s1(i,j));
                    
                end
            end
            
            deg=deg*(180/pi);
            for i=1:size(deg,1)
                for j=1:size(deg,2)
                    if (deg(i,j)<0) 
                        deg(i,j)=360+deg(i,j);
                    end
                end
            end
            
            
            for i = 1  : size(deg,1)
                for j = 1 : size(deg,2)
                    if ((deg(i, j) >= 0 ) && (deg(i, j) < 22.5) || (deg(i, j) >= 157.5) && (deg(i, j) < 202.5) || (deg(i, j) >= 337.5) && (deg(i, j) <= 360))
                        deg(i, j) = 0;
                    elseif ((deg(i, j) >= 22.5) && (deg(i, j) < 67.5) || (deg(i, j) >= 202.5) && (deg(i, j) < 247.5))
                        deg(i, j) = 45;
                    elseif ((deg(i, j) >= 67.5 && deg(i, j) < 112.5) || (deg(i, j) >= 247.5 && deg(i, j) < 292.5))
                        deg(i, j) = 90;
                    elseif ((deg(i, j) >= 112.5 && deg(i, j) < 157.5) || (deg(i, j) >= 292.5 && deg(i, j) < 337.5))
                        deg(i, j) = 135;
                    end
                end
            end
            
            
            
            
            
            
        end
        
         function R = NonMaxSuppression(mag,dir)
             R=zeros(size(mag,1),size(mag,2));
             for i=1:(size(R,1))
                 for j=1:(size(R,2))
                     R(i,j)=mag(i,j);
                 end
             end
             for i=2:(size(R,1)-1)
                 for j=2:(size(R,2)-1)
                     if(dir(i,j)==45)
                         if(mag(i,j)<mag(i-1,j+1) && mag(i,j)<mag(i+1,j-1))
                             R(i,j)=0;
                         end
                     end
                     if(dir(i,j)==90)
                         if(mag(i,j)<mag(i+1,j) && mag(i,j)<mag(i-1,j))
                           R(i,j)=0;
                         end
                     end
                     
                     if(dir(i,j)==135)
                         if(mag(i,j)<mag(i+1,j+1) && mag(i,j)<mag(i-1,j-1))
                           R(i,j)=0;
                         end
                     end
                     
                     if(dir(i,j)==0)
                         if(mag(i,j)<mag(i,j+1) && mag(i,j)<mag(i,j-1))
                           R(i,j)=0;
                         end
                     end
                     
                     
                     
                 end
             end
         end
        
           
         
         
         function R=HysteresisThresholding(I)
             lt=0.075*(max(max(I)));
             ht=0.175*(max(max(I)));
             R=zeros(size(I,1),size(I,2));
             for i = 1  : size(I,1)
                for j = 1 : size(I,2)
                    if (I(i, j) < lt)
                        R(i, j) = 0;
                    elseif (I(i, j) > ht)
                        R(i, j) = 255;
                    
                    end
                end
             end
             
         end
        
        
        
%%%%%%%%%%%%%%%%%%%%% Final Canny function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555       
        
        function R = canny(I)
           if(size(I,3)>1)
               I=CannyAlgo.LOG(I);
           end
           [m,d]=CannyAlgo.add(I);
           R=CannyAlgo.NonMaxSuppression(m,d);
           R=CannyAlgo.HysteresisThresholding(R);    
        end
        
       
    end
end

