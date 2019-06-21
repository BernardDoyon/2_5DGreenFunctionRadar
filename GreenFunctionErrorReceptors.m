function GreenFunctionErrorReceptors(FileNameAnalytic, FileNameNum)
               
load(FileNameNum,'f','ExR','EzR','EyR');

EzNumR = EzR;
ExNumR = ExR;
EyNumR = EyR;

indexFrequency = find(f(:,3));

load(FileNameAnalytic,'ExR','EzR','EyR','f','q','sx','sz','rx','rz','ry');

EzAnaR = EzR;
ExAnaR = ExR;
EyAnaR = EyR;

NbReceptors = length(rx);
%check if analytical solution is available
flagEx = zeros(1,NbReceptors);
flagEz = zeros(1,NbReceptors);
flagEy = zeros(1,NbReceptors);

for n = 1:NbReceptors
   if (isnan(ExAnaR(:,n)))
       flagEx(1,n) = 0;
   else
       flagEx(1,n) = 1;
   end
end

for n = 1:NbReceptors
   if (isnan(EzAnaR(:,n)))
       flagEz(1,n) = 0;
   else
       flagEz(1,n) = 1;
   end
end

for n = 1:NbReceptors
   if (isnan(EyAnaR(:,n)))
       flagEy(1,n) = 0;
   else
       flagEy(1,n) = 1;
   end
end

errPhase_z = mod(atan2(imag(EzNumR),real(EzNumR)),2*pi) - mod(atan2(imag(EzAnaR),real(EzAnaR)),2*pi);
errPhase_z = ModPi(errPhase_z)/pi*100;
errPhase_x = mod(atan2(imag(ExNumR),real(ExNumR)),2*pi) - mod(atan2(imag(ExAnaR),real(ExAnaR)),2*pi);
errPhase_x = ModPi(errPhase_x)/pi*100;
errPhase_y = mod(atan2(imag(EyNumR),real(EyNumR)),2*pi) - mod(atan2(imag(EyAnaR),real(EyAnaR)),2*pi);
errPhase_y = ModPi(errPhase_y)/pi*100;

errMag_z = (abs(EzNumR)-abs(EzAnaR))./abs(EzAnaR)*100;
errMag_x = (abs(ExNumR)-abs(ExAnaR))./abs(ExAnaR)*100;
errMag_y = (abs(EyNumR)-abs(EyAnaR))./abs(EyAnaR)*100;

for n = 1:NbReceptors
    % 3 components (Homogeneous model)
    if (flagEx(n)&&flagEz(n)&&flagEy(n))
        figure(n)
        subplot(2,1,1);
        plot(f(indexFrequency,1),errMag_z(indexFrequency,n),'*b',...
                    f(indexFrequency,1),errMag_x(indexFrequency,n),'or',...
                    f(indexFrequency,1),errMag_y(indexFrequency,n),'xk');
        title({strcat('Green function error at receptor',...
                   num2str(n)),strcat('x = ',num2str(rx(n),4),'m; z = ',num2str(rz(n),4),'m')});
        legend('z component','x component','y component','Location','NorthWest');
        % ylim([-2 2]);
        ylabel({'Error in magnitude';'of Green''s function';'(%)'});
        ax = gca;
        ax.YGrid = 'on';
        xlabel('Frequency (MHz)')
        
        subplot(2,1,2);
        plot(f(indexFrequency,1),errPhase_z(indexFrequency,n),'*b',...
            f(indexFrequency,1),errPhase_x(indexFrequency,n),'or',...
            f(indexFrequency,1),errPhase_y(indexFrequency,n),'xk');
        ylabel({'Error in phase';'of Green''s function';'(%)'});
        %ylim([-1 4]);
        ax = gca;
        ax.YGrid = 'on';
        xlabel('Frequency (MHz)')
        
    elseif(flagEz(n)) % only z components
        figure(n)
        subplot(2,1,1);
        plot(f(indexFrequency,1),errMag_z(indexFrequency,n),'*b');
        title({strcat('Green function error at receptor',num2str(n)),strcat('x = ',num2str(rx(n),4),'m; z = ',num2str(rz(n),4),'m')});
        legend('z component','Location','NorthWest');
        % ylim([-2 2]);
        ylabel({'Error in magnitude';'of Green''s function';'(%)'});
        ax = gca;
        ax.YGrid = 'on';
        xlabel('Frequency (MHz)')
        
        subplot(2,1,2);
        plot(f(indexFrequency,1),errPhase_z(indexFrequency,n),'*b');
        ylabel({'Error in phase';'of Green''s function';'(%)'});
        %ylim([-1 4]);
        ax = gca;
        ax.YGrid = 'on';
        xlabel('Frequency (MHz)') 
    end

end

end

function lambda = ModPi(lambda)
s =size(lambda);

for i =1:s(1)
    for j = 1:s(2)
        if (lambda(i,j)>pi)
            lambda(i,j) = 2*pi-lambda(i,j);
        end
        if (lambda(i,j)<-pi)
            lambda(i,j) = 2*pi+lambda(i,j);
        end
    end
    
end
end