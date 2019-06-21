function GreenFunctionErrorReceptorsComp(FileNameAnalytic, FileNameNum1,FileNameNum2)
               
load(FileNameNum1,'f','ExR','EzR','EyR');

EzRStd  = EzR;
ExRStd  = ExR;
EyRStd  = EyR;

indexFrequency = find(f(:,3));

load(FileNameNum2,'ExR','EzR','EyR');

EzROpt = EzR;
ExROpt = ExR;
EyROpt = EyR;

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
 
errPhaseStd_z = mod(atan2(imag(EzRStd),real(EzRStd)),2*pi) - mod(atan2(imag(EzAnaR),real(EzAnaR)),2*pi);
errPhaseOpt_z = mod(atan2(imag(EzROpt),real(EzROpt)),2*pi) - mod(atan2(imag(EzAnaR),real(EzAnaR)),2*pi);
errPhaseStd_z = ModPi(errPhaseStd_z)/pi*100;
errPhaseOpt_z = ModPi(errPhaseOpt_z)/pi*100;

errPhaseStd_x = mod(atan2(imag(ExRStd),real(ExRStd)),2*pi) - mod(atan2(imag(ExAnaR),real(ExAnaR)),2*pi);
errPhaseOpt_x = mod(atan2(imag(ExROpt),real(ExROpt)),2*pi) - mod(atan2(imag(ExAnaR),real(ExAnaR)),2*pi);
errPhaseStd_x = ModPi(errPhaseStd_x)/pi*100;
errPhaseOpt_x = ModPi(errPhaseOpt_x)/pi*100;

errPhaseStd_y = mod(atan2(imag(EyRStd),real(EyRStd)),2*pi) - mod(atan2(imag(EyAnaR),real(EyAnaR)),2*pi);
errPhaseOpt_y = mod(atan2(imag(EyROpt),real(EyROpt)),2*pi) - mod(atan2(imag(EyAnaR),real(EyAnaR)),2*pi);
errPhaseStd_y = ModPi(errPhaseStd_y)/pi*100;
errPhaseOpt_y = ModPi(errPhaseOpt_y)/pi*100;

errMagStd_z = (abs(EzRStd)-abs(EzAnaR))./abs(EzAnaR)*100;
errMagOpt_z = (abs(EzROpt)-abs(EzAnaR))./abs(EzAnaR)*100;
errMagStd_x = (abs(ExRStd)-abs(ExAnaR))./abs(ExAnaR)*100;
errMagOpt_x = (abs(ExROpt)-abs(ExAnaR))./abs(ExAnaR)*100;
errMagStd_y = (abs(EyRStd)-abs(EyAnaR))./abs(EyAnaR)*100;
errMagOpt_y = (abs(EyROpt)-abs(EyAnaR))./abs(EyAnaR)*100;



for n = 1:NbReceptors  
    if (flagEx(n)&&flagEz(n)&&flagEy(n))% 3 Components
        figure(n)
        subplot(2,1,1);
        plot(f(indexFrequency,1),errMagOpt_z(indexFrequency,n),'*b',f(indexFrequency,1),errMagStd_z(indexFrequency,n),'ob',...
            f(indexFrequency,1),errMagOpt_x(indexFrequency,n),'*r',f(indexFrequency,1),errMagStd_x(indexFrequency,n),'or',...
            f(indexFrequency,1),errMagOpt_y(indexFrequency,n),'*k',f(indexFrequency,1),errMagStd_y(indexFrequency,n),'ok');
        title({strcat('Green function error at receptor',num2str(n)),strcat('x = ',num2str(rx(n),4),'m; z = ',num2str(rz(n),4),'m')});
        legend('Z component, Optimal coefficients','Z component, Standard coefficients',...
                'X component, Optimal coefficients','X component, Standard coefficients',...
                'Y component, Optimal coefficients','Y component, Standard coefficients','Location','NorthWest');
        ylabel({'Error in magnitude';'of Green''s function';'(%)'});
        ax = gca;
        ax.YGrid = 'on';
        xlabel('Frequency (MHz)')
        
        subplot(2,1,2);
        plot(f(indexFrequency,1),errPhaseOpt_z(indexFrequency,n),'*b',f(indexFrequency,1),errPhaseStd_z(indexFrequency,n),'ob',...
                f(indexFrequency,1),errPhaseOpt_x(indexFrequency,n),'*r',f(indexFrequency,1),errPhaseStd_x(indexFrequency,n),'or',...
                f(indexFrequency,1),errPhaseOpt_y(indexFrequency,n),'*k',f(indexFrequency,1),errPhaseStd_y(indexFrequency,n),'ok');
        ylabel({'Error in phase';'of Green''s function';'(%)'});
        ax = gca;
        ax.YGrid = 'on';
        xlabel('Frequency (MHz)')
    
    elseif(flagEz(n)) % only z components
        figure(n)
        subplot(2,1,1);
        plot(f(indexFrequency,1),errMagOpt_z(indexFrequency,n),'*b',f(indexFrequency,1),errMagStd_z(indexFrequency,n),'ob');
        title({strcat('Green function error at receptor',num2str(n)),strcat('x = ',num2str(rx(n),4),'m; z = ',num2str(rz(n),4),'m')});
        legend('Z component, Optimal coefficients','Z component, Standard coefficients','Location','NorthWest');
        ylabel({'Error in magnitude';'of Green''s function';'(%)'});
        ax = gca;
        ax.YGrid = 'on';
        xlabel('Frequency (MHz)')
        
        subplot(2,1,2);
        plot(f(indexFrequency,1),errPhaseOpt_z(indexFrequency,n),'*b',f(indexFrequency,1),errPhaseStd_z(indexFrequency,n),'ob');
        ylabel({'Error in phase';'of Green''s function';'(%)'});
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