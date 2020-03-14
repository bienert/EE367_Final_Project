%Nicole Bienert
%Purpose: Use temperature to calculate conductivity in S/m using the
%concentrations of H+ and salt from Siple Dome. The equation is based upon
%Matsuoka 2012. 
% Use Gradient Descent to solve

%Version History:
%This is a revised version of conductivityTempRelation.m
%ver1: takes a temperature array and imposes a constant salt and H+ concentration
%ver2: takes salt and H+ concentration arrays
%ver3: accepts matricies and reconstructs the temperature using gradient
%descent

%if using single H+ and salt concentrations, use: C1=1.2 uM; C2, from sea salt,=4.1 uM)

function  T = cond2temp_v3_mat(cond,Hconcentration,saltConcentration)

Tr=251; %degreed kelvin for the recerence temperature
k = 0.00008617; %boltzman constant in eV/K

%types = ['pure ice';'acid';'salt']; %what each index corresponds to
condMolar=[9.2e-6, 3.2, 0.43];%molar conductivities in S/m/M
E=[0.51,0.2,0.19]; %activiation energy in eV

%reshape matrix into array
C_sipleDome(:,1)= ones(prod(size(Hconcentration)),1);
C_sipleDome(:,2)=reshape(Hconcentration,[prod(size(Hconcentration)),1]);
C_sipleDome(:,3)=reshape(saltConcentration,[prod(size(Hconcentration)),1]);
condArray=reshape(cond,[prod(size(Hconcentration)),1]);

%% gradient descent
Tguess=-13.2;
stepSize=1;
tol=1e-12;
maxItter=200;
figure()


for ind = 1: length(condArray)
    iter=2;
    residual=300;
     while iter<maxItter & residual(iter-1)>tol
         %compute conductivity
        condNow(iter)=sum(condMolar.*C_sipleDome(ind,:).*exp(-E./k.*(1./(Tguess(iter-1)+273.15)-1/Tr)),2);
        %calculate residual
        residual(iter)=sqrt(abs(condNow(iter)^2-condArray(ind)^2));
        %update step size 
        stepSize(iter)=(condNow(iter)-condArray(ind))*10^4*-1; %subtract derivative
        %update temperature
         Tguess(iter)=Tguess(iter-1)+stepSize(iter);

        
%         stepSize(iter)
%         subplot(2,2,1)
%         plot(residual)
%         title('residual')
%         subplot(2,2,2)
%         plot(Tguess)
%         title('T')
%         subplot(2,2,3)
%         plot(stepSize)
%         title('Step Size')
%         subplot(2,2,4)
% pause(0.01)
        
        iter=iter+1;
     end
    T(ind)=Tguess(iter-1);
end
  %%  

T=reshape(T,size(cond));

end