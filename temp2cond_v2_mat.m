%Nicole Bienert
%Purpose: Use temperature to calculate conductivity in S/m using the
%concentrations of H+ and salt from Siple Dome. The equation is based upon
%Matsuoka 2012

%Version History:
%This is a revised version of conductivityTempRelation.m
%ver1: takes a temperature array and imposes a constant salt and H+ concentration
%ver2: takes salt and H+ concentration arrays

%if using single H+ and salt concentrations, use: C1=1.2 uM; C2, from sea salt,=4.1 uM)

function  condTot = temp2cond_v2_mat(T,Hconcentration,saltConcentration)


T = T+273.15; %switch from C to K


Tr=251; %degreed kelvin for the recerence temperature
k = 0.00008617; %boltzman constant in eV/K

%types = ['pure ice';'acid';'salt']; %what each index corresponds to
cond=[9.2e-6; 3.2; 0.43];%molar conductivities in S/m/M
E=[0.51;0.2;0.19]; %activiation energy in eV

C_sipleDome(:,:,1)= ones(size(Hconcentration));
C_sipleDome(:,:,2)= Hconcentration;
C_sipleDome(:,:,3)= saltConcentration; %Siple Dome

%Ice at Siple Dome
condTot=zeros(size(Hconcentration));
for ind=1:length(cond)
    condTot=condTot+cond(ind).*C_sipleDome(:,:,ind).*exp(-E(ind)./k.*((1./T)-1/Tr));
end



% atten = 0.914*condTot*10^6;
% figure()
% hold on
% plot(T-273.15,atten)
% hXlabel = xlabel('Temperature (Deg C)')
% hYlabel = ylabel('Attenuation (dB/km)')
% Aesthetics_Script

end