clc
clear
FreqMax = pi; % PSD is plotted in the range [0, FreqMax]
OutputFolder = "./Results/Feedback"; % folder where output is stored (PadeDerivatives.txt and DirectG_Estimate.txt)
order = 2; % select order of the Pade Approximant
[NumCoeff, DenCoeff, ValidationScore] = ComputePadePSD(OutputFolder,order);
disp('The validation score is:')
ValidationScore
omega = 0:0.01:FreqMax;
Pade_PSD = zeros(length(omega),1);
for j=1:1:length(omega)
    w = omega(j);
    Pade_PSD(j) = PadeApproximationPSD(w, NumCoeff, DenCoeff);
end
figure;
plot(omega,Pade_PSD,'-','Color','r','LineWidth',3); hold on;
xlim([0 FreqMax]);
xlabel('\omega','FontSize',20) ;
ylabel('PSD','FontSize',20) ;

function psd = PadeApproximationPSD(omega, NumCoeff, DenCoeff)
    order = length(NumCoeff);
    x = 1i*omega;
    num = 0;
    den = x^order;
    for j=1:1:order
        num = num + NumCoeff(j)*(x^(order- j));
        den = den + DenCoeff(j+1)*(x^(order- j));
    end
    psd = 2*real(num/den);
end

function [NumCoeff, DenCoeff, ValidationScore] = ComputePadePSD(Outputfolder,order)
    FileName1 = Outputfolder + "/PadeDerivatives.txt";
    FileName2 = Outputfolder + "/DirectG_Estimate.txt";
    delimiterIn =' ';
    headerlinesIn = 1;
    
    Data_PadeDerivatives = importdata(FileName1,delimiterIn,headerlinesIn);
    Data_DirectG_Estimates = importdata(FileName2,delimiterIn,headerlinesIn);

    PadeDerivatives_Values = Data_PadeDerivatives.data(1:(2*order),2);
    DirectG_Estimates = Data_DirectG_Estimates.data(:,2);
    Svalues = Data_DirectG_Estimates.data(:,1);

    [NumCoeff, DenCoeff] = PadeApproximation(PadeDerivatives_Values);
    %NumCoeff = NumCoeff./PadeDerivatives_Values(1); %Normalization
    disp('Pade Approximation');
    vpa(poly2sym(NumCoeff)/poly2sym(DenCoeff),5)

    ValidationScore =1 - (1/length(Svalues))*norm(Computed_DRMS(Svalues,NumCoeff, DenCoeff) - DirectG_Estimates)/PadeDerivatives_Values(1);
    
    if(order == 2)
        delta_1 = (PadeDerivatives_Values(1)*PadeDerivatives_Values(4) - PadeDerivatives_Values(2)*PadeDerivatives_Values(3))/(PadeDerivatives_Values(1)*PadeDerivatives_Values(3) - PadeDerivatives_Values(2)^2 );
        delta_2 = (PadeDerivatives_Values(2)*PadeDerivatives_Values(4) - PadeDerivatives_Values(3)^2)/(PadeDerivatives_Values(1)*PadeDerivatives_Values(3) - PadeDerivatives_Values(2)^2 );
        beta = PadeDerivatives_Values(1)*delta_1/PadeDerivatives_Values(2) - 1;
        w_max_sq = beta*delta_2*(-1 + sqrt(1 + (beta^(-2))*(1 - beta*( delta_1^2 - 2*delta_2 )/delta_2   ) ));
        w_max = sqrt(w_max_sq) 
    end
    

end

function [NumCoeff, DenCoeff] = PadeApproximation(PadeDerivatives)

    D = PadeDerivatives;
    syms x;
    p = length(D)/2;

    C = sym(zeros(p,p+1));

    for i=1:1:p
        for j=1:1:(p+1)
            C(i,j) = D(i+j-1);
        end
    end
    L_row1 = sym(zeros(1,p+1));
    L_row2 = sym(zeros(1,p+1));

    for j=1:1:p
        for i=(p-j):1:(p-1)
            L_row1(j+1) = L_row1(j+1) + (x^i)*D(i+j-p+1);
            L_row2(j) = x^(p-j+1);
        end
    end
    L_row2(p+1) = 1;

    NumMatrix = [C;L_row1];
    DenMatrix = [C;L_row2];


    NumCoeff = flip(sym2poly(vpa(simplify(det(NumMatrix)),5)));
    DenCoeff = flip(sym2poly(vpa(simplify(det(DenMatrix)),5)));

    NumCoeff = NumCoeff/DenCoeff(1);
    DenCoeff = DenCoeff/DenCoeff(1);

end

function p = Computed_DRMS(Svalues,NumCoeff, DenCoeff)
    order = length(NumCoeff);
    p = zeros(length(Svalues),1);
    for k=1:1:length(Svalues)
        den = Svalues(k)^order;
        num = 0;
        for j=1:1:order
            num = num + NumCoeff(j)*(Svalues(k)^(order- j));
            den = den + DenCoeff(j+1)*(Svalues(k)^(order- j));
        end
        p(k) = num/den;
    end
end
 
