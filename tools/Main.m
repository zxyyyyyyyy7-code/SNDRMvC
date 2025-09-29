clc
clear;
% 参数数组
Eta=[0.1,0.01,0.001,0.0001,0.00001,0.000001];
Lambda=[0.5,0.1,0.05,0.01,0.005,0.0001];
Beta1=[0.5,0.1,0.05,0.01,0.005,0.0001];
Beta2=[0.5,0.1,0.05,0.01,0.005,0.0001];
Beta3=[0.5,0.1,0.05,0.01,0.005,0.0001];
for eta = Eta
    for lambda = Lambda
        for beta1 = Beta1

            for beta2=Beta2
                for beta3 = Beta3
                    fprintf('Eta: %.6f, Lambda: %.6f, Beta1: %.6f, Beta3: %.6f\n', ...
                        eta, lambda, beta1, beta3);
                    
                end
            end
        end
    end
end