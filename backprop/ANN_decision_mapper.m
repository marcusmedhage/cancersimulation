clear all;
close all;

samples = 100000;
p = 4;

ANN = ANN_loader_pressure();

glucose = 0.7;
pH = 8;

features = ["Oxygen", "Glucose", "pH", "Neighbors"]; % Age eventually out

test_x = rand(samples, p); % test x data
test_x(:,2) = glucose;

pH_max = 9;
pH_min = 6;
n_max = 9;
test_x(:,3) = pH; %test_x(:,3)*(pH_max - pH_min) + pH_min;
test_x(:,4) = (test_x(:,4)*n_max); % sets neighbors to 0-4

test_y = cell_resp(test_x(:,1), test_x(:,2), test_x(:,3), test_x(:,4));

preds = zeros(samples,1);

for n = 1:samples
   preds(n) = ANN_predict(ANN, test_x(n,:)'); % prediction on sample 
end

figure(1);

gscatter(test_x(:,1), test_x(:,4), preds);
xlabel("Oxygen");
ylabel("Neighbors");
legend("Apoptosis", "Quiescence", "Proliferate");
title("Cell response from trained neural network");

colors = zeros(samples,1);


%figure(2);
%scatter3(test_x(:,1), test_x(:,2) ,test_x(:,3), 20, preds, 'filled');
%legend(preds);
%xlabel("Oxygen");
%ylabel("Glucose");