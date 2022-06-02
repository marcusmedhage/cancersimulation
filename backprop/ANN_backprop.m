clear all;
close all;

rng(1);

samples = 400000;
p = 4;
epochs = 10000;

train_x = rand(samples, p); % test x data

pH_max = 9;
pH_min = 5;
train_x(:,3) = train_x(:,3)*(pH_max - pH_min) + pH_min;
train_x(:,4) = round(train_x(:,4)*9); % sets neighbors to 0-4
%disp(test_x);


train_y_cat = cell_resp(train_x(:,1), train_x(:,2), train_x(:,3), train_x(:,4));
train_y_ohe = zeros(samples, 3);

for n = 1:samples
    
    if train_y_cat(n) == 1
        train_y_ohe(n,:) = [1 0 0]; 
    elseif train_y_cat(n) == 0
        train_y_ohe(n,:) = [0 1 0]; 
    elseif train_y_cat(n) == -1
        train_y_ohe(n,:) = [0 0 1]; 
    end
end


%disp(train_y);
% Featues in order Oxygen, GLucose, pH, Neighbors
test_x = rand(samples, p); % test x data
test_x(:,3) = round(test_x(:,3)*4); % sets neighbors to 0-4
%disp(test_x);
test_y = cell_resp(test_x(:,1), test_x(:,2), test_x(:,3));
%disp(train_y);


% NETWORK ARCHITECTURE
p = 4; % features
U1 = 5; % hidden layer size
U2 = 3; % outputs
trainFunc = 'trainlm';

%net = fitnet(U1,trainFunc);



% NETWORK INITIATION AND TRAINING
net = patternnet(U1, trainFunc);
init(net);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio =   15/100;
net.divideParam.testRatio = 15/100;
net.trainParam.epochs = epochs;
%net.performParam.regularization = 0.90;


%net.inputs{1}.processFcns = 
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'logsig';

disp("Starting training");
[net,tr] = train(net, train_x', train_y_ohe');
% NETWORK TESTING

pred_y = sim(net, test_x');
%performance = perform(net, test_y, pred_y);

disp("TRAINING COMPLETE");

    

W1 = net.IW{1};
b1 = net.b{1};
W2 = net.LW{2,1};
b2 = net.b{2};


W1_m = W1;
W2_m = [W2; 0 0 0 0 0];
b1_m = b1;
b2_m = [b2; 0];



