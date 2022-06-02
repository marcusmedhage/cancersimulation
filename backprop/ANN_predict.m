function [prediction, metabolism] = ANN_predict(weights, x)
%ANN_PREDICT Takes a weights object and inputs and returns prediction
%   weights - A cell object with dense matrices W1,b1,W2,b2 used for multiplication of
%   neural network output
%   x - single line input of neural network, should be p x 1 shape matrix
% 
%   Returns prediction of response,
%   1   = Proliferate, [1 0 0]
%   0   = Quiescent, [0 1 0]
%   -1  = Apoptosis, [0 0 1]
    %ANN = {W1, b1, W2, b2}; % Matlab-cell object of neural network matrices
    ANN = weights;
    %disp("X");
    %disp(x);
    
    x = inputProcess(x);
    
    for layer=1:(length(ANN)/2)
        W = ANN{2*layer -1}; % weight matrix
        b = ANN{2*layer}; %bias
        
     %   disp(2*layer-1);
     %   disp(2*layer);
%         disp("NEW");
%         disp(x);
%         disp("W");
%         disp(W);
%         disp(b);
        
        q = sigmoid(W*x + b); % x must be transpose [p x 1] vector needed
        x = q;
        %disp(x);
    end % end loop
    %disp(x);
    % Evaluates highest response and returns
    [maxval, maxind] = max(q);
    %disp("OUTPUT");
    %disp(q);
    metabolism = maxval;
    prediction = -(maxind - 2);

    
end

