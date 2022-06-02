function mutated = ANN_mutate(ANN)
%ANN_MUTATE Takes an ANN, cell-matrix object and mutates the weights
%according to rules specified in paper, this results
%   Detailed explanation goes here
    mut_prob = 0.01; %0.01 used in real script 
    mut_std = 25; % standard deviation of mutations, 0.25 used in real

    weights = ANN;
    W1 = weights{1};
    b1 = weights{2};
    W2 = weights{3};
    b2 = weights{4};

    % NOTE WE HAVE TO CONSIDER SIGNS OF MUTATION, we also have to
    % use randn, note in this case we model using binomial distribution
    dW1 = (rand(size(W1)) < mut_prob).*randn(size(W1))*mut_std;
    dW2 = (rand(size(W2)) < mut_prob).*randn(size(W2))*mut_std;
    db1 = (rand(size(b1)) < mut_prob).*randn(size(b1))*mut_std;
    db2 = (rand(size(b2)) < mut_prob).*randn(size(b2))*mut_std;

    newW1 = W1 + dW1;
    newb1 = b1 + db1;
    newW2 = W2 + dW2;
    newb2 = b2 + db2;
    mutated = {newW1, newb1, newW2, newb2};

end

