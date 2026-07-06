function ordered_evals= evals_ordering(HTM_evals)

%Matlab does not compute the eigenvalues following a consistent order
%during each iteration. 

%This function groups or match the eigenvalues generated in each iteration. 
% Then, this allows us to properly track the evolution of the each eigenvalue during each iteration. 
% It uses the 'matchpairs' function built in Matlab (it solves the linear assignment problem by minimizing the total cost) 

%Input: HTM_eval - set of eigenvalues computed each iteration. 
%Output: ordered_evals - set of ordered eigenvalues based on the initial set. 


ordered_evals = zeros(size(HTM_evals)); %Output matrix for the ordered set of eigenvalues
ordered_evals (:,1) =  HTM_evals(:,1);      %Saving the first set of evals. 


prev_eval = ordered_evals(:, 1);

xo = real (prev_eval);
yo = imag(prev_eval);


[nrows, ncolums] = size (HTM_evals);


for k = 2: ncolums

    track_evals = HTM_evals(:, k); 
    x = real(track_evals);
    y = imag(track_evals);

    %Define the cost matrix
    cost= zeros(nrows);

    for j = 1: nrows
        cost (j, :) = vecnorm([x(j)-xo'; y(j)-yo'] , 2,1)';
    end 

    M = matchpairs (cost, 1e4);
    index = M (:,1);


    reorder_evals = zeros(size (track_evals)); 


    for p  = 1: length(index(:,1))
        curr_idx = index (p, 1);
        reorder_evals (p) = track_evals(curr_idx); 
    end 

    prev_eval = reorder_evals;
    xo = real (prev_eval);
    yo = imag(prev_eval);

    ordered_evals (:, k) = reorder_evals;
end 
end