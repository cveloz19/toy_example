function [ordered_evals, ordered_evecs]= evals_and_evec_ordering(HTM_evals, HTM_evec)

%Matlab does not compute the eigenvalues following a consistent order
%during each iteration. 

%This function groups or match the eigenvalues generated in each iteration. 
% Then, this allows us to properly track the evolution of the each eigenvalue during each iteration. 
% It uses the 'matchpairs' function built in Matlab (it solves the linear assignment problem by minimizing the total cost) 

%Input: 
% HTM_eval - set of eigenvalues computed each iteration.
% HTM_evec - set of eigenvectors computed each iteration.

%Output: 
% ordered_evals - set of ordered eigenvalues based on the initial set. 
% ordered_evals - set of ordered eigenvectors based on the initial set. 



%Note: This function order the eigenvectors based on the proximity of
%eigenvalues calculated by matchpair function. 



%% START

% Ordered evals
ordered_evals = zeros(size(HTM_evals)); %Output matrix for the ordered set of eigenvalues
ordered_evals (:,1) =  HTM_evals(:,1);      %Saving the first set of evals. 

% Ordered evecs 
ordered_evecs = cell(length(HTM_evec), 1); %Output cell structure for the ordered set of eigenvalues
ordered_evecs {1,1} = cell2mat(HTM_evec{1,1});         % Saving the first set of eigenvectors 




prev_eval = ordered_evals(:, 1);

xo = real (prev_eval);
yo = imag(prev_eval);


[nrows, ncolums] = size (HTM_evals);






for k = 2: ncolums

    track_evals = HTM_evals(:, k); 
    track_evecs = HTM_evec {k,1};
    track_evecs = cell2mat(track_evecs);

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
    reorder_evecs = zeros(size (track_evecs)); 


    for p  = 1: length(index(:,1))
        %Ordering evals
        curr_idx = index (p, 1);
        reorder_evals (p) = track_evals(curr_idx);

        %Ordering evecs
        reorder_evecs (:,p) = track_evecs (:,curr_idx);
    end 



    prev_eval = reorder_evals;



    % Defining values for the cost function
    xo = real (prev_eval);
    yo = imag(prev_eval);

    ordered_evals (:, k) = reorder_evals;
    ordered_evecs {k} = reorder_evecs; 


end 
end