function [grads,idx,idx2] = remove_duplicated_gradients(grads)
% function [grads,idx,idx2] = remove_duplicated_gradients(grads)
%
%   Removes duplicated gradients from a list grads with size 2Gx3, where it
%   is assumed that for each gi(n,:) there is a gi(m,:)=-gi(n,:). The
%   function returns grads with size Gx3 without their opposites, idx with
%   size 1xG with the indices of the selected gradients among the original
%   set annd idx2 with size 1xG with the indices of their corresponding
%   opposites within the original volume.

moduli = sqrt(sum(grads.*grads,2));
grads  = grads./repmat(moduli,[1,3]);
% Find dot-products of each gradient with the remaining ones:
dp   = grads(:,1)*(grads(:,1)') + grads(:,2)*(grads(:,2)') + grads(:,3)*(grads(:,3)');
dp   = abs(dp); % The closer to 1, the more similar
dp   = dp - diag(diag(dp)); % Avoid self-similarity
idx  = 1:size(grads,1);
idx2 = 1:(size(grads,1)/2);
for n=1:size(grads,1)/2
    pos = idx(n);
    [~,mpos] = max(dp(pos,:)); % mpos is the most similar gradient to n, so...
    idx = setdiff(idx,mpos);   % ... remove it from the list
    idx2(n) = mpos;            % store the most similar gradient to idx(n)
end
grads = grads(idx,:);