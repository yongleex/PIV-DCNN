function Y = euclideanloss(X, c, dzdy)
%EUCLIDEANLOSS Summary of this function goes here
%   Detailed explanation goes here
c= single(c);
assert(numel(X) == numel(c));

%X = X + 1e-6 ;
sz = [size(X,1) size(X,2) size(X,3) size(X,4)] ;

if numel(c) == sz(4)
  % one label per image
  c = reshape(c, [1 1 1 sz(4)]) ;
end
if size(X,3) == size(c,1)
    % multiple output
   c = reshape(c, [1 1 sz(3) sz(4)]) ;
end



d = size(X);

assert(all(d == size(c)));

if nargin == 2 || (nargin == 3 && isempty(dzdy))
    
    Y = 1 / 2 * sum(subsref((X - c) .^ 2, substruct('()', {':'}))); % Y is divided by d(4) in cnn_train.m / cnn_train_mgpu.m.
%     Y = 1 / (2 * prod(d(1 : 3))) * sum(subsref((X - c) .^ 2, substruct('()', {':'}))); % Should Y be divided by prod(d(1 : 3))? It depends on the learning rate.
    
elseif nargin == 3 && ~isempty(dzdy)
    
    assert(numel(dzdy) == 1);
    
    Y = dzdy * (X - c); % Y is divided by d(4) in cnn_train.m / cnn_train_mgpu.m.
%     Y = dzdy / prod(d(1 : 3)) * (X - c); % Should Y be divided by prod(d(1 : 3))? It depends on the learning rate.
    
end

end