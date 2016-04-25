function RM = vrot3( M, r, alpha )
%VROT Rotate an array of vectors about an axis by alpha
%   M: array of (3D) vectors
%   r: axis of rotation (3D) vector. Need not be normalized. Can also be
%       an Nx3 array of rotation vectors.
%   alpha: angle of rotation (rad). Can also be an N-dim array.
RM = [];
[N,col] = size(M);
if col ~= 3
    error('M must by Nx3.'); 
end

[row,col] = size(r);
if col ~= 3
    error('r must by either 1x3 or Nx3.');
end
if row == 1
    r = repmat(r,N,1);
end
[row,col] = size(r);
if row ~= N
    error('r must by either 1x3 or Nx3.');
end

[row,col] = size(alpha);
if col ~= 1
    error('alpha must be either 1x1 or Nx1.');
end
if row == 1
    alpha = repmat(alpha,N,1);
end
[row,col] = size(alpha);
if row ~= N
    error('alpha must be either 1x1 or Nx1.');
end

for i = 1:N
    r(i,:) = r(i,:)/norm(r(i,:));
end
ca = cos(alpha);
sa = sin(alpha);

for i = 1:N
    R = [...
        ca(i)+r(i,1)^2*(1-ca(i)),...
        r(i,1)*r(i,2)*(1-ca(i))-r(i,3)*sa(i),...
        r(i,1)*r(i,3)*(1-ca(i))+r(i,2)*sa(i)
        r(i,2)*r(i,1)*(1-ca(i))+r(i,3)*sa(i),...
        ca(i)+r(i,2)^2*(1-ca(i)),...
        r(i,2)*r(i,3)*(1-ca(i))-r(i,1)*sa(i)
        r(i,3)*r(i,1)*(1-ca(i))-r(i,2)*sa(i),...
        r(i,3)*r(i,2)*(1-ca(i))+r(i,1)*sa(i),...
        ca(i)+r(i,3)^2*(1-ca(i)) ];
    R = transpose(R);
    RM = [ RM; M(i,:)*R ];
end

end

