% gamma plus
g1=5;
% gamma minus
g2=4;

h = 10000;

H = [[0, -1i/2*h];
    [1i/2*h, 0]];

L = [[1/2*sqrt((g1+g2)/2), -1i/2*sqrt((g1+g2)/2)];
    [1/2*sqrt((g1-g2)/2), 1i/2*sqrt((g1-g2)/2)]];

R = [[0, 0, 1, 0];
    [0, 0, 0, 1];
    [-1, 0, 0, 0];
    [0, -1, 0, 0]];

B = [[0, 1i/2*g1, -1i/2*g2, 1/2*g2];
    [-1i/2*g1, 0, 1/2*g2, 1i/2*g2];
    [1i/2*g2, -1/2*g2, 0, 1i/2*g1];
    [-1/2*g2, -1i/2*g2, -1i/2*g1, 0]];

%A = -h*R+B;

M = compute_M(L);
A = compute_A(H, M);

[n, ~] = size(A);
n = n/2;
[eigenvectors, eigenvalues] = eig(A);
eigenvalues = diag(eigenvalues);

[V, num_degen_eigenval] = sort_eigenvalues( eigenvectors, eigenvalues);
% [~, index] = sort(real(eigenvalues));
% V = zeros(2*n, 2*n);
% for i=1:n
%    V(2*i-1, :) = eigenvectors(:, index(i)); % neg eig
%    V(2*i, :) = eigenvectors(:, index(2*n+1-i)); % pos eig
% end

% number of different eigenvalues (up to sign)
num_blocks = size(num_degen_eigenval, 2);
processed_eigenvectors = 0;
for block=1:num_blocks
   num_eigenval = num_degen_eigenval(block);
   T = zeros(num_eigenval/2, num_eigenval/2);
   for i=1:num_eigenval/2
       for j=1:num_eigenval/2
           T(i, j) = V(processed_eigenvectors+2*i, :)*V(processed_eigenvectors+2*j-1, :).';
       end
   end
   U = (inv(T)).';
   V_copy = V;
   for i=1:num_eigenval/2
       V(processed_eigenvectors+2*i-1, :) = 0;
       for j=1:num_eigenval/2
           V(processed_eigenvectors+2*i-1, :) = V(processed_eigenvectors+2*i-1, :) + U(i, j)*V_copy(processed_eigenvectors+2*j-1, :);
       end
   end
   processed_eigenvectors = processed_eigenvectors + num_eigenval;
end

% V = [[g2/g1-1, 1i*(g2/g1+1), -1i*(g2/g1-1), g2/g1+1];
%     [-1/4, -1i/4, -1i/4, 1/4];
%     [g2/g1+1, 1i*(g2/g1-1), -1i*(-g2/g1-1), -g2/g1+1];
%     [1/4, 1i/4, -1i/4, 1/4]];

%w1w2 = 1/2*(V(2, 1)*V(1, 3)-V(2, 2)*V(1, 4)-1i*V(2, 2)*V(1, 3)-1i*V(2, 1)*V(1, 4)+V(4, 1)*V(3, 3)-V(4, 2)*V(3, 4)-1i*V(4, 2)*V(3, 3)-1i*V(4, 1)*V(3, 4));
w1w2 = 0;
for m=1:n
   w1w2 = w1w2 + 4*V(2*m, 1)*V(2*m-1, 3);% - V(2*m, 2)*V(2*m-1, 4) - 1i*V(2*m, 2)*V(2*m-1, 3) - 1i*V(2*m, 1)*V(2*m-1, 4); 
end

w1w2 = w1w2/2;

disp(1/2-1j/2*w1w2);
disp((g1+g2)/(2*g1));

[n, ~] = size(V);
N = zeros(n, n);
for i=1:n
    for j=1:n
        N(i, j) = V(i, :)*V(j, :).';
    end
end
%disp(N);