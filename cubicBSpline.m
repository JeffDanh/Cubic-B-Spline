% Jeff Duong

% Cubic B-spline Interpolation

fprintf(1,'------ Cubic B-spline Interpolation ------\n');

% input n, number of nodes
fprintf(1,'Enter n, number of nodes, for x(1), x(2), ... , x(n): \n');
numberOfNodes = input(' ');
 
% We use zeros() function to make our matrix
% zeros(m,n) gives a m-by-n matrix of zeros
% x is our matrix of 'x' node values
x_values = zeros(1,numberOfNodes);
% From Crout's reduction we create matrix 'A'
% A_coefficient is our values for f(x) or "y" values
A_coefficient = zeros(1,numberOfNodes);

% for loop will take input of each node/point, like 1st 'x', 2nd 'x', etc..
% 'i' will be our index or iteration. 
for i = 1:numberOfNodes 
   % formatting input
   fprintf(1,'Input x(%d), press enter. Then input f(x(%d)), press enter: \n', i, i);
   
   % adding input nodes into a matrix for later use
   % our x_n nodes will go into matrix x
   % x_values(i) indicates our position. We start at i because matrix position starts from 1
   x_values(i) = input(' ');
   % our f(x_n) or "y" value will go into matrix A_coefficient
   A_coefficient(i) = input(' ');
end


% 'h' is the difference or increment(step size) between each 'x' value. We create a
% matrix and store those values
h_stepSize = zeros(1,numberOfNodes-1);
% iterate through h
for i = 1:numberOfNodes-1
   % get difference between each x node
   h_stepSize(i) = x_values(i+1) - x_values(i);
end

% subinterval values for tridiagonal system
subintervalValues = zeros(1,numberOfNodes-1);

for i = 1:numberOfNodes-2
   % process of tridiagonal system
   subintervalValues(i+1) = 3*(A_coefficient(i+2)*h_stepSize(i) - A_coefficient(i+1)*(x_values(i+2)-x_values(i)) + A_coefficient(i)*h_stepSize(i+1))/(h_stepSize(i+1)*h_stepSize(i));
end

% creating section of L from LU decomposition in Crout factorization
L = zeros(1,numberOfNodes);

% upper values from LU decomposition in Crout factorization
U = zeros(1,numberOfNodes);

% intermediate values 
y_intermediateValues = zeros(1,numberOfNodes);

% The lower matrix values, L, has underlying ones. The first value in the lower is 1
L(1) = 1;
% upper matrix values, U, has underlying zeros. The first value in the upper is 0
U(1) = 0;
% If upper matrix first value is 0, then the intermediate value is 0
y_intermediateValues(1) = 0;

% Fill in values after the ones and zeros in L and U
for i = 1:numberOfNodes-2
   % filling in values from second up to last in matrices
   L(i+1) = 2*(x_values(i+2) - x_values(i)) - h_stepSize(i)*U(i);
   U(i+1) = h_stepSize(i+1)/L(i+1);
   y_intermediateValues(i+1) = (subintervalValues(i+1) - h_stepSize(i)*y_intermediateValues(i))/L(i+1);
end

% The last value in the lower is 1
L(numberOfNodes) = 1;

% The upper matrix has ended so the intermediate value is 0
y_intermediateValues(numberOfNodes) = 0;

% coefficient values for interpolant
% create B, C, D coefficients
B_coefficient = zeros(1,numberOfNodes);
C_coefficient = zeros(1,numberOfNodes);
D_coefficient = zeros(1,numberOfNodes);
% C coefficient was already calculated during intermediate calculation
C_coefficient(numberOfNodes) = y_intermediateValues(numberOfNodes);

% Calculating coefficients B, C, D 
for i = 1:numberOfNodes-1
   % Secondary counter to access index before i
   j = numberOfNodes-1-i;
   % Use coefficient equations from conditions of cubic spline
   C_coefficient(j+1) = y_intermediateValues(j+1) - U(j+1)*C_coefficient(j+2);
   B_coefficient(j+1) = (A_coefficient(j+2) - A_coefficient(j+1))/h_stepSize(j+1) - h_stepSize(j+1)*(C_coefficient(j+2) + 2*C_coefficient(j+1))/3;
   D_coefficient(j+1) = (C_coefficient(j+2) - C_coefficient(j+1)) / (3*h_stepSize(j+1));
end
 
% Output data
fprintf(1, '------ Cubic Spline Interpolation Data ------ \n');
fprintf(1, 'Nodes x(1), x(2), ... , x(n): \n');
for i = 1:numberOfNodes
    fprintf(1, 'x(%d) = %7.6f \n', i, x_values(i));
end
fprintf(1, '\n');
fprintf(1, 'Coefficients of the cubic spline, S, on each subinterval: \n\n'); 

fprintf(1, '           A(i)         B(i)         C(i)         D(i) \n');
for i = 1:numberOfNodes
    fprintf(1,'x(%d):    %7.6f     %7.6f     %7.6f     %7.6f \n', i, A_coefficient(i), B_coefficient(i), C_coefficient(i), D_coefficient(i));
end
 
% Graph of spline curve and nodes

% Plot the node points
plot(x_values, A_coefficient, 'r.');
% keep on the same graph
hold on; 

for i = 1:numberOfNodes-1
   % create vector of x values. 100 spaces even between each x.
   vector = linspace(x_values(i), x_values(i+1), 100);
   % Equation of interpolant. Plug in each vector value.
   S_interpolant = A_coefficient(i) + B_coefficient(i)*(vector-x_values(i)) + C_coefficient(i)*(vector-x_values(i)).^2 + D_coefficient(i)*(vector-x_values(i)).^3;
   % Plot the spline curve
   plot(vector, S_interpolant, 'b-');
end


title('Cubic Spline Plot')
xlabel('x')
ylabel('f(x)')
legend('Node points','Spline')    
grid on;    
 
 
 



