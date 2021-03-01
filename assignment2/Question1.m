clear all
close all
clc
set(0, 'DefaultFigureWindowStyle', 'docked')

%Muhammad Shabeh Raza Abbas
%101092004

% QUESTION 1
% 2D Laplace case

Width = 2;
Len = 3;
V0 = 1;

dx = 0.2;
dy = 0.2;
nx = Len/dx;
ny = Width/dy;



C1 = -2*(1/dx^2 + 1/dy^2);
C2 = 1/(dx^2);
C3 = 1/(dy^2);
i = zeros(nx*ny,nx*ny);
j = zeros(nx*ny,1);


for x = 2:(nx-1)
    for y=2:(ny-1)
        c = (y-1).*nx + x;
        i(c,c) = C1;
        i(c,((y-1).*nx + x-1)) = C2;
        i(c,((y-1).*nx + x+1)) = C2;
        i(c,((y-2).*nx + x)) = C3;
        i(c,((y).*nx + x)) = C3;
    end
end


for y=1:ny
    c = ((y-1).*nx + 1);
    i(c,c) = 1;
    
    j(c) = V0;
    
    c = ((y-1).*nx + nx);
    i(c,c) = 1;
end


for x=2:(nx-1)
    c = (1-1).*nx + x;
    i(c,c) = 1;
    i(c,(2-1).*nx + x) = -1;
    
    c = ((ny-1).*nx + x);
    i(c,c) = 1;
    i(c,((ny-2).*nx + x)) = -1;
end

v = i\j;
v = reshape(v,[],ny)';

figure();
surf(linspace(0,Len,nx),linspace(0,Width,ny),v);
xlabel('x');
ylabel('y');
title(sprintf('2-D plot of V(x) = %.2f', dx));
set(gca, 'View', [45 45])

figure();
surf(linspace(0,Len,nx),linspace(0,Width,ny),v);
xlabel('x');
ylabel('y');
title(sprintf('2-D plot of V(x) = %.2f', dx));
set(gca, 'View', [45 45])
shading interp


Soln = zeros(ny, nx);
x1 = repmat(linspace(-Len/2,Len/2,nx),ny,1);
y1 = repmat(linspace(0,Width,ny),nx,1)';
it = 100;
avgError = zeros(it,1);


for c=1:it
    n = 2*c - 1;
    Soln = Soln + 1./n.*cosh(n.*pi.*x1./Width) ...
        ./cosh(n.*pi.*(Len./2)./Width).*sin(n.*pi.*y1./Width);

    avgError(c) = mean(mean(abs(Soln.*4.*V0./pi - v)));
end

Soln = Soln.*4.*V0./pi;

figure();
surf(linspace(0,Len,nx),linspace(0,Width,ny),Soln);
shading interp
xlabel('x');
ylabel('y');
title(sprintf('Analytical Graph with %d iterations', it));