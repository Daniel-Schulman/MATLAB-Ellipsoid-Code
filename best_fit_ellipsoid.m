clear;
clc;

%Reads the point cloud data from inserted PCD file and stores it in variable A.
A = pcread('');

%Computes the centroid of the point cloud, Centers the point cloud around its centroid, storing the result in B.
centroid = mean(A.Location,1);
B = pointCloud(A.Location-mean(A.Location,1));

% figure;
% pcshow(B);

% Fits an ellipsoid to the centered point cloud B and returns parameters such as center, radii, eigenvectors, etc.
[ center, radii, evecs, v, chi2 ] = ellipsoid_fit_new(B.Location);

%Functions to calculate Cartesian coordinates from spherical coordinates.
x = @(r,theta,phi)r.*sind(theta).*cosd(phi);
y = @(r,theta,phi)r.*sind(theta).*sind(phi);
z = @(r,theta,phi)r.*cosd(theta);

%Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
F = @(r,theta,phi)[x(r,theta,phi).^2, y(r,theta,phi).^2, z(r,theta,phi).^2,...
    2*x(r,theta,phi).*y(r,theta,phi), 2*x(r,theta,phi).*z(r,theta,phi), 2*y(r,theta,phi).*z(r,theta,phi),...
    2*x(r,theta,phi), 2*y(r,theta,phi), 2*z(r,theta,phi), 1]*v;

%Range of theta (pi/2 to -pi/2) and phi (0 to 2pi) 100 points in those
%ranges
theta = linspace(90, -90,100);
phi = linspace(0, 360,100);

%Creates a grid of combinations of theta and phi values.
[Theta, Phi] = meshgrid(theta,phi);

r = zeros(size(Theta));

%Uses a nested loop to iterate over each combination of theta and phi.
%For each combination, fzero is used to find the root of the function F to
%get the corresponding radius r.
%Populates the matrix r with these computed radii.
for i = 1:length(theta)
    for j = 1:length(phi)
        F2 = @(r)F(r,theta(i),phi(j));
        r(j,i) = fzero(F2, 1);      
    end
end

%scatter3 is used to plot points in 3D space corresponding to the calculated radii.
%can remove hold on and pcshow(B) on this to just show the ellipsoid on its
%own
figure 
scatter3(x(r,Theta,Phi),y(r,Theta,Phi),z(r,Theta,Phi), "cyan")
hold on
pcshow(B);

%for creating radii of the semi-major axises
%when running this command, radii= are the values of the semi-major axises
radii
m1=evecs(:,1);
m2=evecs(:,2);
m3=evecs(:,3);
X=[0, m1(1)*radii(1), m2(1)*radii(2), m3(1)*radii(3)];
Y=[0, m1(2)*radii(1), m2(2)*radii(2), m3(2)*radii(3)];
Z=[0, m1(3)*radii(1), m2(3)*radii(2), m3(3)*radii(3)];
%WIP;
hold on
plot3(X([1 2]),Y([1 2]),Z([1 2]),...
      X([1 3]),Y([1 3]),Z([1 3]), ...
      X([1 4]),Y([1 4]),Z([1 4]))
% quiver3(xyz,xyz,xyz,evecs(:,1)'*radii(1),evecs(:,2)'*radii(2),evecs(:,3)'*radii(3));


