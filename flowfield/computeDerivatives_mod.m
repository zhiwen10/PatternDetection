function [fx, fy, ft] = computeDerivatives(im1, im2)

if size(im2,1)==0
    im2=zeros(size(im1));
end
im1 = angle(exp(i*(im1+pi)));
im2 = angle(exp(i*(im2+pi)));
% Horn-Schunck original method
fx = angle(exp(i*conv2(im1,0.25* [-1 1; -1 1],'same'))) + angle(exp(i*conv2(im2, 0.25*[-1 1; -1 1],'same')));
fy = angle(exp(i*conv2(im1, 0.25*[-1 -1; 1 1], 'same'))) + angle(exp(i*conv2(im2, 0.25*[-1 -1; 1 1], 'same')));
ft = angle(exp(i*conv2(im1, 0.25*ones(2),'same'))) + angle(exp(i*conv2(im2, -0.25*ones(2),'same')));

fx = angle(exp(i*fx));
fy = angle(exp(i*fy));
ft = angle(exp(i*ft));
% derivatives as in Barron
% fx= conv2(im1,(1/12)*[-1 8 0 -8 1],'same');
% fy= conv2(im1,(1/12)*[-1 8 0 -8 1]','same');
% ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
% fx=-fx;fy=-fy;

% An alternative way to compute the spatiotemporal derivatives is to use simple finite difference masks.
% fx = conv2(im1,[1 -1]);
% fy = conv2(im1,[1; -1]);
% ft= im2-im1;