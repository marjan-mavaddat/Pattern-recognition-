function roc_curve1( mu_1, mu_2, sigma_1, sigma_2) % a function with four arguments

%samples
t_1 = mu_1-3*sigma_1:sigma_1/1000:mu_1+3*sigma_1;
t_2 = mu_2-3*sigma_2:sigma_2/1000:mu_2+3*sigma_2;

%Gaussian Generate
GaussianValues_1=(1./sqrt(2*pi*sigma_1)*exp(-0.5*((t_1-mu_1).^2)/sigma_1));
GaussianValues_2=(1./sqrt(2*pi*sigma_2)*exp(-0.5*((t_2-mu_2).^2)/sigma_2));

%show Gaussian
figure ,plot(t_1,GaussianValues_1);
hold on;
plot(t_2,GaussianValues_2,'r');
hold off
title('show PDF')
%calculate tresh
if mu_1-3*sigma_1 < mu_2-3*sigma_2
   x1= mu_1-3*sigma_1;
else
   x1=mu_2-3*sigma_2;
end   

if mu_1+3*sigma_1 < mu_2+3*sigma_2
   x2= mu_2+3*sigma_2;
else
   x2=mu_1+3*sigma_1;
end

%calculate and show ROC
figure,
i=1;
fa=zeros();
tp=zeros();
for tresh=x1:min(sigma_1/1000,sigma_2/1000):x2
    fa(i)=integral(@(x) normpdf(x, mu_1, sigma_1), tresh, x2);
    tp(i)=integral(@(x) normpdf(x, mu_2, sigma_2), tresh, x2);
    i=i+1;
end 
plot(fa,tp,'r*');
xlabel('false alarm') 
ylabel('TRUE detection')
title('ROC Curve')
% error 'must be mu_1 < mu_2'
if mu_2 < mu_1
error('Error: must be mu_1 < mu_2'); 
end
end

