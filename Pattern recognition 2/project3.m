close all;
clear variables;
clc;

% means and covariances
mean_1 = [1;3;3];
mean_2 = [3;1;3];
mean_3 = [3;3;1];
cov_1 = [2.5 0.5 3.6;0.5 0.5 1.1; 3.6 1.1 8];
cov_2 = [2.5 2 1;2 2.5 0.5;1 0.5 2.5];
cov_3 = [2.5 1.5 1;1.5 2 1.5;1 1.5 2.5];

% samples
xa = mvnrnd(mean_1, cov_1, 100);
xb = mvnrnd(mean_2, cov_2, 100);
xc = mvnrnd(mean_3, cov_3, 100);


% transpose
x1=xa';
x2=xb';
x3=xc';

%show samples
plot3(x1(1,:),x1(2,:),x1(3,:),'ro')
grid on
hold on
plot3(x2(1,:),x2(2,:),x2(3,:),'go')
hold on
plot3(x3(1,:),x3(2,:),x3(3,:),'bo')
legend('x1','x2','x3')
title('x1 x2 x3')
xlabel('x1')
ylabel('x2')
zlabel('x3')

%mean of sample (pca)
x=[xa; xb; xc];
mx=(mean(x))';
figure,
plot3(mx(1,:),mx(2,:),mx(3,:),'yo', 'LineWidth', 5)
hold on
plot3(x1(1,:),x1(2,:),x1(3,:),'ro')
grid on
hold on
plot3(x2(1,:),x2(2,:),x2(3,:),'go')
hold on
plot3(x3(1,:),x3(2,:),x3(3,:),'bo')
title('mean')

%1
xx1=(x-repmat(mx',300,1))/300; %Xi-MeanXi
Sw=xx1'*xx1;

%2
%pca_1
figure,
% A = cov(x) ; 
[V,D] = eig(Sw);
[D, ind] = sort(diag(D), 'descend');
% V = V(:, ind);
% 
% w=zeros(3,2);
w(:,1)=V(:,ind(1,1));
w(:,2)=V(:,ind(2,1));

y1=w'*x1;
plot(y1(1,:),y1(2,:),'ro')
hold on
y2=w'*x2;
plot(y2(1,:),y2(2,:),'go')
hold on
y3=w'*x3;
plot(y3(1,:),y3(2,:),'bo')
grid on
legend('x1','x2','x3')
title('x1 x2 x3 pca 2d')
xlabel('x1')
ylabel('x2')
zlabel('x3')

%3
xz=(x-repmat(mx',300,1))/300; %Xi-MeanXi
Swa=xz'*xz;

%pca
figure,
% A = cov(x) ; 
[Va,Da] = eig(Swa);
[Da, inda] = sort(diag(Da), 'descend');
% V = V(:, ind);
% 
% w=zeros(3,2);
wz(:,1)=Va(:,inda(1,1));
wz(:,2)=Va(:,inda(2,1));
wz(:,3)=Va(:,inda(3,1));

yz1=wz'*x1;
plot3(yz1(1,:),yz1(2,:),yz1(3,:),'ro')
hold on
yz2=wz'*x2;
plot3(yz2(1,:),yz2(2,:),yz2(3,:),'go')
hold on
yz3=wz'*x3;
plot3(yz3(1,:),yz3(2,:),yz3(3,:),'bo')

grid on
legend('x1','x2','x3')
title('x1 x2 x3 pca 3d')
xlabel('x1')
ylabel('x2')
zlabel('x3')

figure,
plot3(yz1(1,:),yz1(2,:),yz1(3,:),'ro')
hold on
plot3(yz2(1,:),yz2(2,:),yz2(3,:),'go')
hold on
plot3(yz3(1,:),yz3(2,:),yz3(3,:),'bo')

grid on
legend('x1','x2','x3')
title('x1 x2 x3 pca 3d')
xlabel('x1')
ylabel('x2')
zlabel('x3')

Xa=xa-repmat(mean_1',100,1); %Xi-MeanXi
Sw1=Xa'*Xa;
[V1,D1] = eig(Sw1);
[D1, ind1] = sort(diag(D1), 'descend');
w1(:,1)=V1(:,ind1(1,1));
w1(:,2)=V1(:,ind1(2,1));
w1(:,3)=V1(:,ind1(3,1));
g1=(w1(:,1)/w1(:,3))'*x1;
plot3(g1(1,:),g1(2,:),g1(3,:),'rp')
hold on

Xb=xb-repmat(mean_2',100,1); %Xi-MeanXi
Sw2=Xb'*Xb;
[V2,D2] = eig(Sw2);
[D2, ind2] = sort(diag(D2), 'descend');
w2(:,1)=V2(:,ind2(1,1));
w2(:,2)=V2(:,ind(2,1));
w2(:,3)=V2(:,ind2(3,1));
g2=(w2(:,1)/w2(:,3))'*x2;
plot3(g2(1,:),g2(2,:),g2(3,:),'gp')
hold on

Xc=xc-repmat(mean_3',100,1); %Xi-MeanXi
Sw3=Xc'*Xc;
[V3,D3] = eig(Sw3);
[D3, ind3] = sort(diag(D3), 'descend');
w3(:,1)=V3(:,ind3(1,1));
w3(:,2)=V3(:,ind3(2,1));
w3(:,3)=V3(:,ind3(3,1));
g3=(w3(:,1)/w3(:,3))'*x3;
plot3(g3(1,:),g3(2,:),g3(3,:),'bp')
hold on

grid on
legend('x1','x2','x3')
title('x1 x2 x3 pca 3d')
xlabel('x1')
ylabel('x2')
zlabel('x3')





%mda
n1=100;
n2=100;
n3=100;
d=3;
c=3;
%% within-class scatter matrix
Xa=xa-repmat(mean_1',100,1); %Xi-MeanXi
Sw1=Xa'*Xa;
Xb=xb-repmat(mean_2',100,1); %Xi-MeanXi
Sw2=Xb'*Xb;
Xc=xc-repmat(mean_3',100,1); %Xi-MeanXi
Sw3=Xc'*Xc;
Sw=Sw1+Sw2+Sw3;
%% between -class scatter matrix
SB1=(n1)*(mean_1-mx)*(mean_1-mx)';
SB2=(n2)*(mean_2-mx)*(mean_2-mx)';
SB3=(n3)*(mean_3-mx)*(mean_3-mx)';
SB=SB1+SB2+SB3;
%% Computing the LDA projection W=Sw \ SB
W1=(inv(Sw))*SB;
% W=Sw \ SB;
%% getting the projection vectors
[V1,D1]=eig(W1);
[D1, ind1] = sort(diag(D1), 'descend');

w1=zeros(d,c-1);
w1(:,1)=V1(:,ind1(1,1));
w1(:,2)=V1(:,ind1(2,1));
figure,
ya=(w1'*x1);
plot(ya(1,:),ya(2,:),'rp')
hold on
yb=(w1'*x2);
plot(yb(1,:),yb(2,:),'gp')
hold on
yc=(w1'*x3);
plot(yc(1,:),yc(2,:),'bp')
grid on
title('mda')
legend('x1','x2','x3')
xlabel('x1')
ylabel('x2')
zlabel('x3')

%5
%% within-class scatter matrix
Xa1=xa-repmat(mean_1',100,1); %Xi-MeanXi
Sw1=Xa1'*Xa1;
Xb1=xb-repmat(mean_2',100,1); %Xi-MeanXi
Sw2=Xb1'*Xb1;

Swb=Sw1+Sw2;
%% Computing the LDA projection W=Sw \(m1-m2)
W_1=(inv(Swb))*(mean_1-mean_2);
% W=Sw \ SB;
%% getting the projection vectors
figure,
y_1=(W_1'*x1);
plot(y_1(1,:),0,'rp')
hold on
y_2=(W_1'*x2);
plot(y_2(1,:),0,'bp')
grid on
title('FLD')

%pdf
da=y_1';
m1=mean(da);
c1=da-repmat(m1,100,1); %Xi-MeanXi
cov1=c1'*c1;
figure
plot(da,normpdf(da,m1,cov1),'r*')
hold on
db=y_2';
m2=mean(db);
c2=db-repmat(m2,100,1); %Xi-MeanXi
cov2=c2'*c2;
plot(db,normpdf(db,m2,cov2),'b*')
grid on
title('PDF')