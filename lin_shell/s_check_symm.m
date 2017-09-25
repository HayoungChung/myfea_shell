u = textread('displacement.txt');
node = textread('node.txt');
bs = textread('bpts_sens_test.txt');

[~,sortx] = sortrows(node,[1,2,3]);
node1 = node(sortx,:);
u1 = u(sortx,:);

node_ry = node;
node_ry(:,2) = -node(:,2)+20;

[~,sortx] = sortrows(node_ry,[1,2,3]);
node2 = node(sortx,:);
u2 = u(sortx,:);
u2(:,1) = -u2(:,1);

% a y-axis symmtry is preserved
figure;
plot(node2(:,2),-u2(:,1),'o');
hold on
plot(node1(:,2),u1(:,1),'r+');

figure;
plot(node2(:,2),u2(:,2),'o');
hold on
plot(node1(:,2),u1(:,2),'r+');

figure;
plot(node2(:,2),u2(:,6),'o');
hold on
plot(node1(:,2),u1(:,6),'r+');

% figure(100); clf;
% subplot(1,3,1);
% scatter(node1(:,1),node1(:,2),40,u1(:,1),'filled')
% subplot(1,3,2);
% scatter(node1(:,1),node1(:,2),40,u1(:,2),'filled')
% subplot(1,3,3);
% scatter(node1(:,1),node1(:,2),40,u1(:,6),'filled')
% 
% figure(1); clf;
% scatter(node1(:,1),node1(:,2),40,sqrt(abs(sum(u1.^2,2))),'filled')
% figure(2); clf;
% scatter(node2(:,1),node2(:,2),40,sqrt(abs(sum(u2.^2,2))),'filled')

% absq = u1-u2;
% 
% figure(3); clf;
% scatter(node2(:,1),node2(:,2),40,sqrt(abs(sum(absq.^2,2))),'filled')