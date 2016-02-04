function xprime = dxdt(t,x)

global Kvec Uini DX;

% G_Val = [];

% for i = 1:size(GroupIndices,1);
%     G_Val(i,1) = sum(x(GroupIndices{i,1}));
% end

xprime = DX(x,Kvec',1,Uini);% G_Val);

end