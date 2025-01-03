% Function to compute fixed points of a 2d system
function [U,eval,v1,v2] = compute_fp(f,IC1,IC2,delta1,delta2)
X1=IC1(1):delta1:IC1(2);
X2=IC2(1):delta2:IC2(2);
[E, I] = meshgrid(X1,X2);
gridpts=[E(:), I(:)];
FixedPTS=NaN(length(gridpts),2);
FixedPTS2=NaN(length(gridpts),2);
Exit=NaN(length(gridpts),1);
options = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);
for i = 1:length(gridpts)
    [FixedPTS(i,:), fval, Exit(i)] = fsolve(@(x) f(x), gridpts(i,:), options);
    if Exit(i)==1
        FixedPTS2(i,:)=FixedPTS(i,:);
    end
end
out = NaN(1,1);
FixedPTS3=FixedPTS2(~any(isnan(FixedPTS2), 2),:);
UniqueEig = unique(round(FixedPTS3,6), 'rows');
for i=1:size(UniqueEig,1)
    U(i,1)= UniqueEig(i,1);
    U(i,2)= UniqueEig(i,2);
    A = jacobianest(@(x) f(x), UniqueEig(i,:));
    [v, evalue] = eig(A);
    eval(i,:) = [evalue(1,1);evalue(2,2)];
    v1(i,:) = [v(:,1)];
    v2(i,:) = [v(:,2)];
end
end
