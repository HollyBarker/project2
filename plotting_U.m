load ('timesteppedElementData.txt')
X=timesteppedElementData(:,1);
U=timesteppedElementData(:,2);
plot(X,U,'g');
hold on