load ('initialElementData2.txt')
X=initialElementData2(:,1);
U=initialElementData2(:,2);
plot(X,U,'g');
hold on

load ('timesteppedElementData2.txt')
X=timesteppedElementData2(:,1);
U=timesteppedElementData2(:,2);
plot(X,U,'k');
hold on