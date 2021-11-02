function [targets] = find_target_phase(EEG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% EEG: 1 x time points
% target has 5 type of elememnt 1=>0 degree 2=.90 degree 3=>180 degree
% 4=>270 degree 0=> otherwise
targets=zeros(size(EEG));
[~,loc2]=findpeaks(EEG);
disp('loc2');
disp(loc2);
[~,loc4]=findpeaks(-EEG);
disp('loc4');
disp(loc4);
grad=gradient(EEG);
[~,loc1]=findpeaks(grad);
[~,loc3]=findpeaks(-grad);
for i=loc1
    targets(i)=1;
end
for i=loc2
    targets(i)=2;
end
for i=loc3
    targets(i)=3;
end
for i=loc4
    targets(i)=4;
end
end

