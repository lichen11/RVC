% function [ TwoClsA, TwoClsLab ] = makeTwoCls( originalA, originalLab, deleteLab )
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% %deleteLab = originalLab(1:length(originalLab) ~= WantedTwoCls);
% deletedId = find(originalLab ~= deleteLab);
% TwoClsA = [originalLab' ; originalA]; %this is the row bind step
% [d1, ~] = size(TwoClsA);
% %[~, d2] = sort(combinedA(1,:));
% TwoClsA = TwoClsA(:, 1:d1 ~= deletedId);
% %combinedA = combinedA(:,d2);
% TwoClsA = TwoClsA(2:d1,:);
% TwoClsA = [originalLab TwoClsA];
% %combinedA = sortrows(combinedA, 1);
% TwoClsA = TwoClsA(1:d1 ~= deletedId, :);
% TwoClsA = TwoClsA(:, 2:d1);
% TwoClsLab = originalLab(1:length(originalLab) ~= deletedId);
% 
% 
% 
% end
% 
