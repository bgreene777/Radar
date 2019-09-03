function [meanx,stdx] = rvstats(x)
%
% rvstats.m
%
% [meanx,stdx] = rvstats(x)
%
% This program will calculate the mean and standard deviation
% of a random variable.
%
meanx=mean(x);
stdx=std(x);