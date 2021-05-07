% Calculate the significance (in a Bayesian framework)
% for the differences between WR -> SD -> PRN 
% in total integration during the resting state condition

% Whole cortex: WR -> SD
sig(1,1) = (sum(HI.Control{21, 1}.int_total.samples > HI.PreNap{21, 1}.int_total.samples))/1000
% 7 Networks: WR -> SD
sig(2,1) = (sum(HI.Networks.Control.N1{21, 1}.int_total.samples > HI.Networks.PreNap.N1{21, 1}.int_total.samples))/1000;
sig(3,1) = (sum(HI.Networks.Control.N2{21, 1}.int_total.samples > HI.Networks.PreNap.N2{21, 1}.int_total.samples))/1000;
sig(4,1) = (sum(HI.Networks.Control.N3{21, 1}.int_total.samples > HI.Networks.PreNap.N3{21, 1}.int_total.samples))/1000;
sig(5,1) = (sum(HI.Networks.Control.N4{21, 1}.int_total.samples > HI.Networks.PreNap.N4{21, 1}.int_total.samples))/1000;
sig(6,1) = (sum(HI.Networks.Control.N5{21, 1}.int_total.samples > HI.Networks.PreNap.N5{21, 1}.int_total.samples))/1000;
sig(7,1) = (sum(HI.Networks.Control.N6{21, 1}.int_total.samples > HI.Networks.PreNap.N6{21, 1}.int_total.samples))/1000;
sig(8,1) = (sum(HI.Networks.Control.N7{21, 1}.int_total.samples > HI.Networks.PreNap.N7{21, 1}.int_total.samples))/1000;

% Whole cortex: SD -> PRN
sig(1,2) = (sum(HI.PreNap{21, 1}.int_total.samples > HI.PostNap{21, 1}.int_total.samples))/1000;
% 7 Networks: SD -> PRN
sig(2,2) = (sum(HI.Networks.PreNap.N1{21, 1}.int_total.samples > HI.Networks.PostNap.N1{21, 1}.int_total.samples))/1000;
sig(3,2) = (sum(HI.Networks.PreNap.N2{21, 1}.int_total.samples > HI.Networks.PostNap.N2{21, 1}.int_total.samples))/1000;
sig(4,2) = (sum(HI.Networks.PreNap.N3{21, 1}.int_total.samples > HI.Networks.PostNap.N3{21, 1}.int_total.samples))/1000;
sig(5,2) = (sum(HI.Networks.PreNap.N4{21, 1}.int_total.samples > HI.Networks.PostNap.N4{21, 1}.int_total.samples))/1000;
sig(6,2) = (sum(HI.Networks.PreNap.N5{21, 1}.int_total.samples > HI.Networks.PostNap.N5{21, 1}.int_total.samples))/1000;
sig(7,2) = (sum(HI.Networks.PreNap.N6{21, 1}.int_total.samples > HI.Networks.PostNap.N6{21, 1}.int_total.samples))/1000;
sig(8,2) = (sum(HI.Networks.PreNap.N7{21, 1}.int_total.samples > HI.Networks.PostNap.N7{21, 1}.int_total.samples))/1000