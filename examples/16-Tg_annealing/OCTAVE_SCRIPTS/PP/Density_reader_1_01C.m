function [ datarray ] = Density_reader_1_01C(filename)
%Loads density/temperature data from the GROMACS custom Anneal script;
%In the format output GROMACS file, the first column is density data, whereas the second column is temperature data.


dum=csvread(filename,1,0);
datarray=[dum(:,1),dum(:,2)];

