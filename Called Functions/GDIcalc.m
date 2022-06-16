function [GDIs] = GDIcalc(SubjKinematics)


[M,N] = size(SubjKinematics);

if [M,N] == [51 24] | [M,N] == [101 24]
    SubjKine = Column2Matrix(SubjKinematics);
else [M,N] == [918 1];
    SubjKine = SubjKinematics;
end

filename = 'GDI-GPS calculator v 3.2.xlsx';
xlswrite(filename,SubjKine,'Subject Kinematics','E2:E919');

GDIs = xlsread(filename,'Subject Analysis','E2:E3');



end