%Tyler Yamori-Little
%Geodynamics Final Project
%4/26/21


%%
%import data
Antarctica_fingerprint = importdata('Antarctica_fingerprint.dat');
FP_lat = importdata('FP_lat.dat');
FP_lon = importdata('FP_lon.dat');
Greenland_fingerprint = importdata('Greenland_fingerprint.dat');
lats = importdata('lats.dat');
lons = importdata('lons.dat');
tidedata = importdata('tidedata.dat');
time = importdata('time.dat');

%%
%create a cell array that houses each tide site to a single row
TideSite = cell(100,1);                     %declare a cell array
for i = 1:1:100                             %interate through it
 TideSite{i} = [time(1,:);tidedata(i,:)];   %create a matrix that has time vs data
end


%%
%plot an individual tide guage to make sure it works
figure(1);
plot(TideSite{1,1}(1,:),TideSite{1,1}(2,:))
title('Sea Level Change at Site ', 1);
xlabel('time (s)');
ylabel('Sea Level Change (mm)'); 


%%
%editing the data in preparation of least squares fitting
EditedTideSite = cell(100,1);                       %declare a cell for future editing
RowTideSite = zeros(1,200);                         %declare a matrix for future editing
Row2TideSite = zeros(1,200);                        %declare a matrix for future editing
for j = 1:1:100                         
    RowTideSite = TideSite{j,1}(1,:);               %Set the edited row to the interated site (time values)
    Row2TideSite = TideSite{j,1}(2,:);              %Set the edited row to the interated site (water level values)
    isNaN = isnan(Row2TideSite);                    %find data points that are NaN
    RowTideSite(isNaN) = [];                        %get rid of those data points
    Row2TideSite(isNaN) = [];                       %get rid of those data points
 EditedTideSite{j,1} = [RowTideSite;Row2TideSite];  %Concatenate the changed data into a new cell array
end

%%
%Determine the rate at each site

TideSiteRate = cell(100,1);                                                                 %declare a cell for interation
TideSiteError = cell(100,1);                                                                %declare a cell for interation
weighted = zeros(100,1);
for l = 1:1:100                                                                 
    [TideSiteRate{l},TideSiteError{l}] = polyfit(EditedTideSite{l,1}(1,:),EditedTideSite{l,1}(2,:),1);       %rates determied by the Matlab linear squares regression function 
    if size(TideSiteError{l}.R) == 2                                                                            %also creates a structure to determine error
     C = (inv(TideSiteError{l}.R)*(inv(TideSiteError{l}.R))')*(TideSiteError{l}.normr)^2/(TideSiteError{l}.df); %calculate the covariance matrix
     sigma = sqrt(C(1,1));      %calculate sigma
     weighted(l,:) = 1/(sigma^2); %create a weighted matrix
    else
        weighted(l,:) = 0;      %some values for sigma are uncalculatable because the covariance matrix is 1x2
    end
    
end
W = eye(100);               %create an identity matrix
for i = 1:1:100
W(i,i) = weighted(i,:);     %the weighted matrix is equal to the weights in reduced row echelon form
end
%%
%Use the calculated rate and error to find a 95% level of condifence for
%the first site

[y_fit,delta] = polyval(TideSiteRate{1},EditedTideSite{1,1}(1,:),TideSiteError{1});                 %calculate delta off of the error function using polyval
figure(2);
plot(EditedTideSite{1,1}(1,:),EditedTideSite{1,1}(2,:),'b')                                        %plot the edited data
hold on
plot(EditedTideSite{1,1}(1,:),y_fit,'r-')                                                           %plot the linear regression
plot(EditedTideSite{1,1}(1,:),y_fit+2*delta,'m--',EditedTideSite{1,1}(1,:),y_fit-2*delta,'m--')     %plot the bounds of interval of confidence
title('Linear Fit of Data with 95% Prediction Interval')
legend('Data','Linear Fit','95% Prediction Interval')


%%
%create the A matrix by correlating lats and lons of the fingerprint to
%the tide sites
AntarcticaIndex = zeros(100,1); %declare a matrix for editing
GreenlandIndex = zeros(100,1);  %declare a matrix for editing
for i = 1:1:100
      indexlat = find(abs(FP_lat(:,1)-lats(i)) == min( abs(FP_lat(:,1)-lats(i)))); %correlate latitudes through finding the absolute difference between the site matrix and lats matrix
      indexlons = find(abs(FP_lon(:,1)-lons(i)) == min( abs(FP_lon(:,1)-lons(i))));%correlate longitudes through finding the absolute difference between the site matrix and lons matrix
      AntarcticaIndex(i,1) = Antarctica_fingerprint(indexlat,indexlons); %Build the fingerprint matrix for antarctica at the sites
     GreenlandIndex(i,1) = Greenland_fingerprint(indexlat,indexlons);    %Build the fingerprint matrix for greenland at the sites
end


%%
%Prepare to solve the matrix equation (Transpose(A)WA)^-1 * transpose(A)Wy
A = [AntarcticaIndex, GreenlandIndex]; %concatenate the A matrix using the fingerprint matricies

y = zeros(100,1);                      %decalre a matrix for the rates
for i=1:1:100
    y(i) = TideSiteRate{i}(1);         %pull out values from the cell and put into a matrix for y

end

%%
%solve with weight and without it


solution = inv(transpose(A)*A)*transpose(A)*y;
alpha = solution(1,:);
beta = solution(2,:);

solution_weighted = inv(transpose(A)*W*A)*transpose(A)*W*y;
alpha_weighted = solution_weighted(1,:);
beta_weigted  = solution_weighted(2,:);
