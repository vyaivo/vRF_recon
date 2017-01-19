function fwhm = rad2fwhm(rad)
% given size constant used to describe point at which basis function
% reaches zero from center of fcn, gives the FWHM of that basis function

pow = 7;

fwhm = 2*rad*acos((0.5^(1/pow)-0.5)/0.5)/pi;

end