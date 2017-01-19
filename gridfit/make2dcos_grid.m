function ff = make2dcos_grid(eval_at,pp)
% pp(1) is mux, pp(2) is muy, pp(3) is size, pow is 7

myr = ((eval_at(:,1)-pp(1)).^2+(eval_at(:,2)-pp(2)).^2).^0.5;

ff = ( ((0.5*(1 + cos(myr*pi/pp(3)) )).^7) .* (myr<=pp(3)) );


return