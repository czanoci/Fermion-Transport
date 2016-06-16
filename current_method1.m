function [ current ] = current_method1(w, gamma, t)

current = 2*w*(1-t)/((t+1)*(4*w^2/gamma+gamma*(t+1)^2/t));
