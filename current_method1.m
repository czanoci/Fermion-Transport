function [ current ] = current_method1(w, gamma, t)

current = 2*w*(gamma-t)/((t+gamma)*(4*w^2/gamma+(t+gamma)^2/t));
