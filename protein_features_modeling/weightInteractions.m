function y = weightInteractions(N,effect,weight,medianIntxs)
y = max(weight.*(N-medianIntxs).*effect,0.99.*effect);
end