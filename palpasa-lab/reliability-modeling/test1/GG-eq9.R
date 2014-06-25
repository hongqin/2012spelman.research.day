#nIn = n*q
nIn = c(50:5000)
t = 1:1000
k = 0.1
l = nIn

my.sum = 0; 
for ( i in 1:nIn) {
   numerator = l^(i-1) * (1 - exp(-k*t))^(i-1)
   denominator = factorial(i-1) * (1 - (1-exp(-k*t))^i)
  my.sum = my.sum + numerator / denominator 
} #i loop

  
  