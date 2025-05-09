p=10
n11=200
n22=200
M=10
n1=n11*M
n2=n22*M
data=generate_model_1(p,n1,n2)
X=data$data_I
Y=data$data_II
DM=data$Omega_II-data$Omega_I

D_est=DDME(X_data = X, Y_data=Y, M)