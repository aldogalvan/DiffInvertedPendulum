import pandas as pd
import numpy as np
from scipy.optimize import minimize

def quadratic_function(x):
    
    return

file_path = 'trial.csv'

# Read the CSV file into a DataFrame
df = pd.read_csv(file_path)

# extract data
t = df.iloc[:, 0].to_numpy();
x = df.iloc[:,1:14].to_numpy();
xh = df.iloc[:,14:27].to_numpy();
fh = df.iloc[:,27:30].to_numpy();
th = df.iloc[:,30:33].to_numpy();

# perform the learning from demonstration


# test the new controller
dt = 0.01;
T = 10000;

# run simulation
s_k = np.zeros([T,13]);
u_k = np.zeros([T,13]);
l = 0.2; // the length
m = 1; // the mass
J; // the moment of inertia

for k in range(0,T):
    P = I * (M_Cart + M_Arm) + (M_Cart * M_Arm * length**2)

    A22 = (-(I + M_Arm * (length**2)) * b) / P
    A23 = ((M_Arm**2) * g * (length**2)) / P
    A42 = (-(M_Arm * length * b)) / P
    A43 = (M_Arm * g * length * (M_Cart + M_Arm)) / P
    
    B2 = (I + M_Arm * (length**2)) / P
    B4 = (M_Arm * length) / P
    
    
    A = np.array([[0, 1, 0, 0], [0, A22, A23, 0], [0, 0, 0, 1], [0, A42, A43, 0]])
    
    B = np.array([[0], [B2], [0], [B4]])
    
    C = np.array([[1, 0, 0, 0], [0, 0, 1, 0]])
    
    D = np.array([[0], [0]])

    
    