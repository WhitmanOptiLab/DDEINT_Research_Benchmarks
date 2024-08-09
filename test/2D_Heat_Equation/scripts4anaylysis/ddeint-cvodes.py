import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_comparison(array1, array2, time_array):
  """Plots a comparison of two arrays against time.

  Args:
    array1: The first array of values.
    array2: The second array of values.
    time_array: The array of time values.

  Returns:
    None
  """

  # Check if arrays have the same length
  if len(array1) != len(array2) or len(array1) != len(time_array):
    raise ValueError("Arrays must have the same length")

  # Create the plot
  plt.figure(figsize=(10, 6))
  plt.plot(time_array, array1, label="DDEint", linestyle='--', color='red', linewidth=2)
  plt.plot(time_array, array2, label="CVODES", linestyle='-', color='blue', linewidth=1)     
  plt.xlabel("t")
  plt.ylabel("u")
  plt.title("DDEINT-CVODES Comparison")
  plt.legend()
  plt.grid(True)
  plt.show()

# Example usage:
ddeint = pd.read_csv("data/heat_equation_output.csv")
DDEint = ddeint["||u||_rms"].values  
time_array = ddeint["time"].values

cvodes =   np.array([   
  3.551986624101131e-01,    
  3.296219563422419e-01 ,   
 2.894444105452381e-01   , 
  2.387095096822459e-01  ,  
 1.824250183424770e-01  ,  
 1.261154243551791e-01  ,  
  7.529884566729866e-02 ,   
 3.495180707843489e-02  ,  
9.024097530846483e-03   , 
 5.539756664495849e-05  ,  
  8.918960792100292e-03 ,   
 3.475300957548041e-02 ,   
 7.502652959091248e-02 ,   
1.257971846174588e-01  ,  
 1.820956996622664e-01 ,   
 2.384117838797803e-01 ,   
 2.892283942350709e-01 ,   
 3.295758083798311e-01 ,   
 3.555036931175249e-01 ,   
 3.644737367140586e-01])


plot_comparison(DDEint, cvodes, time_array)


 # Create the plot



