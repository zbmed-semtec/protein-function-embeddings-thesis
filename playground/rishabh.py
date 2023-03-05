import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# create data
# x = [2, 3, 4]
# cbow_200 = [0.00348, 0.00362, 0.00359]
# cbow_400 = [0.00334, 0.00338, 0.00351]
# sg_200 = [0.000038, 0.000039, 0.000040]
# sg_400 = [0.000053, 0.000056, 0.000058]
# # plot lines
# plt.plot(x, cbow_200, label="cbow 200", linestyle="-")
# plt.plot(x, cbow_400, label="cbow 400", linestyle="-.")
# plt.plot(x, sg_200, label="sg 200", linestyle=":")
# plt.plot(x, sg_400, label="sg 400", linestyle="--")
# plt.ylabel("Precision")
# plt.xlabel("Min count")
# plt.locator_params(axis='x', nbins=3)
# plt.legend()
# plt.show()

df = pd.read_excel("./playground/HeightvsDisplacement.xlsx")
print(df)
