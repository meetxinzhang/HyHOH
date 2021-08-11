import pandas as pd

df = pd.DataFrame({"Name": ["blue",
                            "delta",
                            "echo",
                            "charlie",
                            "alpha"],
                   "Type": ["Raptors",
                            "Raptors",
                            "Raptors",
                            "Raptors",
                            "Tyrannosaurus rex"]
                   })

print(df)
print('\n')

print(df.index[df['Name'].str.contains('ha')].tolist())
print('\n')

print(df.loc[df['Name'].str.contains('ha')])
print('\n')

print(df.loc[(df['Name'].str.contains('ha')) & (df['Type'].str.contains('Rex'))])
