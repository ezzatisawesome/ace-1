import pandas as pd

# Read the CSV file
drag_table_df = pd.read_csv('Drag Table.csv')

# Convert to dictionary where each row is a dictionary with Component as the key
drag_table_dict = drag_table_df.set_index('Component').to_dict('index')

# Print the dictionary to verify
print("Drag Table Dictionary:")
for component, values in drag_table_dict.items():
    print(f"\n{component}:")
    for key, value in values.items():
        print(f"  {key}: {value}")

# Example of how to access values:
# drag_table_dict['Fuselage']['L']  # Gets the length of the fuselage
# drag_table_dict['Main Wing']['S_Wet']  # Gets the wetted area of the main wing 