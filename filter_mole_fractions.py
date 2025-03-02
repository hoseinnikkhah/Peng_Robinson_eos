import pandas as pd
import os
import sys

def process_mole_fractions(input_file='mole_fraction_comparison.csv', output_file='data_point.csv'):
    try:
        # Read the input CSV file
        print(f"Reading data from {input_file}...")
        df_input = pd.read_csv(input_file)
        
        # Filter rows where RelativeDifferencePercent is less than or equal to 5
        df_filtered = df_input[df_input['RelativeDifferencePercent'] <= 5]
        
        if df_filtered.empty:
            print("No data points with RelativeDifferencePercent <= 5 found.")
            return
        else:
            print(f"Found {len(df_filtered)} data points with RelativeDifferencePercent <= 5.")
        
        # Check if output file already exists
        if os.path.exists(output_file):
            # Read existing output file to avoid duplicates
            df_existing = pd.read_csv(output_file)
            
            # Use InitialMoleFraction as a key to identify already processed data
            existing_mole_fractions = set(df_existing['InitialMoleFraction'].values)
            
            # Filter out data that's already in the output file
            df_new = df_filtered[~df_filtered['InitialMoleFraction'].isin(existing_mole_fractions)]
            
            if df_new.empty:
                print("No new data points to add. All filtered data is already in the output file.")
                return
            
            print(f"Appending {len(df_new)} new data points to {output_file}...")
            
            # Append new data to the output file without header
            df_new.to_csv(output_file, mode='a', header=False, index=False)
            
        else:
            # Create a new output file with filtered data
            print(f"Creating new file {output_file} with {len(df_filtered)} data points...")
            df_filtered.to_csv(output_file, index=False)
        
        print("Processing completed successfully!")
        
    except Exception as e:
        print(f"Error occurred: {e}")
        return
    
if __name__ == "__main__":
    # Get input and output file names from command line arguments if provided
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'mole_fraction_comparison.csv'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'data_point.csv'
    
    # Process the data
    process_mole_fractions(input_file, output_file)
    
    # Exit the script
    sys.exit(0)
