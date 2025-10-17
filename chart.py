import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Configure matplotlib for better display of Chinese characters (if needed)
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']  # or ['SimHei'] for Chinese
plt.rcParams['axes.unicode_minus'] = False

def plot_csv_data(csv_path):
    """
    Read CSV file and plot time series data
    
    Parameters:
    csv_path (str): Path to the CSV file
    """
    try:
        # Read CSV file
        df = pd.read_csv(csv_path)
        
        # Print column names for verification
        print(f"Columns in CSV: {df.columns.tolist()}")
        print(f"Data shape: {df.shape}")
        print(f"\nFirst few rows:")
        print(df.head())
        
        # Get column names (assuming first column is Time)
        time_col = df.columns[0]
        data_cols = df.columns[1:]
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        
        # Plot 1: Both signals on same axis
        ax1 = axes[0]
        for col in data_cols:
            ax1.plot(df[time_col], df[col], label=col, linewidth=2, marker='o', markersize=4)
        ax1.set_ylabel('Amplitude')
        ax1.set_title('Time Series Data - Combined View')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='best')
        
        # Plot 2: Separate traces
        ax2 = axes[1]
        colors = plt.cm.Set1(np.linspace(0, 1, len(data_cols)))
        for i, col in enumerate(data_cols):
            # Normalize for better comparison if needed
            data = df[col]
            offset = i * 0.002  # Add small offset for visibility
            ax2.plot(df[time_col], data + offset, label=f'{col} (offset: {offset:.3f})', 
                    color=colors[i], linewidth=2, alpha=0.8)
        
        ax2.set_xlabel(f'{time_col}')
        ax2.set_ylabel('Amplitude (with offset)')
        ax2.set_title('Time Series Data - Separated View')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='best')
        
        plt.tight_layout()
        
        # Additional statistics
        print(f"\n=== Data Statistics ===")
        for col in data_cols:
            print(f"\n{col}:")
            print(f"  Min: {df[col].min():.6f}")
            print(f"  Max: {df[col].max():.6f}")
            print(f"  Mean: {df[col].mean():.6f}")
            print(f"  Std: {df[col].std():.6f}")
        
        # Show the plot
        plt.show()
        
        # Optional: Save the figure
        save_option = input("\nSave figure? (y/n): ").strip().lower()
        if save_option == 'y':
            output_path = Path(csv_path).stem + "_plot.png"
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved as: {output_path}")
        
    except FileNotFoundError:
        print(f"Error: File '{csv_path}' not found!")
    except Exception as e:
        print(f"Error processing file: {e}")

def main():
    """
    Main function to run the CSV plotter
    """
    # Option 1: Hardcoded path (modify this to your CSV file path)
    csv_file_path = "mixer_torque_comparison.csv"
    
    # Option 2: Interactive input (uncomment to use)
    # print("CSV Time Series Data Plotter")
    # print("-" * 30)
    # csv_file_path = input("Enter CSV file path: ").strip()
    
    # Check if file exists
    if not Path(csv_file_path).exists():
        print(f"Error: File '{csv_file_path}' not found!")
        print("Please ensure the CSV file exists at the specified path.")
        return
    
    # Plot the data
    plot_csv_data(csv_file_path)

if __name__ == "__main__":
    main()