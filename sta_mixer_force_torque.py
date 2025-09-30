#!/usr/bin/env python3
"""
Post-processing script for mixer contact CSV files.
Calculates mesh forces and torques from particle contact data.

Can be used both as a command-line tool and as an importable module.
"""

import os
import pandas as pd
import numpy as np
import glob
from pathlib import Path
import argparse

def process_single_csv(filepath):
    """
    Process a single mixer_contact CSV file.
    
    Args:
        filepath: Path to the CSV file
    
    Returns:
        tuple: (frame_number, mesh_force, mesh_torque)
    """
    # Extract frame number from filename (xxxx part)
    filename = os.path.basename(filepath)
    frame_number = int(filename.replace('mixer_contact_', '').replace('.csv', ''))
    
    # Read CSV file
    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None
    
    # Filter for SM contact type
    sm_data = df[df['contact_type'] == 'SM']
    
    if sm_data.empty:
        # No SM contacts, return zero forces
        return frame_number, np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])
    
    # Sum forces and torques for all SM contacts
    total_force = sm_data[['f_x', 'f_y', 'f_z']].sum().values
    total_torque = sm_data[['torque_x', 'torque_y', 'torque_z']].sum().values
    
    # Mesh receives opposite forces (Newton's third law)
    mesh_force = -total_force
    mesh_torque = -total_torque
    
    return frame_number, mesh_force, mesh_torque

def process_directory(directory, output_file='mesh_forces_summary.csv', verbose=True):
    """
    Process all mixer_contact CSV files in a directory.
    
    Args:
        directory: Directory containing mixer_contact_xxxx.csv files
        output_file: Output CSV filename
        verbose: Print progress information
    
    Returns:
        DataFrame: Processed results
    """
    # Find all mixer_contact_*.csv files
    pattern = os.path.join(directory, 'mixer_contact_*.csv')
    csv_files = glob.glob(pattern)
    
    if not csv_files:
        if verbose:
            print(f"No mixer_contact_*.csv files found in {directory}")
        return None
    
    # Sort files by frame number
    csv_files.sort(key=lambda x: int(os.path.basename(x).replace('mixer_contact_', '').replace('.csv', '')))
    
    if verbose:
        print(f"Found {len(csv_files)} CSV files to process...")
    
    # Process each file and collect results
    results = []
    for filepath in csv_files:
        result = process_single_csv(filepath)
        if result:
            frame, force, torque = result
            results.append({
                'frame': frame,
                'mesh_force_x': force[0],
                'mesh_force_y': force[1],
                'mesh_force_z': force[2],
                'mesh_force_magnitude': np.linalg.norm(force),
                'mesh_torque_x': torque[0],
                'mesh_torque_y': torque[1],
                'mesh_torque_z': torque[2],
                'mesh_torque_magnitude': np.linalg.norm(torque)
            })
            
            # Print progress
            if verbose and frame % 100 == 0:
                print(f"  Processed frame {frame:04d}: Force magnitude = {np.linalg.norm(force):.6e} N, "
                      f"Torque magnitude = {np.linalg.norm(torque):.6e} N·m")
    
    # Create DataFrame and save to CSV
    if results:
        df_output = pd.DataFrame(results)
        
        # Save to file if output_file is specified
        if output_file:
            df_output.to_csv(output_file, index=False, float_format='%.8e')
            if verbose:
                print(f"\nResults saved to {output_file}")
        
        # Print summary statistics
        if verbose:
            print("\nSummary Statistics:")
            print("-" * 50)
            print(f"Total frames processed: {len(results)}")
            print(f"Average force magnitude: {df_output['mesh_force_magnitude'].mean():.6e} N")
            print(f"Max force magnitude: {df_output['mesh_force_magnitude'].max():.6e} N")
            print(f"Average torque magnitude: {df_output['mesh_torque_magnitude'].mean():.6e} N·m")
            print(f"Max torque magnitude: {df_output['mesh_torque_magnitude'].max():.6e} N·m")
        
        return df_output
    else:
        if verbose:
            print("No data to process")
        return None

def plot_results(df, output_dir='.', show_plot=True):
    """
    Create plots of forces and torques over time.
    
    Args:
        df: DataFrame with results
        output_dir: Directory to save plot
        show_plot: Whether to display the plot
    
    Returns:
        matplotlib figure object if successful, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        
        # Plot forces
        axes[0, 0].plot(df['frame'], df['mesh_force_x'], label='F_x', linewidth=1)
        axes[0, 0].plot(df['frame'], df['mesh_force_y'], label='F_y', linewidth=1)
        axes[0, 0].plot(df['frame'], df['mesh_force_z'], label='F_z', linewidth=1)
        axes[0, 0].set_xlabel('Frame')
        axes[0, 0].set_ylabel('Force [N]')
        axes[0, 0].set_title('Mesh Force Components')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot force magnitude
        axes[0, 1].plot(df['frame'], df['mesh_force_magnitude'], 'b-', linewidth=1)
        axes[0, 1].set_xlabel('Frame')
        axes[0, 1].set_ylabel('Force Magnitude [N]')
        axes[0, 1].set_title('Mesh Force Magnitude')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot torques
        axes[1, 0].plot(df['frame'], df['mesh_torque_x'], label='T_x', linewidth=1)
        axes[1, 0].plot(df['frame'], df['mesh_torque_y'], label='T_y', linewidth=1)
        axes[1, 0].plot(df['frame'], df['mesh_torque_z'], label='T_z', linewidth=1)
        axes[1, 0].set_xlabel('Frame')
        axes[1, 0].set_ylabel('Torque [N·m]')
        axes[1, 0].set_title('Mesh Torque Components')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot torque magnitude
        axes[1, 1].plot(df['frame'], df['mesh_torque_magnitude'], 'r-', linewidth=1)
        axes[1, 1].set_xlabel('Frame')
        axes[1, 1].set_ylabel('Torque Magnitude [N·m]')
        axes[1, 1].set_title('Mesh Torque Magnitude')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_dir:
            plot_file = os.path.join(output_dir, 'mesh_forces_plot.png')
            plt.savefig(plot_file, dpi=150)
            print(f"\nPlot saved to {plot_file}")
        
        if show_plot:
            plt.show()
        
        return fig
        
    except ImportError:
        print("\nNote: matplotlib not available, skipping plot generation")
        return None

def analyze_mixer_contacts(data_directory=None, 
                          output_file='mesh_forces_summary.csv',
                          generate_plot=False,
                          show_plot=True,
                          verbose=True):
    """
    Main analysis function that can be called directly from Python code.
    
    Args:
        data_directory: Directory containing mixer_contact_xxxx.csv files
                       If None, uses current directory
        output_file: Output CSV filename (None to skip saving)
        generate_plot: Whether to generate plots
        show_plot: Whether to display plots (if generate_plot is True)
        verbose: Print progress and statistics
    
    Returns:
        DataFrame with analysis results
    
    Example:
        # Basic usage
        df = analyze_mixer_contacts('path/to/data')
        
        # With plotting
        df = analyze_mixer_contacts('path/to/data', generate_plot=True)
        
        # Without saving output file
        df = analyze_mixer_contacts('path/to/data', output_file=None)
    """
    if data_directory is None:
        data_directory = '.'
    
    # Process the directory
    if verbose:
        print(f"Processing directory: {data_directory}")
    
    df = process_directory(data_directory, output_file, verbose)
    
    # Generate plots if requested
    if generate_plot and df is not None:
        plot_results(df, data_directory, show_plot)
    
    return df

def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(description='Post-process mixer contact CSV files')
    parser.add_argument('directory', nargs='?', default='.',
                        help='Directory containing mixer_contact_xxxx.csv files (default: current directory)')
    parser.add_argument('-o', '--output', default='mesh_forces_summary.csv',
                        help='Output CSV filename (default: mesh_forces_summary.csv)')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='Generate plots of forces and torques')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress verbose output')
    
    args = parser.parse_args()
    
    # Use the main analysis function
    df = analyze_mixer_contacts(
        data_directory=args.directory,
        output_file=args.output,
        generate_plot=args.plot,
        show_plot=True,
        verbose=not args.quiet
    )

# For direct code usage examples
if __name__ == "__main__":
    # Check if running interactively (e.g., in Jupyter or Python interpreter)
    import sys
    if len(sys.argv) == 1 and hasattr(sys, 'ps1'):
        print("Running in interactive mode. Examples:")
        print()
        print("# Example 1: Basic analysis with default settings")
        print("df = analyze_mixer_contacts('path/to/your/data')")
        print()
        print("# Example 2: Analysis with plotting")
        print("df = analyze_mixer_contacts('path/to/your/data', generate_plot=True)")
        print()
        print("# Example 3: Analysis without saving file")
        print("df = analyze_mixer_contacts('path/to/your/data', output_file=None)")
        print()
        print("# Example 4: Silent analysis")
        print("df = analyze_mixer_contacts('path/to/your/data', verbose=False)")
        print()
        print("# Example 5: Direct specification in code")
        print("# You can also directly specify the directory in the code:")
        print("# DATA_DIR = '/home/user/simulation/results'")
        print("# df = analyze_mixer_contacts(DATA_DIR, 'my_results.csv', generate_plot=True)")
    else:
        # Run as command-line tool
        main()