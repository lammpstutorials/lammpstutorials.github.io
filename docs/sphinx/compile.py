import subprocess
import re

# Define ANSI escape codes for colors
RED = '\033[91m'
GREEN = '\033[92m'
GRAY = '\033[90m'
RESET = '\033[0m'

def filter_make_output():

    # Define multiple patterns to ignore
    ignore_patterns = [
        "Pygments lexer name 'bw' is not known",
        "Pygments lexer name 'lammps' is not known",
        ".. label:: start_",
        ".. label:: end_",
        "Unknown directive type"
    ]
    
    # Define a pattern to identify warnings (example pattern, adjust as needed)
    warning_pattern = re.compile('|'.join(re.escape(p) for p in ["WARNING:", "ERROR:"]))
    
    # Combine ignore patterns into a single regex pattern
    ignore_pattern = re.compile('|'.join(re.escape(p) for p in ignore_patterns))
    
    # Run 'make clean'
    subprocess.run(['make', 'clean'], check=True)
    
    # Run 'make html' and capture output
    process = subprocess.Popen(['make', 'html'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True)
    
    # Read and filter output
    output_lines = []
    for line in process.stdout:
        if len(line) > 1:
            if not ignore_pattern.search(line):
                # Determine the color based on whether the line matches the warning pattern
                if warning_pattern.search(line):
                    output_lines.append(RED + line + RESET)
                else:
                    output_lines.append(GRAY + line + RESET)

    # Wait for the process to complete
    process.wait()

    # Print the filtered output
    print(''.join(output_lines), end='')

if __name__ == "__main__":
    filter_make_output()
