import re

def convert_citations(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Convert LaTeX citations to Sphinx citations
    converted_content = re.sub(r'~\\cite{(.*?)}', r' :cite:`\1`', content)
    
    # Convert \flecmd{...} to *...*
    converted_content = re.sub(r'\\flecmd{(.*?)}', r'*\1*', converted_content)

    # Convert \lmpcmd{...} to *...* with underscore replaced by space
    converted_content = re.sub(r'\\lmpcmd{(.*?)}', lambda m: '*' + m.group(1).replace('_', ' ') + '*', converted_content)
    
    # Convert \guicmd{...} to *...*
    converted_content = re.sub(r'\\guicmd{(.*?)}', r'*\1*', converted_content)

    # Convert \begin{lstlisting} to .. code-block:: lammps
    converted_content = re.sub(r'\\begin{lstlisting}', r'.. code-block:: lammps', converted_content)
    
    # Convert \end{lstlisting} to a blank line
    converted_content = re.sub(r'\\end{lstlisting}', r'\n', converted_content)

    # Convert \begin{note} to .. admonition:: Note
    converted_content = re.sub(r'\\begin{note}', r'.. admonition:: Note\n    :class: non-title-info', converted_content)
    
    # Convert \end{note} to a blank line
    converted_content = re.sub(r'\\end{note}', r'\n', converted_content)

    # Convert inline math expressions to Sphinx math notation
    converted_content = re.sub(r'\$(.*?)\$', r':math:`\1`', converted_content)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(converted_content)

input_tex_file = "../../../../../article/lammps-tutorials.tex"  # Change this to your actual input file
output_tex_file = "main.rst"  # Change this to your desired output file

convert_citations(input_tex_file, output_tex_file)
