import rstparse, os, shutil
import numpy as np
from PIL import Image

from utilities import count_line, replace_special_character, fix_math, fix_link, fix_italic

class WriteTex:
    """Write Tex file."""
    def __init__(self, file_name, RST, git_path, *args, **kwargs,):
        """Initialize"""
        super().__init__(*args, **kwargs)
        self.file_name = file_name
        self.RST = RST
        self.git_path = git_path

    def convert_file(self):
        """Main convert function."""
        self.f = open(self.file_name, "w") 
        self.write_main_title()

        for block_id in np.unique(self.RST.main_block):
            raw_block, block_type, block_lines = clean_block(self.RST.file_content,
                                                             self.RST.main_block_type,
                                                             self.RST.main_block, block_id)
            filtered_block, n_subblock = filter_block(raw_block, block_type)

            # Look for title
            if block_lines[0] in self.RST.title_positions:
                title = self.RST.file_content[block_lines[0]]
                filtered_block = title
                block_type = np.array(self.RST.title_types)[np.array(self.RST.title_positions) == block_lines[0]]

            ids_subblock, types_subblock = identify_subblock_id(n_subblock, block_lines, self.RST)
            filtered_subblock, sub_block_number = read_sublock(n_subblock, filtered_block, ids_subblock, self.RST)

            self.write_paragraph(filtered_block, block_type, filtered_subblock, ids_subblock, types_subblock, sub_block_number)
            self.write_equation(filtered_block, block_type)
            self.write_title(filtered_block, block_type)
        self.f.close()

    def write_main_title(self):
        title_position = self.RST.title_positions[np.where(np.array(self.RST.title_types) == "main")[0][0]]
        self.f.write('\chapter{'+self.RST.file_content[title_position]+'}')
        self.f.write('\n')

    def write_paragraph(self, filtered_block, block_type, filtered_subblock, ids_subblock, types_subblock, sub_block_number):
        if ("text" in block_type):
            for line in filtered_block:
                line = replace_special_character(line, '#', r'$\#$')
                line = replace_special_character(line, '*->*', r'$\rightarrow$')
                line = fix_link(self.RST, line)
                line = fix_math(line)
                line = fix_italic(line, replace_underscore=True)
                self.f.write(line)
                self.f.write('\n')
        elif ("admonition" in block_type):
            cpt = 0
            caption = block_type[11:]
            self.f.write(r'\begin{tcolorbox}[colback=mylightblue!5!white,colframe=mylightblue!75!black,title='+caption+']'+'\n')
            #self.f.write(r'\noindent \textbf{' + caption + '} -- ')
            for line in filtered_block:
                line = fix_math(line)
                line = fix_link(self.RST, line)
                line = fix_italic(line, replace_underscore=True)
                if line == '[insert-sub-block]': 
                    filtered_subblock_0 = []
                    for line, n in zip(filtered_subblock, sub_block_number):
                        if n == cpt:
                            filtered_subblock_0.append(line)
                    if 'figure' not in types_subblock[cpt]:
                        self.write_equation(filtered_subblock_0, types_subblock[cpt])
                    cpt += 1
                else:
                    self.f.write(line)
                    self.f.write('\n')
            self.f.write(r'\end{tcolorbox}')
        elif "hatnote" in block_type:
            for line in filtered_block:
                line = fix_math(line)
                self.f.write(r'\vspace{-1cm} '
                + r'\noindent \textcolor{graytitle}{\textit{{\Large '+line+
                '}}'
                + r'\vspace{0.5cm} }')
                self.f.write('\n')
        elif "figure" in block_type:
            if "dark" in block_type:
                pass
            else:
                align = None
                for line in filtered_block:
                    if "height" in line:
                        height = line[9:]
                    elif "align" in line:
                        align = line[8:]
                path = block_type[9:]
                figure_path = self.git_path+'/docs/sphinx/source/tutorials/'+path[3:]
                if os.path.exists(figure_path) is False:
                    print("Figure not found", figure_path)
                figure_format = figure_path.split('.')[-1]
                level = figure_path.split('/')[-3]
                tutorial = figure_path.split('/')[-2]
                name = figure_path.split('/')[-1].split('.')[0]
                if os.path.exists(self.git_path+'/ebook/tutorials/'+level+'/'+tutorial) is False:
                    os.mkdir(self.git_path+'/ebook/tutorials/'+level+'/'+tutorial+'/')
                alternative_figure = figure_path[:-len(figure_format)]+'png'
                new_figure = self.git_path+'/ebook/tutorials/'+level+'/'+tutorial+'/'+name+'.png'
                if os.path.exists(alternative_figure):
                    shutil.copyfile(alternative_figure, new_figure)
                else:      
                    print("webp convert into png")  
                    print(new_figure)            
                    im = Image.open(figure_path).convert('RGB')
                    im.save(new_figure, 'png')

                if align is None:
                    self.f.write(r'\begin{figure}'+'\n')
                    self.f.write(r'\includegraphics[width=\linewidth]{tutorials/'+level+'/'+tutorial+'/'+name+'.png}'+'\n')
                    self.f.write(r'\end{figure}'+'\n')  
                elif 'right' in align:
                    self.f.write(r'\hspace{-0.45cm}')
                    self.f.write(r'\begin{wrapfigure}{r}{4cm}'+'\n')
                    self.f.write(r'\includegraphics[width=4cm]{tutorials/'+level+'/'+tutorial+'/'+name+'.png}'+'\n')
                    self.f.write(r'\end{wrapfigure}'+'\n')
        self.f.write('\n')

    def write_title(self, filtered_block, block_type):
        if ("subtitle" in block_type):
            self.f.write('\section{' + filtered_block + '}')
            self.f.write('\n')
        elif ("subsubtitle" in block_type):
            self.f.write('\subsection{' + filtered_block + '}')
            self.f.write('\n')

    def write_equation(self, filtered_block, block_type):
        if ("lammps" in block_type):
            self.f.write(r'\begin{lcverbatim}'+'\n')
            for line in filtered_block:
                if ':caption:' not in line:
                    self.f.write(line)
                    self.f.write('\n')
            self.f.write(r'\end{lcverbatim}'+'\n')
        elif ("bw" in block_type) | ("bash" in block_type) | ("python" in block_type):
            self.f.write(r'\begin{lcverbatim}'+'\n')
            for line in filtered_block:
                if ':caption:' not in line:
                    self.f.write(line)
                    self.f.write('\n')
            self.f.write(r'\end{lcverbatim}'+'\n')
        elif ("math" in block_type):
            for line in filtered_block:
                    if len(line) > 0:
                        self.f.write('$$' + line + '$$')
        self.f.write('\n')

def clean_block(content, types, ids, id):
    line_numbers = np.where(np.array(ids) == id)[0]
    if len(line_numbers) > 0:
        block_lines = []
        block_types = []
        for n in line_numbers:
            line = content[n]
            type = types[n]
            block_lines.append(line)
            block_types.append(type)
        block_type = np.unique(block_types)
        assert len(block_type) == 1, """Block of mixed type ?"""
        return block_lines, block_type[0], line_numbers
    else:
        return None, None, None
    
def filter_block(block_text, block_type):
    """Convert a block of text into clean lines"""
    first_non_zero_indentation = None
    itemize_indentation = -1
    indentation = -1
    n_subblock = 0
    if (block_text is not None) & (block_type != 'unknown'):
        # put the lines into a list 
        filtered_text = []
        for line in block_text[1:]:
            indentation = count_line(line)
            if (first_non_zero_indentation is None) & (indentation > 0) & (len(line) > 0):
                first_non_zero_indentation = indentation

            # detect itemized
            if first_non_zero_indentation is not None:
                if len(line)> first_non_zero_indentation:
                    if line[first_non_zero_indentation] == '-':
                        itemize_indentation = indentation

            if (':class:' not in line) & (len(line) > 0) & (first_non_zero_indentation != None):
                if (indentation == first_non_zero_indentation) | (indentation == itemize_indentation+2):
                    if ('..' in line) & ('...' not in line) & ('../' not in line):
                        filtered_text.append('[insert-sub-block]')
                        n_subblock += 1
                    else:
                        filtered_text.append(line[first_non_zero_indentation:])
        return filtered_text, n_subblock
    else:
        return None, 0
    
def identify_subblock_id(n_subblock, block_lines, RST):
    if n_subblock > 0:
        
        sub_block = []
        for n in block_lines:
            type = RST.sub_block_type[n]
            id = RST.sub_block[n]
            if type != 'unknown':
                sub_block.append(id)
        id = np.unique(sub_block)
        assert len(id) == n_subblock, print(block_lines, n_subblock, id)

        sub_block_type = []
        for id0 in id:
            type0 = np.unique(np.array(RST.sub_block_type)[np.array(RST.sub_block) == id0]).tolist()[0]
            sub_block_type.append(type0)

        return id, sub_block_type
    else:
        return None, None 
    
def read_sublock(n_subblock, filtered_block, ids_block, RST):
    if n_subblock > 0:
        cpt = 0
        sub_block = []
        sub_block_number = []
        for line in filtered_block:
            if line == '[insert-sub-block]':
                id = ids_block[cpt]
                raw_subblock, subblock_type, subblock_lines = clean_block(RST.file_content,
                                                            RST.sub_block_type,
                                                            RST.sub_block, id)   
                filtered_subblock, n_subsubblock = filter_block(raw_subblock, subblock_type)
                assert n_subsubblock == 0, """WARNING, subsub blocks were used"""
                
                for subline in filtered_subblock:
                    sub_block_number.append(cpt)
                    sub_block.append(subline)
                cpt += 1  
        return sub_block, sub_block_number
    else:
        return None, None