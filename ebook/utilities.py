import rstparse, os, shutil
import numpy as np
from PIL import Image

class FixDocument:
    """Fix Tex file."""
    def __init__(self, file_name, RST, git_path, *args, **kwargs,):
        """Initialize"""
        super().__init__(*args, **kwargs)
        self.file_name = file_name
        self.RST = RST
        self.git_path = git_path

    def fix_document(self):
        """Improve the final document"""
        self.remove_space()
        keywords = [r'\end{lcverbatim}', r'\Large', r'\section{',
                    r'\subsection{', r'\end{wrapfigure}', r'\end{tcolorbox}']
        self.add_non_indent(keywords)
        self.convert_itemize()

    def add_non_indent(self, keywords):
        initial_tex_file = self.import_tex_file()
        new_tex_file_name = []
        add_noindent = False
        for line in initial_tex_file:
            if line != '\n':
                prev_add_noindent = add_noindent
                add_noindent = False
                for word in keywords:
                    if word in line:
                        add_noindent = True
            if prev_add_noindent: 
                new_tex_file_name.append(r'\noindent '+line)
                add_noindent = False
                prev_add_noindent = False
            else:
                new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)

    def remove_space(self):
        initial_tex_file = self.import_tex_file()
        # remove extra space
        new_tex_file_name = []
        consecutive_empty = 0
        for line in initial_tex_file:
            # remove double space
            if line == '\n':
                consecutive_empty += 1
            else:
                consecutive_empty = 0
            if consecutive_empty <= 1:
                new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)

    def convert_itemize(self):
        initial_tex_file = self.import_tex_file()
        in_verbatim = False
        add_start = False
        add_end = False
        within_item = False
        new_tex_file_name = []
        consecutive_item = 0
        for line in initial_tex_file:
            in_verbatim = is_in_verbatim(in_verbatim, line)
            add_start = False
            add_end = False
            if in_verbatim is False:
                if (line[0] == '-'):
                    line = '\item' + line[1:]
                    if consecutive_item == 0:
                        add_start = True
                    consecutive_item += 1
                else:
                    if (consecutive_item > 0) & (line[0] != ' '):
                        add_end = True
                    consecutive_item = 0
                if (add_start) & (within_item is False):
                    new_tex_file_name.append(r'\begin{itemize}'+'\n')
                    add_start = False
                    within_item = True
                elif (add_end) & (within_item):
                    new_tex_file_name.append(r'\end{itemize}'+'\n')
                    add_end = False
                    within_item = False
                    consecutive_item = 0
            new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)

    def import_tex_file(self):
        f = open(self.file_name, "r") 
        initial_tex_file = []
        for line in f:
            initial_tex_file.append(line)
        f.close()
        return initial_tex_file

    def write_file(self, new_tex_file_name):
        # write filtered file        
        f = open(self.file_name, "w") 
        for line in new_tex_file_name:
            f.write(line)
        f.close()


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
        elif ("bw" in block_type):
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

class ReadRST:
    """Read RST file."""
    def __init__(self, file_name, *args, **kwargs,):
        """Initialize"""
        super().__init__(*args, **kwargs)
        self.file_name = file_name

    def convert_file(self):
        """Main convert function."""
        self.read_rst()
        self.detect_blocks()
        self.detect_sub_blocks()
        #self.detect_label()
        self.detect_title()

    def read_rst(self):
        """Convert the rst file into a list of strings"""
        rst = rstparse.Parser()
        with open(self.file_name) as f:
            rst.read(f)
        rst.parse()
        file_content = []
        for line in rst.lines:
            file_content.append(line)
        self.file_content = file_content

    def detect_blocks(self):
        """Detect the blocks in the file"""
        main_block_type, main_block = read_block(self.file_content)
        self.main_block_type = main_block_type
        self.main_block = main_block 
     
    def detect_sub_blocks(self):
        """Detect the sub-blocks in the file"""
        sub_block_type, sub_block = read_subblock(self.file_content)
        self.sub_block_type = sub_block_type
        self.sub_block = sub_block 

    def detect_label(self):
        """Detect the labels in the file"""
        positions, types, ids = read_block(self.file_content, ['-label'])
        self.label_ids = ids
        self.label_types = types
        self.label_positions = positions
            
    def detect_title(self):
            self.detect_title_position()
            assert np.sum(np.array(self.title_types) == "main") == 1, """More than one main title was found"""

    def detect_title_position(self):
        self.title_positions = []
        self.title_types = []
        for n, line in enumerate(self.file_content):
            if line[:3] == "***":
                self.title_positions.append(n-1)
                self.title_types.append("main")
            elif line[:3] == "===":
                self.title_positions.append(n-1)
                self.title_types.append("subtitle")
            elif line[:3] == "---":
                self.title_positions.append(n-1)
                self.title_types.append("subsubtitle")

def fix_caption(line):
    caption = None
    if ':caption:' in line:
        caption = line.split('*')[1]+':\n'
    return caption

def read_block(file_content):
    main_block = []
    main_block_type = []
    cpt_main_block = 0
    type = 'start'
    for line in file_content:
        new_block = False
        if ('.. ' in line) & (' .. ' not in line) & ('...' not in line):
            # new main block with no indentation
            cpt_main_block += 1
            new_block = True
        elif len(line) > 0:
            if line[0] != ' ':
                cpt_main_block += 1
                new_block = True
        if new_block:
            if ('container:: justify' in line) | ('container:: abstract' in line):
                type = 'text'
            elif 'container:: hatnote' in line: 
                 type = 'hatnote'
            elif 'admonition::' in line:
               type = 'admonition' + line.split('::')[1]
            elif ('code-block' in line) & ('bw' in line):
                type = 'bw-equation'
            elif ('code-block' in line) & ('lammps' in line):
                type = 'lammps-equation'
            elif ('figure:: ' in line):
                type = 'figure::' + line.split('::')[1]
            elif 'math::' in line:
                type = 'math'
            else:
                type = 'unknown'
                # print("unknown type", line)
        main_block_type.append(type)
        main_block.append(cpt_main_block)
    return main_block_type, main_block

def count_line(line):
    space_number = 0
    empty_start = True
    while empty_start:
        if len(line)>0:
            if line[0] == ' ':
                line = line[1:]
                space_number += 1
            else:
                empty_start = False
        else:
            empty_start = False
    return space_number

def read_subblock(file_content):
    sub_block = []
    sub_block_type = []
    cpt_sub_block = 0
    ref_space_number = 0
    type = 'start'
    for line in file_content:
        new_block = False
        space_number = count_line(line)
        if (' .. ' in line) & ('...' not in line):
            # new main block with no indentation
            cpt_sub_block += 1
            new_block = True
            ref_space_number = space_number
        elif (space_number <= ref_space_number) & (len(line) > 0):
            cpt_sub_block += 1
            ref_space_number = space_number
            new_block = True
            type = 'unknown'
        if new_block:
            if len(line) > 0:
                if line[0] != ' ':
                    type = 'unknown'
                else:
                    if 'container:: justify' in line:
                        type = 'text'
                    elif 'container:: hatnote' in line: 
                        type = 'hatnote'
                    elif 'admonition::' in line:
                        type = 'admonition'
                    elif ('code-block' in line) & ('bw' in line):
                        type = 'bw-equation'
                    elif ('code-block' in line) & ('lammps' in line):
                        type = 'lammps-equation'
                    elif 'math::' in line:
                        type = 'math'
                    elif ('figure:: ' in line):
                        type = 'figure'
                    else:
                        type = 'unknown'
                        # print("unknown type", line)
        sub_block_type.append(type)
        sub_block.append(cpt_sub_block)
    return sub_block_type, sub_block

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

def identify_subblock_id(n_subblock, block_lines, RST):
    if n_subblock > 0:
        
        sub_block = []
        for n in block_lines:
            type = RST.sub_block_type[n]
            id = RST.sub_block[n]
            if type != 'unknown':
                sub_block.append(id)
        id = np.unique(sub_block)
        assert len(id) == n_subblock

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
                    if ('..' in line) & ('...' not in line):
                        filtered_text.append('[insert-sub-block]')
                        n_subblock += 1
                    else:
                        filtered_text.append(line[first_non_zero_indentation:])
        return filtered_text, n_subblock
    else:
        return None, 0

def identify_subblock(block_lines, subblock_positions, subblock_types, subblock_ids):
    """Identify the subblock that are within a given block"""
    sub_type_in = []
    sub_id_in = []
    sub_pos_in = []
    for sub_pos in subblock_positions:
        if sub_pos in block_lines:
            sub_type = subblock_types[sub_pos]
            sub_id = subblock_ids[sub_pos]
            sub_type_in.append(sub_type)
            sub_id_in.append(sub_id)
            sub_pos_in.append(sub_pos)
    return sub_type_in, sub_id_in, sub_pos_in

def is_the_block_main(file_content, position):
    main_block = False
    first_line = file_content[position]
    cpt = 0
    while first_line[0] == ' ':
        cpt += 1
        first_line = first_line[1:]
    if cpt == 0:
        main_block = True                
    return main_block

def count_occurence(main_string, sub_string):
    count=0
    start_index=0
    positions = []
    for _ in range(len(main_string)):
        j = main_string.find(sub_string,start_index)
        if(j!=-1):
            start_index = j+1
            count+=1
            positions.append(j)
    return count, positions

def extract_link(RST, label):
    n = 0
    position_link = []
    for m, line in enumerate(RST.file_content):
        if label in line:
            n += 1
            position_link.append(m)
    assert n == 2
    position_link = position_link[1]
    for sub in RST.file_content[position_link+2].split():
        if 'href' in sub:
            link = sub[6:-1]
    
    for sub in RST.file_content[position_link+2].split('>'):
        if '</a' in sub:
            text = sub[:-3]
    return link, text

def fix_link(RST, line):
    """Deal with links"""
    if '|' in line:
        _, pos_bar = count_occurence(line, '|')
        links = []
        texts = []
        labels = []
        for ini, end in zip(pos_bar[::2], pos_bar[1::2]):
            label = line[ini:end+1]
            labels.append(label[1:-1])
            link, text = extract_link(RST, label)
            links.append(link)
            texts.append(text)
        l = 0
        sentence = line.split('|')
        new_line = ''
        for part in sentence:
            if part not in labels:
                new_line += part
                try:
                    replacement = '\href{'+links[l]+'}{'+texts[l]+'}'
                    new_line += replacement
                except:
                    pass                        
                l += 1
        return new_line
    else:
        return line

def fix_italic(line, caracter = '*', replace_with = [r'\textit{', '}'], replace_underscore=False):
    """Deal with special characters"""
    if caracter in line:
        _, pos_bar = count_occurence(line, caracter)
        italics = []
        for ini, end in zip(pos_bar[::2], pos_bar[1::2]):
            italics.append(line[ini:end+1][1:-1])
        l = 0
        sentence = line.split(caracter)
        new_line = ''
        for part in sentence:
            if part not in italics:
                new_line += part
                try:
                    italic = italics[l]
                    if replace_underscore:
                        italic = replace_special_character(italic, '_', '$\_$')
                    replacement = replace_with[0]+italic+replace_with[1]
                    new_line += replacement
                except:
                    pass                        
                l += 1
        return new_line
    else:
        return line
    
def fix_math(line):
    if ':math:' in line:
        _, pos_math = count_occurence(line, ':math:')
        _, pos_eq = count_occurence(line, "`")
        equations = []
        for ini, end in zip(pos_eq[::2], pos_eq[1::2]):
            equations.append(line[ini:end+1][1:-1])
        l = 0
        sentence = line.split("`")
        new_line = ''
        for part in sentence:
            if part not in equations:
                new_line += part
            else:
                try:
                    new_line += '$'+equations[l]+'$'             
                    l += 1
                except:
                    pass
        rests = new_line.split(':math:')
        new_line = ''
        for rest in rests:
            new_line += rest
        return new_line
    else:
        return line

def replace_special_character(line, to_replace, replace_with):
    if to_replace in line:
        sentence = line.split(to_replace)
        if len(line.split(to_replace)) == 2:
            new_line = sentence[0] + replace_with + sentence[1]
        else:
            print('ERROR ----', len(line.split(to_replace)))          
        return new_line
    else:
        return line

def is_in_verbatim(in_verbatim, line):
    if 'begin{lcverbatim}' in line:
        in_verbatim = True
    elif 'end{lcverbatim}' in line:
        in_verbatim = False
    return in_verbatim