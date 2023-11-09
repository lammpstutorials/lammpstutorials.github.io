import rstparse
import numpy as np


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
            if 'container:: justify' in line:
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
                type = 'figure::'
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

##########################################






















def detect_block(n, file_content, keep_line_break=False):
    within_block = True
    words_in_block = []
    all_jump = []
    while within_block:
        try:
            line = file_content[n]
            jump = 0
            for letter in line:
                if letter == ' ':
                    jump += 1
                else:
                    break
            all_jump.append(jump)
            n += 1
            spitted_line = line.split()
            if line == '':
                words_in_block.append('\n')
            elif line[:2] == '  ':
                within_block = True
                for word in spitted_line:
                    words_in_block.append(word)
            else:
                within_block = False
            if keep_line_break:
                words_in_block.append('\n')
        except:
            within_block = False
    return words_in_block

def block_to_sentence(lines):
    sentences = []
    sentence = ''
    for word in lines:
        if word == '\n':
            if (sentence != '\n') & (sentence != ''):
                sentences.append(sentence[:-1])
            sentence = ''
        else:
            sentence += word + ' '
    if (sentence != '\n') & (sentence != ''):
        sentences.append(sentence[:-1])
    return sentences

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
        print(line)
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
        print(new_line)
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
            print(len(line.split(to_replace)))          
            stop
        return new_line
    else:
        return line

def is_in_verbatim(in_verbatim, line):
    if 'begin{lcverbatim}' in line:
        in_verbatim = True
    elif 'end{lcverbatim}' in line:
        in_verbatim = False
    return in_verbatim