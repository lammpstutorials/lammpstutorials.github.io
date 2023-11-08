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
        self.detect_label()
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
        positions, types, ids = read_block(self.file_content,
                                           ['container:: justify', 'container:: hatnote', 'admonition::'])
        self.block_ids = ids
        self.block_types = types 
        self.block_positions = positions
     
    def detect_sub_blocks(self):
        """Detect the sub-blocks in the file"""
        positions, types, ids = read_block(self.file_content,
                                           ['figure::', ':: lammps', ':: bw', ':: python'])
        self.subblock_ids = ids
        self.subblock_types = types
        self.subblock_positions = positions

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

def read_block(file_content, possible_types):
    positions = []
    types = []
    ids = []
    last_type = 'start'
    id = 0
    for n, line in enumerate(file_content):
        if ('..' in line) & ('...' not in line):
            part = line.split('..')
            for type in possible_types:
                if type in part[1]:
                    id += 1
                    last_type = type
                    positions.append(n)
        ids.append(id)
        types.append(last_type)
    return positions, types, ids

def clean_block(content, types, ids, id):
    line_numbers = np.where(np.array(ids) == id)[0]
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

def filter_block(block_text):
    """Convert a block of text into clean lines"""
    filtered_text = []
    keep_reading = True
    for line in block_text[1:]:
        if '..' in line:
            keep_reading = False
        if (keep_reading) & (len(line)>0):
            try:
                while line[0] == ' ':
                    line = line[1:]
                if ':class:' not in line:
                    filtered_text.append(line)
            except:
                pass
    return filtered_text

def identify_subblock(block_lines, subblock_positions, subblock_types, subblock_ids):
    """Identify the subblock that are within a given block"""
    sub_pos_in = []
    sub_id_in = []
    for sub_pos in subblock_positions:
        if sub_pos in block_lines:
            sub_type = subblock_types[sub_pos]
            sub_id = subblock_ids[sub_pos]
            sub_pos_in.append(sub_type)
            sub_id_in.append(sub_id)
    return sub_pos_in, sub_id_in

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
            print(len(line.split(to_replace)))          
            stop
        return new_line
    else:
        return line
