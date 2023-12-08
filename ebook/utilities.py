import rstparse, os, shutil
import numpy as np
from PIL import Image

def is_in_verbatim(in_verbatim, line):
    """Search if line is in verbatim or not"""
    if 'begin{lcverbatim}' in line:
        in_verbatim = True
    elif 'end{lcverbatim}' in line:
        in_verbatim = False
    return in_verbatim

def count_line(line):
    """Count empty space at the beginning of a line"""
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

def replace_special_character(line, to_replace, replace_with):
    if to_replace in line:
        sentence = line.split(to_replace)
        if len(line.split(to_replace)) == 2:
            new_line = sentence[0] + replace_with + sentence[1]
        elif len(line.split(to_replace)) == 3:
            new_line = sentence[0] + replace_with + sentence[1] + replace_with + sentence[2]
        else:
            print(line)
            print('ERROR ----', len(line.split(to_replace)))          
        return new_line
    else:
        return line

def fix_math(line):
    """Fix equation"""
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

def read_label(sub_line):
    split = sub_line.split('`')
    rest = []
    label = None
    for sub in split:
        if '-label' in sub:
            label = sub
        else:
            rest.append(sub)
    return label, rest

def read_link(sub_line):
    cpt_link = 0
    cpt_rest = 0
    rest = ["", "", ""]
    link = ["", ""]
    in_link = False
    print(sub_line)
    for letter in sub_line:
        if (letter=='{') & (in_link==False):
            in_link = True
            cpt_rest += 1
        elif (letter=='}') & (in_link):
            in_link=False
            cpt_link += 1
        elif in_link:
            link[cpt_link] += letter
        elif in_link is False:
            rest[cpt_rest] += letter
    return rest, link






'''

    def fix_caption(line):
        caption = None
        if ':caption:' in line:
            caption = line.split('*')[1]+':\n'
        return caption

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



'''
